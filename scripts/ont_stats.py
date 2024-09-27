import pandas as pd
import requests
import statistics
import sys
import glob
import os
from collections import OrderedDict
import re

# TODO get barcode info from lims
# check if the run is pooled
def if_pooled(sequencing_summary_df):
    pooled = False
    if "barcode_kit" in sequencing_summary_df.columns:
        pooled = True
    return pooled

# get stats metric if the run is not pooled
def get_read_length_and_summary(sequencing_summary_df, flowcell, position):
    read_length = sequencing_summary_df[sequencing_summary_df["passes_filtering"]]["sequence_length_template"].tolist()
    if len(read_length) != 0:
        read_length.sort(reverse = True)
        median = statistics.median(read_length)
        N50_value = sum(read_length) / 2
        total = 0
        for item in read_length:
            total += item
            if total >= N50_value:
                N50 = item
                break
    else:
        median = 0
        N50_value = 0
        N50 = 0
    return(len(read_length), N50_value * 2 / 1000000000, N50, median, flowcell, position)

# get stats metric if the run is pooled
def get_read_length_and_summary_pooled(sequencing_summary_df, sample_name, flowcell, position):
    sample_dict = {}
    samples = sequencing_summary_df["barcode_arrangement"].unique()
    for sample in samples:
        sample_df = sequencing_summary_df.loc[sequencing_summary_df['barcode_arrangement'] == sample]
        sample_sub = sample_name + "_" + sample
        stats = get_read_length_and_summary(sample_df, flowcell, position)
        # only record barcodes with more than 10000 reads
        if stats[0] > 10000:
            sample_dict[sample_sub] = get_read_length_and_summary(sample_df, flowcell, position)
    return sample_dict

def extract_flowcell(text):
    # Regular expression to match the characters after 'sequencing_summary_' and before the next '_'
    match = re.search(r'sequencing_summary_([^_]+)', text)
    if match:
        return match.group(1)
    else:
        return None

def write_to_csv(sample_dict, file_name):
    print("Writing stats file: " + file_name)
    with open(file_name,'w') as file:
        file.write("sample_id, Reads, Bases, N50, Median Read Length, Flowcell, Position\n")
        for key, value in sample_dict.items():
            file.write("{}, {}, {}, {}, {}, {}, {}\n".format(key, value[0], value[1], value[2], value[3], value[4], value[5]))

if __name__ == '__main__':
    # Usage: python ont_stats.py [project_directory]
    # example: python ont_stats.py /igo/staging/promethion/Project_14607
    
    project_directory = sys.argv[1]
    os.chdir(project_directory)
    sample_list = next(os.walk("."))[1]
    sample_dict = {}
    sample_list.sort()
    for sample in sample_list:
        print("Processing sample: " + sample)
        destination = project_directory + "/" + sample
        file = glob.glob(destination + "/*/sequencing_summary_*")
        if len(file) != 0:
            file_count = 0
            for i in file:
                file_count += 1
                position = i.split("/")[-2].split("_")[2]
                flowcell = extract_flowcell(i)
                print("Processing file: " + i + " from flowcell: " + flowcell + " at position:" + position)
                summary_matrix = pd.read_csv(i, delimiter = "\t")
                pooled = if_pooled(summary_matrix)
                # give different sample name for multi runs on one flow cell
                if file_count != 1:
                    sample = sample + "_" + str(file_count)
                if pooled:
                    sample_dict_sub = get_read_length_and_summary_pooled(summary_matrix, sample, flowcell, position)
                    sample_dict.update(sample_dict_sub)
                else:
                    sample_dict[sample] = get_read_length_and_summary(summary_matrix, flowcell, position)
                print(sample_dict)

    print(sample_dict)
    write_to_csv(sample_dict, "summary.csv")
    print("ONT stats .csv complete for: " + project_directory)
    
    # List of parameter names corresponding to the values skipping columns "estimatedCoverage", "bamCoverage", "sequencerName"
    parameter_names = ["reads", "bases", "N50", "medianReadLength", "flowcell", "sequencerPosition", "igoId"]

    # Convert initial dictionary to a nested dictionary with parameter names
    converted_sample_dict = {}
    for key, values in sample_dict.items():
        # Create a nested dictionary by zipping parameter names and values
        values = values + (key,)
        converted_sample_dict[key] = dict(zip(parameter_names, values))
    print(converted_sample_dict)

    # Write to LIMS endpoint with a GET:
    # /LimsRest/updateLimsSampleLevelSequencingQcONT?igoId=04540_U_26_1_1_1_1_1&flowcell=PAY61078&reads=19775442&bases=9103016668&N50=16508&medianReadLength=766&estimatedCoverage=0&bamCoverage=0&sequencerPosition=1A&sequencerName=zeppelin
    LIMS_ENDPOINT="https://igo-lims02.mskcc.org:8443/LimsRest/updateLimsSampleLevelSequencingQcONT"
    for sample_id, params in converted_sample_dict.items():
        # Send GET request for each set of parameters
        print("Sending LIMS get request for: " + params)
        response = requests.get(LIMS_ENDPOINT, params=params, verify=False)

        # Check the response status and print the output
        if response.status_code == 200:
            print(f"Request for {sample_id} successful!")
            print("Response Data:", response.json())
        else:
            print(f"Request for {sample_id} failed with status code {response.status_code}")
            print("Error details:", response.text)