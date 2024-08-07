import pandas as pd
import statistics
import sys
import glob
import os
from collections import OrderedDict

# TODO get barcode info from lims
# check if the run is pooled
def if_pooled(sequencing_summary_df):
    pooled = False
    if "barcode_kit" in sequencing_summary_df.columns:
        pooled = True
    return pooled

# get stats metric if the run is not pooled
def get_read_length_and_summary(sequencing_summary_df):
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
    return(len(read_length), N50_value * 2 / 1000000000, N50, median)

# get stats metric if the run is pooled
def get_read_length_and_summary_pooled(sequencing_summary_df, sample_name):
    sample_dict = {}
    samples = sequencing_summary_df["barcode_arrangement"].unique()
    for sample in samples:
        sample_df = sequencing_summary_df.loc[sequencing_summary_df['barcode_arrangement'] == sample]
        sample_sub = sample_name + "_" + sample
        stats = get_read_length_and_summary(sample_df)
        # only record barcodes with more than 10000 reads
        if stats[0] > 10000:
            sample_dict[sample_sub] = get_read_length_and_summary(sample_df)
    return sample_dict

def write_to_csv(sample_dict):
    file_name = "summary.csv"
    print("Writing stats file: " + file_name)
    with open(file_name,'w') as file:
        file.write("sample_id, Reads, Bases, N50, Meidan Read Length\n")
        for key, value in sample_dict.items():
            file.write("{}, {}, {}, {}, {}\n".format(key, value[0], value[1], value[2], value[3]))

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
                summary_matrix = pd.read_csv(i, delimiter = "\t")
                pooled = if_pooled(summary_matrix)
                # give different sample name for multi runs on one flow cell
                if file_count != 1:
                    sample = sample + "_" + str(file_count)
                if pooled:
                    sample_dict_sub = get_read_length_and_summary_pooled(summary_matrix, sample)
                    sample_dict.update(sample_dict_sub)
                else:
                    sample_dict[sample] = get_read_length_and_summary(summary_matrix)

    write_to_csv(sample_dict)
    print("ONT stats complete for: " + project_directory)
