# function to run dragen pipeline for tumor/normal pair of GLP samples. 
# input of igo_id and get tumor/normal list of the same patient then generate dragen commnad
# output should be under /igo/staging/GLP_pipeline/ followed by patient id folder then normal_tumor_igo_id.

import os
import requests
import json
import pandas as pd
import sys


OUTPUT_PATH = "/igo/staging/GLP_pipeline/"
COMMAND_PREFIX = "/opt/edico/bin/dragen -f -r /igo/work/igo/dragen_hash_tables/hg38-alt_masked.cnv.graph.hla.rna-9-r3.0-1"
COMMAND_PARAMETER = "--enable-duplicate-marking true --enable-map-align true --enable-sv true --enable-map-align-output true --output-format bam --enable-variant-caller true --enable-hla true"
FASTQ_PATH = "/igo/staging/FASTQ/"
ENDPOINT = "https://igolims.mskcc.org:8443/LimsRest/getSamplePairs?igoSampleId="
FASTQ_LIST_HEADER = "RGID,RGSM,RGLB,Lane,Read1File,Read2File"

# find fastq file under staging drive by list of fastq sample folder name eg: [GLP-CCV-0001-WB_IGO_16606_B_1]
def find_fastq_file(sample_ID_list):
    # get whole list of all fastq files that available with project folder as tag
    run_list = os.listdir(FASTQ_PATH)
    # dictionary of run_ID->project_list
    run_project_dict = {}
    for run_ID in run_list:
        current_path = FASTQ_PATH + run_ID + "/"
        if os.path.isdir(current_path):
            run_project_dict[run_ID] = []
            file_list = os.listdir(current_path)
            for item in file_list:
                if "Project_" in item:
                    run_project_dict[run_ID].append(item)
    
    fastq_file_list_dict = {}
    # get fastq file path list for given sample_ID
    for sample_ID in sample_ID_list:
        fastq_file_list = []
        project_ID = "Project_" + "_".join(sample_ID.split("_")[sample_ID.split("_").index("IGO") + 1:-1])
        sample_folder_name = "Sample_" + sample_ID
        for run_ID, project_list in run_project_dict.items():
            if project_ID in project_list:
                fastq_file_path_prefix = FASTQ_PATH + run_ID + "/" + project_ID + "/"
                sample_folder_list = os.listdir(fastq_file_path_prefix)
                if sample_folder_name in sample_folder_list:
                    fastq_file_list.append(fastq_file_path_prefix + sample_folder_name)
        fastq_file_list_dict[sample_ID] = fastq_file_list
    return fastq_file_list_dict

def gather_sample_list(IGO_id):
    response = requests.get(ENDPOINT + IGO_id , auth = ("pms", "tiagostarbuckslightbike"), verify = False)
    if response.status_code == 200:
        response_data = json.loads(response.text.encode("utf8"))
        patient_id = response_data["patientId"]
        normal_lst = response_data["normalIgoIdSampleName"]
        tumor_lst = response_data["tumorIgoIdSampleName"]
        patient_dict = {patient_id: {"normal": normal_lst, "tumor": tumor_lst}}
    else:
        print("endpoint failed for sample {}, please check".format(IGO_id))
        return
    return patient_dict

def create_fastq_list(fastq_file_list_dict, output_file):
    # get the line from previous fastq file list file
    merged_data = []
    for rgsm_keyword, value in fastq_file_list_dict.items():
        # Load the CSV file
        for file_path in value:
            csv_file = "/".join(file_path.split("/")[0:5]) + "/Reports/fastq_list.csv"
            print("gathering info from {}".format(csv_file))
            data = pd.read_csv(csv_file)
            # Filter rows based on the RGSM keyword
            filtered_data = data[data['RGSM'] == rgsm_keyword]
            # Append the filtered data to the list
            merged_data.append(filtered_data)

    # Concatenate all DataFrames
    merged_result = pd.concat(merged_data, ignore_index=True)
    # Save the merged results with the header line to a file
    with open(output_file, "w") as f:
        f.write(FASTQ_LIST_HEADER + "\n") 
        merged_result.to_csv(f, index=False, header=False)  

    return merged_result

# generate cmd and create folder for output
def generate_cmd(fastq_file_list, patient_id, tumor_name, normal_name):
    tumor_igo_id = "_".join(tumor_name.split("_")[tumor_name.split("_").index("IGO") + 1:])
    normal_igo_id = "_".join(normal_name.split("_")[normal_name.split("_").index("IGO") + 1:])
    sample_output_path = "{}{}/{}_{}".format(OUTPUT_PATH, patient_id, normal_igo_id, tumor_igo_id)
    if not os.path.exists(sample_output_path):
        os.makedirs(sample_output_path)
    job_name = "{}_{}_GLP_pipeline".format(tumor_name, normal_name)
    cmd = "bsub -J {} -o {}.out -n48 -q dragen {} --tumor-fastq-list {} --fastq-list {} --tumor-fastq-list-sample-id {} --fastq-list-sample-id {} --output-directory {} --output-file-prefix {} {}".format(job_name, job_name, COMMAND_PREFIX, fastq_file_list, fastq_file_list, tumor_name, normal_name, sample_output_path, patient_id, COMMAND_PARAMETER)
    return cmd

def main(igo_id):
    # gather sample list and corresponding fastq list
    patient_dict = gather_sample_list(igo_id)
    print(patient_dict)
    fastq_file_dict = {}
    for key, value in patient_dict.items():
        # only start process when tumor and normal are all exist
        if len(value["normal"]) != 0 and len(value["tumor"]) != 0:
            # find fastq files for all the samples with same patient id
            fastq_file_dict.update(find_fastq_file(value["normal"]))
            fastq_file_dict.update(find_fastq_file(value["tumor"]))
            # create a fastq list file for all samples
            output_file_folder = "{}{}".format(OUTPUT_PATH, key)
            if not os.path.exists(output_file_folder):
                os.makedirs(output_file_folder)
            output_file = "{}/fastq_list.csv".format(output_file_folder)
            create_fastq_list(fastq_file_dict, output_file)
            for normal in value["normal"]:
                for tumor in value["tumor"]:
                    cmd = generate_cmd(output_file, key, tumor, normal)
                    print(cmd)

    return

if __name__ == '__main__':
    # launch GLP WES dragen commands by igo_id
    # Usage: python dragen_GLP.py [IGO_id]
    # example: python dragen_GLP.py 16606_B_1
    igo_id = sys.argv[1]
    main(igo_id)