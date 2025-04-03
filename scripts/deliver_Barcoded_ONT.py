# For barcoded ONT samples, make new sub sample folder and mv the corresponding barcoded fastq.
import os
import requests
import json
import re
import sys
import subprocess

ACCESS = 0o775

ONT_path = "/igo/staging/promethion/"

def get_numbers_from_string(input_string):
    # Use regular expression to find numbers at the end of the string
    match = re.search(r'\d+$', input_string)
    # Return the matched numbers or None if no match is found
    return match.group() if match else None

def get_project_ID_from_sample(sample_name):
    match = re.match(r"^(\d+)(_([A-Za-z]))?", sample_name)
    if match:
        if match.group(3):  # If the second part is a letter
            return f"{match.group(1)}_{match.group(3)}"
        else:  # Otherwise, return only the first part
            return match.group(1)
    return None

# get barcode information from lims
def get_sample_barcode(pool_id):
    sample_dict = {}
    lims_endpoint = "https://igolims.mskcc.org:8443/LimsRest/getPoolsBarcodes?poolId="
    response = requests.get(lims_endpoint + pool_id, auth = ("pms", "tiagostarbuckslightbike"), verify = False)
    response_data = json.loads(response.text.encode("utf8"))
    for sample in response_data:
        project_id = "Project_" + str(get_project_ID_from_sample(sample["librarySample"]))
        sample_dict[sample["librarySample"]] =[project_id ,"barcode" + str(get_numbers_from_string(sample["sampleBarcode"]["barcodId"]))]
    print(sample_dict)
    return sample_dict

# given sample:[project,barcode] dictionary per pool and pool path, create sample folder if not exist and generate cmd for mv and rename folder
def mv_fastq(sample_dict, parent_folder_name, pool_name):
    for key, value in sample_dict.items():
        project = value[0]
        barcode = value[1]
        # create sub sample folder if not exist
        sample_folder = "/igo/staging/promethion/{}/{}".format(project, key)
        destination_fastq_pass = "{}/fastq_pass".format(sample_folder)
        destination_fastq_fail = "{}/fastq_fail".format(sample_folder)
        destination_pod5_pass = "{}/pod5_pass".format(sample_folder)
        destination_pod5_fail = "{}/pod5_fail".format(sample_folder)

        if not os.path.exists(sample_folder):
            print("creating folder: " + sample_folder)
            os.makedirs(sample_folder, ACCESS)
            os.makedirs(destination_fastq_pass, ACCESS)
            os.makedirs(destination_fastq_fail, ACCESS)
            os.makedirs(destination_pod5_pass, ACCESS)
            os.makedirs(destination_pod5_fail, ACCESS)
        
        source_fastq_pass = "{}/{}/*/fastq_pass/{}".format(parent_folder_name, pool_name, barcode)
        source_fastq_fail = "{}/{}/*/fastq_fail/{}".format(parent_folder_name, pool_name, barcode)
        source_pod5_pass = "{}/{}/*/pod5_pass/{}".format(parent_folder_name, pool_name, barcode)
        source_pod5_fail = "{}/{}/*/pod5_fail/{}".format(parent_folder_name, pool_name, barcode)

        cmd1 = "mv {}/* {}/".format(source_fastq_pass, destination_fastq_pass)
        cmd2 = "mv {}/* {}/".format(source_fastq_fail, destination_fastq_fail)
        cmd3 = "mv {}/* {}/".format(source_pod5_pass, destination_pod5_pass)
        cmd4 = "mv {}/* {}/".format(source_pod5_fail, destination_pod5_fail)

        print(cmd1)
        print(cmd2)
        print(cmd3)
        print(cmd4)
        subprocess.run(cmd1, shell=True)
        subprocess.run(cmd2, shell=True)
        subprocess.run(cmd3, shell=True)
        subprocess.run(cmd4, shell=True)

if __name__ == '__main__':
    # Usage: python deliver_Barcoded_ONT.py [project_directory]
    # example: python deliver_Barcoded_ONT.py Project_15710_D_16975
    project_name = sys.argv[1]
    parent_folder_path = ONT_path + project_name
    for i in os.listdir(parent_folder_path):
        sample_dict = get_sample_barcode(i)
        mv_fastq(sample_dict, parent_folder_path, i)
