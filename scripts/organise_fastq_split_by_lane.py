import re
import sys
import glob 
import os
import linecache
from subprocess import call

"""
Immitate the bcl2fastq option which created a sub-directory per sample from the sample sheet since IGO has been delivering fastq.gz files that way for many years.
"""
def create_fastq_folders(run_demux_dir):
    if run_demux_dir.endswith("/"):
        run_demux_dir = run_demux_dir[:len(run_demux_dir) - 1]
    os.chdir(run_demux_dir)

    for project in glob.iglob("Project_*"):
        # print(project)
        project_dir = run_demux_dir + "/" +  project
        os.chdir(project_dir)
        for fastq in glob.iglob("*fastq.gz"):
            sample_header = "Sample_"
            sample_name = re.search("(.+?)_IGO_(.+?)_(\d+)", fastq)[0]
            print("Processing: " + sample_name)
            fastq_folder = sample_header + sample_name
            fastqs_folders = next(os.walk("."))[1]
            if fastq_folder not in fastqs_folders:
                os.mkdir(fastq_folder, 0o775)
                print(fastq_folder)
                all_fastqs_of_one_sample = sample_name + "_*" + ".fastq.gz"
                move_fastq_2_folder = "mv " + all_fastqs_of_one_sample + " " + fastq_folder
                print(move_fastq_2_folder)
                call(move_fastq_2_folder, shell = True)
                # Special behavior for 08822 projects, create merged R1 fastq.gz & merged R2 fastq.gz for bwamem2 consumption
                if "_IGO_08822" in fastq_folder and not "_RNA_IGO_" in fastq_folder:
                    ppg_dir = run_demux_dir +"_PPG/" + project + "/" + fastq_folder + "/"
                    print("Creating PPG fastq dirs: " + ppg_dir)
                    os.makedirs(ppg_dir)
                    cat_r1_cmd = "cat " + fastq_folder + "/*_R1_001.fastq.gz > " + ppg_dir + "/" + sample_name + "_S01_R1_001.fastq.gz"
                    print(cat_r1_cmd)
                    call(cat_r1_cmd, shell = True)
                    cat_r2_cmd = "cat " + fastq_folder + "/*_R2_001.fastq.gz > " + ppg_dir + "/" + sample_name + "_S01_R2_001.fastq.gz"
                    print(cat_r2_cmd)
                    call(cat_r2_cmd, shell = True)


# regex. regular expression
# ^ = start at beginning of string. $ = start at end of string
# findall gives the output in a list, so [0] gets the 0th element out of the list
# look below
# SAMPLE = re.findall(r'^(.+?)___',FILE)[0]

"""
Add Sample_ prefix for sample folder for 10X atac demux(bcl2fastq)
"""
def correct_sample_folder_name(run_demux_dir):
    
    if run_demux_dir.endswith("/"):
        run_demux_dir = run_demux_dir[:len(run_demux_dir) - 1]
    os.chdir(run_demux_dir)
    
    # loop through each project folder
    for project in glob.iglob("Project_*"):
        # print(project)
        project_dir = run_demux_dir + "/" +  project
        os.chdir(project_dir)
        sample_list = os.listdir(project_dir)
        for sample in sample_list:
            sample_updated = "Sample_" + sample
            cmd = "mv {} {}".format(sample, sample_updated)
            print(cmd)
            call(cmd, shell = True)


"""
Creating the "Sample_" directories for each fastq made the DRAGEN generated Reports/fastq_list.csv paths incomplete, fix those paths.
"""
def correct_fastq_list_csv(demux_reports_dir):
    # Read the fastq_list.csv file in the demux "Reports" directory
    filename = demux_reports_dir + "/fastq_list.csv"
    print("Reading " + filename)
    with open(filename, 'r+') as f:
        text = f.read()
        f.seek(0)
        updated = re.sub('/([a-zA-Z0-9_-]+)_IGO_([0-9]{5}(_[A-Z]+)?(_[0-9]+))',r'/Sample_\1_IGO_\2/\1_IGO_\2', text)
        f.write(updated)

# optional invoke directly to move fastq file to sample folder or change sample folder for bcl2fastq
# usage: python organise_fastq_split_by_lane.py [create/correct] [demux_dir_full_path]
# example: python organise_fastq_split_by_lane.py create /igo/staging/FASTQ/MICHELLE_0552_AHYF7NDSX3_v2
if __name__ == '__main__':
    demux_type = sys.argv[1]
    demux_dir = sys.argv[2]
    if demux_type == "create":
        create_fastq_folders(demux_dir)
        correct_fastq_list_csv(demux_dir+"/Reports")
    elif demux_type == "correct":
        correct_sample_folder_name(demux_dir)
    else:
        print("Demux type has to be either create(for bclconvert) or correct(for bcl2fastq)")
