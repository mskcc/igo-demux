import re
import sys
import glob	
import os
import linecache
from subprocess import call

"""
Immitate the bcl2fastq option which created a sub-directory per sample from the sample sheet.
"""

def create_fastq_folders(run_demux_dir):
    run = run_demux_dir
    
    os.chdir(run)

    x = 0
    for project in glob.iglob("Project_*"):
        # print(project)
        project_dir = run + "/" +  project
        os.chdir(project_dir)
        for fastq in glob.iglob("*fastq.gz"):
            sample_header = "Sample_"
            sample_name = re.search("(.+?)_IGO_(.+?)_(\d+)", fastq)[0]
            fastq_folder = sample_header + sample_name
            fastqs_folders = next(os.walk("."))[1]
            if fastq_folder not in fastqs_folders:
                os.mkdir(fastq_folder, 0o775)
                print(fastq_folder)
                all_fastqs_of_one_sample = sample_name + "*" + ".fastq.gz"
                move_fastq_2_folder = "mv " + all_fastqs_of_one_sample + " " + fastq_folder
                print(move_fastq_2_folder)
                x += 1
                print(x)
                call(move_fastq_2_folder, shell = True)

# regex. regular expression
# ^ = start at beginning of string. $ = start at end of string
# findall gives the output in a list, so [0] gets the 0th element out of the list
# look below
# SAMPLE = re.findall(r'^(.+?)___',FILE)[0]

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
        updated = re.sub('/([a-zA-Z0-9_]+)_IGO_([0-9]{5}(_[A-Z]+)?(_[0-9]+))',r'/Sample_\1_IGO_\2/\1_IGO_\2', text)
        f.write(updated)