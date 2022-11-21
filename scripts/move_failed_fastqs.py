"""
When a fast is marked as "Failed" on the QC website:
    - find and move failed fastq in /igo/staging/FASTQ/{run} folder
    - move .bam in the /igo/staging/stats folder
"""

import sys
import os
import glob
import shutil
import re

FASTQ_DIR = "/igo/staging/FASTQ/"
FAILED_FASTQ_DIR = "/igo/staging/failed_fastq/"
STATS_DIR = "/igo/staging/stats/"

def move_failed_fastqs(igo_id, run):
    if not igo_id or not run:
        return "igo_id and run are required arguments."

    print("Removing failed fastqs and .bams for " + igo_id + " on run " + run)
    project_id = get_project_id(igo_id)
    igo_id = get_base_igo_id(igo_id)

    error_msg = "" # store error messages at any step

    # first try to move any fastq.gz files
    # for example: "/igo/staging/FASTQ/MICHELLE_0523_BHM52WDSX3/Project_13220_B"
    fastq_path = FASTQ_DIR + run + "/Project_" + project_id
    fastq_dest = FAILED_FASTQ_DIR + run + "/Project_" + project_id
    if not os.path.exists(fastq_path):
        error_msg = "Path does not exist:" + fastq_path
        print(error_msg)
    else:
        if not os.path.exists(fastq_dest):
            os.makedirs(fastq_dest)
        for folder in glob.glob(fastq_path + "/*_IGO_" + igo_id, recursive=False):
            # if the run is already archived this will mail run as the 'igo' user
            print("Moving folder: " + folder)
            dest = shutil.move(folder, fastq_dest)
            print("Moved to " + dest)

    # then remove any .bam files for the fastq.gz(s)   
    bam_path = STATS_DIR + run
    if not os.path.exists(bam_path):
        msg = "Bam path does not exist:" + bam_path
        print(msg)
    else:
        for file in glob.glob(bam_path + "/*_IGO_" + igo_id, recursive=True):
            print("Moving .bam " + file)
            shutil.move(file, fastq_dest)

    if error_msg != "":
        raise Exception(error_msg)
    else:
        return "Completed move_failed_fastqs"

def get_project_id(igo_id):
    # 13220_B_31 to project ID only
    return re.split(r"_(\d)+", igo_id)[0]

# defined for valid IGO IDs as input
def get_base_igo_id(igo_id):
    parts = re.split(r"_(\d)+", igo_id)
    return parts[0] +"_" + parts[1]

#optionally invoke directly, for example:
#python move_failed_fastqs.py igo_id run
if __name__ == '__main__':
    igo_id = sys.argv[1]
    run = sys.argv[2]
    move_failed_fastqs(igo_id, run)