"""
At time of delivery for all RNASeq projects:

- Given input arguments of project, pi and recipe find all per run .bams for that project if it is RNASeq recipe
- Merge each .bam for the same sample when multiple .bams and write to the delivery/pipeline directory
- Log all input .bams merged in the delivery/pipeline directory
- Re-run setaccess.py (on a separate server)

At time of delivery for all 10X projects:
- Search under folder /igo/stats/CELLRANGER/ for any possible cell ranger output
- If existing, then copy to delivery/pipeline/cellranger directory
"""

from distutils.log import error
import os
import sys
import shutil
import logging
import glob
from subprocess import call
import requests
import re
import scripts.deliver_cellranger

LAB_SHARE_DIR = "/igo/delivery/share"
STATS_DIR = "/igo/staging/stats"
PICARD = "java -jar /igo/home/igo/resources/picard2.23.2/picard.jar "
NGS_STATS_FASTQ_ENDPOINT = "http://delphi.mskcc.org:8080/ngs-stats/permissions/getRequestPermissions/"

def deliver_pipeline_output(project, pi, recipe):
    if not project or not pi or not recipe:
        return "Project, pi and recipe are all required arguments."
    # change pi to all lowercase
    pi = pi.lower()
    delivery_folder = LAB_SHARE_DIR + "/" + pi + "/Project_" + project + "/pipeline"

    if recipe.startswith("RNASeq"):
        print("Delivering all RNASeq .bams for {} {} {}".format(project, pi, recipe))
        bamdict = find_bams(project, STATS_DIR)
        bsub_commands =  write_bams_to_share(bamdict, delivery_folder)
        reconcile_bam_fastq_list(project, bamdict)
        return bsub_commands
    
    # if 10X recipe or SCRI project starting with 12437, copy cell ranger result to project folder
    elif recipe.startswith("10XGenomics") or project.startswith("12437_"):
        folder_list = scripts.deliver_cellranger.find_cellranger(project)
        if len(folder_list) == 0:
            print("No cellragner result available")
        else:
            # create pipeline folder if not exists
            cellranger_delivery_folder = delivery_folder + "/cellranger"
            if not os.path.exists(cellranger_delivery_folder):
                print("Creating pipeline delivery folder {}".format(cellranger_delivery_folder))
                os.makedirs(cellranger_delivery_folder)

            # copy each sample folder to the delivery folder
            for folder in folder_list:
                sample_name = folder.split("/")[-1]
                sample_delivery_name = cellranger_delivery_folder + "/" + sample_name
                print("copy {}".format(folder))
                shutil.copytree(folder, sample_delivery_name)

    else:
        # TODO automate delivery of pipelines that are copied to the delivery share manually
        print("Pipeline delivery is not yet automated for recipe {} and project {}".format(recipe, project))
    return "Completed pipeline delivery"

def find_bams(project, stats_base_dir):
    """
    Find all bams for a project and return a dictionary of "igo_id" -> "bam list"
    """
    
    bam_unix_regex = stats_base_dir + '/**/*_IGO_' + project + '_*.bam'
    print("Searching for all .bams for project {} starting in folder {} matching glob {}".format(project, stats_base_dir, bam_unix_regex))
    # search for all .bams named like /igo/staging/stats/DIANA_0479_BHM2NVDSX3/RNA/GA28_ot_IGO_12785_H_1.bam
    project_bams = glob.glob(bam_unix_regex, recursive=True)  # recurisive=True causes the search to be a bit slow (more than 1 min)
    print("Total bams found {}".format(len(project_bams)))

    bamdict = {}
    for bam in project_bams:
        igo_id = bam.split("_IGO_")[1]
        igo_id = igo_id.replace(".bam","")
        print("Adding IGO ID {} bam {} to bam dictionary".format(igo_id, bam))
        if igo_id in bamdict: # add the bam to the list of bams for that igoId
            bamdict[igo_id].append(bam)
        else:
            bamdict.update({igo_id:[bam]}) # add dictionary entry with 1 list item
    
    return bamdict

def write_bams_to_share(bamdict, delivery_folder):
    """
    For each sample in the dictionary write the merged .bam and .log file to the pipeline folder.
    """
    print("Merging .bams if necessary and writing the output to {}".format(delivery_folder))
    if not os.path.exists(delivery_folder):
        print("Creating pipeline delivery folder {}".format(delivery_folder))
        os.makedirs(delivery_folder)
    
    bsub_commands = []
    # setup log file to record all .bams copied to the delivery folder
    log_file = delivery_folder + "/bamCreation.log"
    logging.basicConfig(filename=log_file)

    for igo_id in bamdict:
        bamlist = bamdict[igo_id]
        # TODO .bam names are either "DIANA_0479_BHM2NVDSX3___P12785_H___GA28_ot_IGO_12785_H_1.bam" or just "GA28_ot_IGO_12785_H_1.bam" ?
        # dest_filename = next(reversed(bamlist[0].split("___"))) # for DIANA_0479_BHM2NVDSX3___P12785_H___GA28_ot_IGO_12785_H_1.bam
        dest_filename = os.path.basename(bamlist[0])
        print("Writing delivery .bam {} to folder {}".format(dest_filename, delivery_folder))
        if len(bamlist) == 0: # actually never skip the merge even when there is just 1 .bam
            # Copy GA28_ot_IGO_12785_H_1.bam" to delivery folder "GA28_ot_IGO_12785_H_1.bam"
            shutil.copy(bamlist[0], delivery_folder)
            msg = "Copied {} to {}".format(bamlist[0], delivery_folder)
            print(msg)
            logging.info(msg)
        else:
            print("Merging .bams {}".format(bamlist))
            bsub_merge = "bsub -J merge_bam_files_to_deliver_" + igo_id + " -o merge_bam_files___" + igo_id + ".out -n 40 -M 8 "
            merge_bams = PICARD + " MergeSamFiles O=" + delivery_folder + "/" + dest_filename + " " + " ".join("I=" + i for i in bamlist)
            bsub_merge_bams = bsub_merge + merge_bams
            print(bsub_merge_bams)
            logging.info(bsub_merge_bams)
            call(bsub_merge_bams, shell = True)
            bsub_commands.append(bsub_merge_bams)
    
    return bsub_commands

def reconcile_bam_fastq_list(project, bam_dict):
    """
    Confirm there is a .bam for every fastq.gz file delivered.  
    """
    # Fastq naming like: "/igo/delivery/FASTQ/KIM_0682_BHJVG7BCX2/Project_08822_C/Sample_XPRO_0034_T_IGO_08822_C_1/XPRO_0034_T_IGO_08822_C_1_S43_L001_R1_001.fastq.gz",
    #   bam naming like: "/igo/staging/stats/DIANA_0479_BHM2NVDSX3/RNA/GA28_ot_IGO_12785_H_1.bam"
    fastq_list = get_request_fastqs(project)

    for fastq in fastq_list:
        fastq_name = os.path.basename(fastq)
        fastq_igo_id = get_igo_id(fastq_name)
        error_list = []
        if fastq_igo_id not in bam_dict.keys():
            error_list.append(fastq_name + " has no matching .bam")

        if len(error_list) > 0:
            raise Exception('\n'.join(error_list))

def get_igo_id(fastq_name):
    igo_id = fastq_name.split("_IGO_")[1]
    # 12958_B_1_S18_L001_R1_001.fastq.gz
    igo_id = re.sub("_S(\d)+_L(\d)+_R[1|2]_001.fastq.gz", '', igo_id)
    return igo_id

def get_request_fastqs(request):
    url = NGS_STATS_FASTQ_ENDPOINT + request
    print("Sending request {}".format(url))
    response = requests.get(url).json()
    # 'status': 500, 'error': 'Internal Server Error',
    if 'status' in response.keys() and response['status'] == 500:
        return None
    return response['fastqs']

#optionally invoke directly, for example:
#python deliver_pipeline.py 13097 abdelwao RNASeq-TruSeqPolyA
if __name__ == '__main__':
    project = sys.argv[1]
    pi = sys.argv[2]
    recipe = sys.argv[3]
    deliver_pipeline_output(project, pi, recipe)
