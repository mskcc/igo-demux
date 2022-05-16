"""
At time of delivery for all RNASeq projects:

- Given input arguments of project, pi and recipe find all per run .bams for that project if it is RNASeq recipe
- Merge each .bam for the same sample when multiple .bams and write to the delivery/pipeline directory
- Log all input .bams merged in the delivery/pipeline directory
- Re-run setaccess.py (on a separate server)
"""

import os
import sys
import shutil
import logging
import glob
from subprocess import call

LAB_SHARE_DIR = "/igo/delivery/share"
STATS_DIR = "/igo/staging/stats"
PICARD = "java -jar /igo/home/igo/resources/picard2.23.2/picard.jar "

def deliver_pipeline_output(project, pi, recipe):
    delivery_folder = LAB_SHARE_DIR + "/" + pi + "/Project_" + project + "/pipeline"

    if recipe.startswith("RNASeq"):
        print("Delivering all RNASeq .bams for {} {} {}".format(project, pi, recipe))
        bamdict = find_bams(project, STATS_DIR)
        write_bams_to_share(bamdict, delivery_folder)
    else:
        # TODO automate delivery of pipelines that are copied to the delivery share manually
        print("Pipeline delivery is not yet automated for recipe {} and project {}".format(recipe, project))

def find_bams(project, stats_base_dir):
    """
    Find all bams for a project and return a dictionary of "igo_id" -> "bam list"
    """
    
    bam_unix_regex = stats_base_dir + '/**/*_IGO_' + project + '_*.bam'
    print("Searching for all .bams for project {} starting in folder {} matching glob {}".format(project, stats_base_dir, bam_unix_regex))
    # search for all .bams named like /igo/staging/stats/DIANA_0479_BHM2NVDSX3/RNA/DIANA_0479_BHM2NVDSX3___P12785_H___GA28_ot_IGO_12785_H_1.bam
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
        dest_filename = next(reversed(bamlist[0].split("___")))
        print("Writing delivery .bam {} to folder {}".format(dest_filename, delivery_folder))
        if len(bamlist) == 1: # skip merge if only one .bam
            # Copy and Rename "DIANA_0479_BHM2NVDSX3___P12785_H___GA28_ot_IGO_12785_H_1.bam" to just "GA28_ot_IGO_12785_H_1.bam"
            shutil.copy(bamlist[0], delivery_folder)
            print("Copied {} to {}".format(bamlist[0], delivery_folder))
            shutil.move(delivery_folder + "/" + os.path.basename(bamlist[0]), delivery_folder + "/" + dest_filename)
            msg = "Copied {} to pipeline delivery folder".format(bamlist[0])
            logging.info(msg)
            print(msg)
        else:
            print("Merging .bams {}".format(bamlist))
            bsub_merge = "bsub -J merge_bam_files_to_deliver_" + igo_id + " -o merge_bam_files___" + igo_id + ".out -n 40 -M 8 "
            merge_bams = PICARD + "MergeSamFiles --OUTPUT " + dest_filename + " " + " ".join("--INPUT " + i for i in bamlist)
            bsub_merge_bams = bsub_merge + merge_bams
            print(bsub_merge_bams)
            logging.info(bsub_merge_bams)
            call(bsub_merge_bams, shell = True)
            bsub_commands.append(bsub_merge_bams)
    
    return bsub_commands

#optionally invoke directly, for example:
#python deliver_pipeline.py 13097 abdelwao RNASeq-TruSeqPolyA
if __name__ == '__main__':
    project = sys.argv[1]
    pi = sys.argv[2]
    recipe = sys.argv[3]
    deliver_pipeline_output(project, pi, recipe)