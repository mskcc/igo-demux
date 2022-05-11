"""
At time of delivery for all RNASeq projects:

- Given input arguments of project, pi and recipe find all per run .bams for that project if it is RNASeq recipe
- Merge each .bam for the same sample when multiple .bams and write to the delivery/pipeline directory
- Log all input .bams merged in the delivery/pipeline directory
- Re-run setaccess.py (on a separate server)
"""

import os
import shutil
import logging

LAB_SHARE_DIR = "/igo/delivery/share"
STATS_DIR = "/igo/staging/stats"

def deliver_pipeline_output(project, pi, recipe):
    delivery_folder = LAB_SHARE_DIR + "/" + pi + "/Project_" + project + "/pipeline"

    if recipe.startswith("RNASeq"):
        print("Delivering all RNASeq .bams for {} {}".format(project, recipe))
        bamdict = find_bams(project)
        write_bams_to_share(bamdict, delivery_folder)
    else:
        # TODO automate delivery of pipelines that are copied to the delivery share manually
        print("Pipeline delivery is not yet automated for recipe {} and project {}".format(recipe, project))

def find_bams(project):
    """
    Find all bams for a project and return a dictionary of "igo_id" -> "bam list"
    """
    # search for all .bams named like DIANA_0479_BHM2NVDSX3___P12785_H___GA28_ot_IGO_12785_H_1___GRCh38.bam
    print("Searching for all .bams for project {}".format(project))
    bamdict = {}
    return bamdict

def write_bams_to_share(bamdict, delivery_folder):
    """
    For each sample in the dictionary write the merged .bam and .log file to the pipeline folder.
    """
    print("Merging .bams if necessary and writing the output to {}".format(delivery_folder))
    if not os.path.exists(delivery_folder):
        print("Creating pipeline delivery folder {}".format(delivery_folder))
        os.makedirs(delivery_folder)
    
    # setup log file to record all .bams copied to the delivery folder
    log_file = delivery_folder + "/bamCreation.log"
    logging.basicConfig(filename=log_file)

    for igo_id in bamdict:
        bamlist = bamdict[igo_id]
        if len(bamlist) == 1: # skip merge if only one .bam
            # TODO Rename "DIANA_0479_BHM2NVDSX3___P12785_H___GA28_ot_IGO_12785_H_1___GRCh38.bam" to just "GA28_ot_IGO_12785_H_1.bam"
            shutil.copy(bamlist[0], delivery_folder)
            logging.info("Copied {} to pipeline delivery folder".format(bamlist[0]))
        else:
            print("Merging .bams {}".format(bamlist))
            #TODO Picard MergeSamFiles https://github.com/mskcc/igo-demux/blob/bf27d2a69d0ea640074ce4081dc035fe7178d9e6/scripts/alignment_and_picard.py#L268