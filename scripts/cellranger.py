# launch cell ranger pipeline (GE, VDJ, ATAC....) for 10X samples by recipe
# put result in /igo/stats/CELLRANGER/<run_ID>
# TODO make this function callable by run ID/project?

import pandas as pd
import re
import sys
import os
import glob
import json
import subprocess
from os.path import join
from os.path import basename
from os.path import abspath
from os.path import isdir
from subprocess import call

"""
input: sample_sheet object(for sample list and essential info), sequencer_and_run(for stats folder and fastq file location)
output: running cmd for cellranger by sample
"""

# work folder
STATS_AREA = '/igo/stats/CELLRANGER/'

# config info 
ACCESS = 0o775
config_dict = {
    "count" : {
        "tool" : ' /igo/work/nabors/tools/cellranger-6.1.2/cellranger count ',
        "genome" : {
            "Human" : ' --transcriptome=/igo/work/nabors/genomes/10X_Genomics/GEX/refdata-gex-GRCh38-2020-A ',
            "Mouse" : ' --transcriptome=/igo/work/nabors/genomes/10X_Genomics/GEX/refdata-gex-mm10-2020-A '
        }
    },
    "vdj" : {
        "tool" : ' /igo/work/nabors/tools/cellranger-6.1.2/cellranger vdj ',
        "genome" : {
            "Human" : ' --reference=/igo/work/nabors/genomes/10X_Genomics/VDJ/refdata-cellranger-vdj-GRCh38-alts-ensembl-2.0.0 ',
            "Mouse" : ' --reference=/igo/work/nabors/genomes/10X_Genomics/VDJ/refdata-cellranger-vdj-GRCm38-alts-ensembl-2.2.0 '
        }
    },
    "atac_count" : {
        "tool" : ' /igo/work/nabors/tools/cellranger-atac-2.1.0/cellranger-atac count ',
        "genome" : {
            "Human" : ' --reference=/igo/work/nabors/genomes/10X_Genomics/ATAC/refdata-cellranger-atac-GRCh38-1.0.1 ',
            "Mouse" : ' --reference=/igo/work/nabors/genomes/10X_Genomics/ATAC/refdata-cellranger-atac-mm10-1.1.0 '
        }
    },
    "cnv" : {
        "tool" : ' /igo/work/nabors/tools/cellranger-dna-1.1.0/cellranger-dna cnv ',
        "genome" : {
            "Human" : ' --reference=/igo/work/nabors/10X_Genomics/CNV/refdata-GRCh38-1.0.0 ',
            "Mouse" : ' --reference=/igo/work/nabors/10X_Genomics/CNV/refdata-GRCm38-1.0.0 '
        }
    },
    "multi" : {
        "tool" : ' /igo/work/nabors/tools/cellranger-6.1.2/cellranger multi '
    }
}

# cellranger command line options
OPTIONS = ' --nopreflight --jobmode=lsf --mempercore=64 --disable-ui --maxjobs=200'

# 10X recipe list for different pipelines TODO 10X_Genomics_Visium, 10X_Genomics_Multiome
COUNT_FLAVORS = ['10X_Genomics_GeneExpression-3', '10X_Genomics_GeneExpression-5']
VDJ_FLAVORS = ['10X_Genomics_VDJ']
ATAC_FLAVORS = ['10X_Genomics_ATAC']
CNV_FLAVORS = ['10X_Genomics_CNV']

"""
steps:
1. check whether there is previous fastq existing under /igo/staging/FASTQ (find_fastq_file)
2. get tag by recipe, if recipe not in the list above, skip for now (get_tag)
3. generate corresponding commands based on tag and genome (generate_cellranger_cmd)
4. run stats and create josn file for each project and call endpoint to push to qc website (launch_cellranger)
"""
# return a dictionary sample_ID -> list of fastq file path by given sample name list
def find_fastq_file(sample_ID_list):
    # get whole list of all fastq files that available with project folder as tag
    path_prefix = "/igo/staging/FASTQ/"
    run_list = os.listdir(path_prefix)
    # dictionary of run_ID->project_list
    run_project_dict = {}
    for run_ID in run_list:
        current_path = path_prefix + run_ID + "/"
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
                fastq_file_path_prefix = path_prefix + run_ID + "/" + project_ID + "/"
                sample_folder_list = os.listdir(fastq_file_path_prefix)
                if sample_folder_name in sample_folder_list:
                    fastq_file_list.append(fastq_file_path_prefix + sample_folder_name)
        fastq_file_list_dict[sample_ID] = fastq_file_list
    return fastq_file_list_dict

def get_tag(recipe):
    tag = "Skip"
    if recipe in COUNT_FLAVORS:
        tag = "count"
    if recipe in CNV_FLAVORS:
        tag = "cnv"    
    if recipe in VDJ_FLAVORS:
        tag = "vdj"
    if recipe in ATAC_FLAVORS:
        tag = "atac_count"
    return tag

# return tag and genome according to sample_ID for SCRI samples, all SCRI samples are starting with Project_12437
# eg: SD-1680_Patient_D_nucseq_H_VDJ_IGO_12437_AN_5 will given tag as vdj, genome as Human
# eg: SDtest_IGO_12437_AN_4 will given tag as skip, genome as na
# _H: Human, _M: Mouse
# _VDJ: vdj, _GE: count, _ATAC: "atac_count"
def get_SCRI_tag(sample_ID):
    tag_orig = sample_ID.split("_")[sample_ID.split("_").index("IGO") - 1]
    tag = "Skip"
    if tag_orig == "VDJ":
        tag = "vdj"
    if tag_orig == "GE":
        tag = "count"
    if tag_orig == "ATAC":
        tag = "atac_count"
    
    genome = "na"
    if tag != "Skip":
        genome_orig = sample_ID.split("_")[sample_ID.split("_").index("IGO") - 2]
        if genome_orig == "H":
            genome = "Human"
        if genome_orig == "M":
            genome = "Mouse"
    # if genome parameter couldn't detected, set tag back to skip
    if genome == "na":
        tag = "Skip"

    return tag, genome

def generate_cellranger_cmd(sample_ID, tag, genome, fastq_file_path, sequencer_and_run):
    tool = config_dict[tag]["tool"]
    transcriptome = config_dict[tag]["genome"][genome]
    project_ID = "Project_" + "_".join(sample_ID.split("_")[sample_ID.split("_").index("IGO") + 1:-1])
    cellranger_cmd = "{}--id=Sample_{}__{}".format(tool, sample_ID, tag) + transcriptome + '--fastqs=' + ','.join(fastq_file_path) + OPTIONS
    job_name = "{}_{}_{}_{}_cellranger".format(sequencer_and_run, project_ID, sample_ID, tag)
    bsub_cmd = "bsub -J {} -o {}.out{}".format(job_name, job_name, cellranger_cmd) 
    return bsub_cmd
        
def create_json(send_json, sequencer_and_run, project, tag, work_area):  
    job_id = sequencer_and_run + "_" + project            
    json_data_file = 'cellranger_json___' + sequencer_and_run + "__" + project + '.json'
    with open(json_data_file, 'w') as jfile:
        json.dump(send_json, jfile)
        
    bsub_json = 'bsub -J create_json___{} -o create_json___{}.log -w \"done({}*)\" sh /home/igo/Scripts/PicardScripts/send_json_data.sh {} {}'.format(job_id, job_id, job_id, work_area, json_data_file)
    print(bsub_json)
    subprocess.run(bsub_json, shell = True)

# Main function: launch cellranger cmd by given samplesheet object and sequencer_and_run
def launch_cellranger(sample_sheet, sequencer_and_run):
    # get parameters from sample_sheet
    # dictionary of Sample_ID->Project
    sample_project_dict = pd.Series(sample_sheet.df_ss_data['Sample_Project'].values,index=sample_sheet.df_ss_data['Sample_ID']).to_dict()
    # dictionary of project->sample_ID
    project_sample_dict = {}
    for sample_ID, project_ID in sample_project_dict.items():
        if project_ID in project_sample_dict.keys():
            project_sample_dict[project_ID].append(sample_ID)
        else:
            project_sample_dict[project_ID] = [sample_ID]
    # dictionary of sample_ID->recipe
    sample_recipe_dict = pd.Series(sample_sheet.df_ss_data['Sample_Well'].values,index=sample_sheet.df_ss_data['Sample_ID']).to_dict()
    # dictionary of sample_ID->genome
    sample_genome_dict = pd.Series(sample_sheet.df_ss_data['Sample_Plate'].values,index=sample_sheet.df_ss_data['Sample_ID']).to_dict()
    # dictionary of sample_ID->fastq_list
    sample_ID_list = list(sample_project_dict.keys())
    sample_fastqfile_dict = find_fastq_file(sample_ID_list)

    for project in project_sample_dict.keys():
        send_json = {}
        send_json['samples'] = []
        # CREATE RUN FOLDER AND PROJECT FOLDER IF NOT ALREADY THERE
        os.chdir(STATS_AREA)
        runs = next(os.walk('.'))[1]
        if sequencer_and_run not in runs:
            os.mkdir(sequencer_and_run, ACCESS)
                    
        stats_and_run = STATS_AREA + sequencer_and_run
        os.chdir(stats_and_run)
        projects = next(os.walk('.'))[1]
        if project not in projects:
            os.mkdir(project, ACCESS)
        work_area = stats_and_run + "/" + project + '/' 
        # GO TO project ID LOCATION to start cellranger command
        os.chdir(work_area)

        # SCRI samples don't need to be pushed onto qc website
        if "Project_12437" not in project:
            sample_list = project_sample_dict[project]
            # call cellranger for each sample and append info to json dict
            for sample in sample_list:
                tag = get_tag(sample_recipe_dict[sample])
                # if recipe within the tool being set up, lanuch cellranger
                if tag != "Skip":
                    if sample_genome_dict[sample] != "Human" and sample_genome_dict[sample] != "Mouse":
                        sample_genome_dict[sample] = "Mouse"
                    cmd = generate_cellranger_cmd(sample, tag, sample_genome_dict[sample], sample_fastqfile_dict[sample], sequencer_and_run)
                    print(cmd)
                    subprocess.run(cmd, shell=True)
                    send_json['samples'].append({'sample':'Sample_' + sample, 'type':tag, 'project':project, 'run':sequencer_and_run})
            if send_json['samples']:
                create_json(send_json, sequencer_and_run, project, tag, work_area)
        else:
            sample_list = project_sample_dict[project]
            # call cellranger for each sample
            for sample in sample_list:
                tag_genome = get_SCRI_tag(sample)
                tag = tag_genome[0]
                genome = tag_genome[1]
                # if recipe within the tool being set up, lanuch cellranger
                if tag != "Skip" and genome != "na":
                    cmd = generate_cellranger_cmd(sample, tag, genome, sample_fastqfile_dict[sample], sequencer_and_run)
                    cmd = cmd + " --include-introns=true"  # SCRI samples always have include-introns true
                    print(cmd)
                    subprocess.run(cmd, shell=True)

# sample_ID_list = ["06265_8869_1_IGO_06265_AG_3","Third-Transcriptome_IGO_11969_E_3", "Second_IGO_11969_E_2"]
# fastq_file_list_dict = {'06265_8869_1_IGO_06265_AG_3': ['/igo/staging/FASTQ/DIANA_0453_AHFKJ5DRXY/Project_06265_AG/Sample_06265_8869_1_IGO_06265_AG_3'], 'Third-Transcriptome_IGO_11969_E_3': ['/igo/staging/FASTQ/DIANA_0450_AH3JL3DSX3/Project_11969_E/Sample_Third-Transcriptome_IGO_11969_E_3', '/igo/staging/FASTQ/DIANA_0454_BH555MDMXY/Project_11969_E/Sample_Third-Transcriptome_IGO_11969_E_3'], 'Second_IGO_11969_E_2': ['/igo/staging/FASTQ/DIANA_0453_AHFKJ5DRXY/Project_11969_E/Sample_Second_IGO_11969_E_2', '/igo/staging/FASTQ/DIANA_0450_AH3JL3DSX3/Project_11969_E/Sample_Second_IGO_11969_E_2']}
# genome_dict = {"06265_8869_1_IGO_06265_AG_3":"Human_GeneticallyModified","Third-Transcriptome_IGO_11969_E_3":"Human", "Second_IGO_11969_E_2":"Mouse"}
# cmd = []
# for sample in sample_ID_list:
#     if genome_dict[sample] != "Human" and genome_dict[sample] != "Mouse":
#         genome_dict[sample] = "Mouse"
#     cmd.append(generate_cellranger_cmd(sample, "count", genome_dict[sample], fastq_file_list_dict[sample], "DIANA_0453_AHFKJ5DRXY"))
# print(cmd)
# test_result = ["bsub -J DIANA_0453_AHFKJ5DRXY_06265_8869_1_IGO_06265_AG_3_count_cellranger -o DIANA_0453_AHFKJ5DRXY_06265_8869_1_IGO_06265_AG_3_count_cellranger.out /igo/work/nabors/tools/cellranger-6.1.2/cellranger count --id=Sample_06265_8869_1_IGO_06265_AG_3__count --transcriptome=/igo/work/nabors/genomes/10X_Genomics/GEX/refdata-gex-GRCh38-2020-A --fastqs=/igo/staging/FASTQ/DIANA_0453_AHFKJ5DRXY/Project_06265_AG/Sample_06265_8869_1_IGO_06265_AG_3 --nopreflight --jobmode=lsf --mempercore=64 --disable-ui --maxjobs=200",
# "bsub -J DIANA_0453_AHFKJ5DRXY_Third-Transcriptome_IGO_11969_E_3_count_cellranger -o DIANA_0453_AHFKJ5DRXY_Third-Transcriptome_IGO_11969_E_3_count_cellranger.out /igo/work/nabors/tools/cellranger-6.1.2/cellranger count --id=Sample_Third-Transcriptome_IGO_11969_E_3__count --transcriptome=/igo/work/nabors/genomes/10X_Genomics/GEX/refdata-gex-mm10-2020-A --fastqs=/igo/staging/FASTQ/DIANA_0450_AH3JL3DSX3/Project_11969_E/Sample_Third-Transcriptome_IGO_11969_E_3,/igo/staging/FASTQ/DIANA_0454_BH555MDMXY/Project_11969_E/Sample_Third-Transcriptome_IGO_11969_E_3 --nopreflight --jobmode=lsf --mempercore=64 --disable-ui --maxjobs=200",
# "bsub -J DIANA_0453_AHFKJ5DRXY_Second_IGO_11969_E_2_count_cellranger -o DIANA_0453_AHFKJ5DRXY_Second_IGO_11969_E_2_count_cellranger.out /igo/work/nabors/tools/cellranger-6.1.2/cellranger count --id=Sample_Second_IGO_11969_E_2__count --transcriptome=/igo/work/nabors/genomes/10X_Genomics/GEX/refdata-gex-mm10-2020-A --fastqs=/igo/staging/FASTQ/DIANA_0453_AHFKJ5DRXY/Project_11969_E/Sample_Second_IGO_11969_E_2,/igo/staging/FASTQ/DIANA_0450_AH3JL3DSX3/Project_11969_E/Sample_Second_IGO_11969_E_2 --nopreflight --jobmode=lsf --mempercore=64 --disable-ui --maxjobs=200"]

# for i in range (3): 
#     assert(cmd[i] == test_result[i])
