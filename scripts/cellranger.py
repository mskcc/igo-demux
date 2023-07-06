# launch cell ranger pipeline (GE, VDJ, ATAC....) for 10X samples by recipe
# put result in /igo/stats/CELLRANGER/<run_ID>

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
import scripts.get_sequencing_read_data
import scripts.cellranger_spatial

"""
input: sample_sheet object(for sample list and essential info), sequencer_and_run(for stats folder and fastq file location)
output: running cmd for cellranger by sample
"""

# work folder
STATS_AREA = "/igo/stats/CELLRANGER/"

# config info 
ACCESS = 0o775
config_dict = {
    "count": {
        "tool": " /igo/work/nabors/tools/cellranger-7.0.0/cellranger count ",
        "genome": {
            "Human": " --transcriptome=/igo/work/nabors/genomes/10X_Genomics/GEX/refdata-gex-GRCh38-2020-A ",
            "Mouse": " --transcriptome=/igo/work/nabors/genomes/10X_Genomics/GEX/refdata-gex-mm10-2020-A "
        }
    },
    "vdj": {
        "tool": " /igo/work/nabors/tools/cellranger-7.0.0/cellranger vdj ",
        "genome": {
            "Human": " --reference=/igo/work/genomes/10X_Genomics/VDJ/refdata-cellranger-vdj-GRCh38-alts-ensembl-7.0.0 ",
            "Mouse": " --reference=/igo/work/genomes/10X_Genomics/VDJ/refdata-cellranger-vdj-GRCm38-alts-ensembl-7.0.0 "
        }
    },
    "atac_count": {
        "tool": " /igo/work/nabors/tools/cellranger-atac-2.1.0/cellranger-atac count ",
        "genome": {
            "Human": " --reference=/igo/work/nabors/genomes/10X_Genomics/ATAC/refdata-cellranger-atac-GRCh38-1.0.1 ",
            "Mouse": " --reference=/igo/work/nabors/genomes/10X_Genomics/ATAC/refdata-cellranger-atac-mm10-1.1.0 "
        }
    },
    "cnv": {
        "tool": " /igo/work/nabors/tools/cellranger-dna-1.1.0/cellranger-dna cnv ",
        "genome": {
            "Human": " --reference=/igo/work/nabors/10X_Genomics/CNV/refdata-GRCh38-1.0.0 ",
            "Mouse": " --reference=/igo/work/nabors/10X_Genomics/CNV/refdata-GRCm38-1.0.0 "
        }
    },
    "multi": {
        "tool": " /igo/work/nabors/tools/cellranger-7.0.0/cellranger multi "
    },
    "arc": {
        "tool": " /igo/work/bin/cellranger-arc-2.0.0/cellranger-arc count ",
        "genome": {
            "Human": " --reference=/igo/work/nabors/genomes/10X_Genomics/ARC/refdata-cellranger-arc-GRCh38-2020-A-2.0.0 ",
            "Mouse": " --reference=/igo/work/nabors/genomes/10X_Genomics/ARC/refdata-cellranger-arc-mm10-2020-A-2.0.0 "
        }
    },
    "spaceranger": {
        "tool": " /igo/work/nabors/tools/spaceranger-2.0.0/spaceranger count ",
        "genome": {
            "Human": " --reference=/igo/work/nabors/genomes/10X_Genomics/GEX/refdata-gex-GRCh38-2020-A ",
            "Mouse": " --reference=/igo/work/nabors/genomes/10X_Genomics/spatial_gex/refdata-gex-mm10-2020-A"
        },
        "probe": {
            "Human": "/igo/work/nabors/genomes/10X_Genomics/spatial_gex/Visium_Human_Transcriptome_Probe_Set_v1.0_GRCh38-2020-A.csv",
            "Mouse": "/igo/work/nabors/tools/spaceranger-2.0.0/external/tenx_feature_references/targeted_panels/Visium_Mouse_Transcriptome_Probe_Set_v1.0_mm10-2020-A.csv"
        }
    }
}

# cellranger command line options
OPTIONS = " --nopreflight --jobmode=lsf --mempercore=64 --disable-ui --maxjobs=200"

# 10X recipe list for different pipelines
COUNT_FLAVORS = ["10X_Genomics_GeneExpression-3", "10X_Genomics_GeneExpression-5"]
VDJ_FLAVORS = ["10X_Genomics_VDJ"]
ATAC_FLAVORS = ["10X_Genomics_ATAC"]
CNV_FLAVORS = ["10X_Genomics_CNV"]
ARC_FLAVORS = ["10X_Genomics_Multiome", "10X_Genomics_Multiome_ATAC", "10X_Genomics_Multiome_GeneExpression"]
SPATIAL_FLAVORS = ["10X_Genomics_Visium"]

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
    if recipe in ARC_FLAVORS:
        tag = "arc"
    if recipe in SPATIAL_FLAVORS:
        tag = "spatial"
    return tag

# return tag and genome according to sample_ID for SCRI samples, all SCRI samples are starting with Project_12437
# eg: SD-1680_Patient_D_nucseq_H_VDJ_IGO_12437_AN_5 will given tag as vdj, genome as Human
# eg: SDtest_IGO_12437_AN_4 will given tag as Skip, genome as na
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
    cellranger_cmd = "{}--id=Sample_{}__{}".format(tool, sample_ID, tag) + transcriptome + "--fastqs=" + ",".join(fastq_file_path) + OPTIONS
    job_name = "{}_{}_{}_{}_cellranger".format(sequencer_and_run, project_ID, sample_ID, tag)
    bsub_cmd = "bsub -J {} -o {}.out{}".format(job_name, job_name, cellranger_cmd) 
    return bsub_cmd
        
def create_json(send_json, sequencer_and_run, project, tag, work_area):  
    job_id = sequencer_and_run + "_" + project            
    json_data_file = "cellranger_json___" + sequencer_and_run + "__" + project + ".json"
    with open(json_data_file, "w") as jfile:
        json.dump(send_json, jfile)
        
    bsub_json = "bsub -J create_json___{} -o create_json___{}.log -w \"done({}*)\" sh /home/igo/Scripts/PicardScripts/send_json_data.sh {} {}".format(job_id, job_id, job_id, work_area, json_data_file)
    print(bsub_json)
    subprocess.run(bsub_json, shell = True)

def create_library_csv_file(ge_sample_path, atac_sample_path, sample_ID):
    """
    ge_sample_path and atac_sample_path will be list just in case if top up is performed
    """
    with open("Sample_{}.csv".format(sample_ID),"a") as file:
        file.write("fastqs,sample,library_type\n")
        for ge in ge_sample_path:
            file.write("{},{},Gene Expression\n".format(ge, sample_ID))
        for atac in atac_sample_path:
            file.write("{},{},Chromatin Accessibility\n".format(atac, sample_ID))

def get_sequencer_runID(fastq_path):
    runID = fastq_path.split("/")[4]
    sequencer = runID.split("_")[0].lower()
    return sequencer, runID

def multiome_valid(fastq_list):
    """
    check whether the list of fastq files contain both GE and ATAC
    return a list which first item will be YES or NO for validation result
    if yes, second item will be GE fastq list and third item will be ATAC fastq list
    """
    sequencer_prefix = "/igo/sequencers/"
    is_valid = "NO"
    ge_list = []
    atac_list = []
    if len(fastq_list) == 1:
        return [is_valid, ge_list, atac_list]
    else:
        for fastq in fastq_list:
            # find corresponding run folder and check if is atac run
            # append the fastq file to corresponding list
            sequencer, runID = get_sequencer_runID(fastq)
            # find sequenceing path
            if sequencer == "pepe" or sequencer == "amelie":
                sequencer_path_pre = sequencer_prefix + sequencer + "/output"
            else:
                sequencer_path_pre = sequencer_prefix + sequencer
            run_list = os.listdir(sequencer_path_pre)
            for run in run_list:
                if re.match(".*" + runID, run):
                    sequencer_path = sequencer_path_pre + "/" + run
                    atac = scripts.get_sequencing_read_data.get_sequencing_read_data(sequencer_path)[0]
                    if atac:
                        atac_list.append(fastq)
                    else:
                        ge_list.append(fastq)
                    break

    if len(ge_list) > 0 and len(atac_list) > 0:
        is_valid = "YES"
    
    return [is_valid, ge_list, atac_list]

# Main function: launch cellranger cmd by given samplesheet object and sequencer_and_run
def launch_cellranger(sample_sheet, sequencer_and_run):
    # get parameters from sample_sheet
    # dictionary of Sample_ID->Project
    sample_project_dict = pd.Series(sample_sheet.df_ss_data["Sample_Project"].values,index=sample_sheet.df_ss_data["Sample_ID"]).to_dict()
    # dictionary of project->sample_ID
    project_sample_dict = {}
    for sample_ID, project_ID in sample_project_dict.items():
        if project_ID in project_sample_dict.keys():
            project_sample_dict[project_ID].append(sample_ID)
        else:
            project_sample_dict[project_ID] = [sample_ID]
    # dictionary of sample_ID->recipe
    sample_recipe_dict = pd.Series(sample_sheet.df_ss_data["Sample_Well"].values,index=sample_sheet.df_ss_data["Sample_ID"]).to_dict()
    # dictionary of sample_ID->genome
    sample_genome_dict = pd.Series(sample_sheet.df_ss_data["Sample_Plate"].values,index=sample_sheet.df_ss_data["Sample_ID"]).to_dict()
    # dictionary of sample_ID->fastq_list
    sample_ID_list = list(sample_project_dict.keys())
    sample_fastqfile_dict = find_fastq_file(sample_ID_list)

    for project in project_sample_dict.keys():
        send_json = {}
        send_json["samples"] = []
        # CREATE RUN FOLDER AND PROJECT FOLDER IF NOT ALREADY THERE
        os.chdir(STATS_AREA)
        runs = next(os.walk("."))[1]
        if sequencer_and_run not in runs:
            os.mkdir(sequencer_and_run, ACCESS)
                    
        stats_and_run = STATS_AREA + sequencer_and_run
        os.chdir(stats_and_run)
        projects = next(os.walk("."))[1]
        if project not in projects:
            os.mkdir(project, ACCESS)
        work_area = stats_and_run + "/" + project + "/" 
        # GO TO project ID LOCATION to start cellranger command
        os.chdir(work_area)

        # SCRI samples don't need to be pushed onto qc website
        if "Project_12437" not in project:
            sample_list = project_sample_dict[project]
            # call cellranger for each sample and append info to json dict
            for sample in sample_list:
                if sample_genome_dict[sample] != "Human" and sample_genome_dict[sample] != "Mouse":
                    sample_genome_dict[sample] = "Mouse"
                tag = get_tag(sample_recipe_dict[sample])
                # if recipe within the tool being set up, lanuch cellranger
                if tag == "arc":
                    validation = multiome_valid(sample_fastqfile_dict[sample])
                    if validation[0] == "YES":
                        create_library_csv_file(validation[1], validation[2], sample)
                        tool = config_dict[tag]["tool"]
                        transcriptome = config_dict[tag]["genome"][sample_genome_dict[sample]]
                        cmd = "{}--id=Sample_{}{}".format(tool, sample, transcriptome) + "--libraries={}Sample_{}.csv".format(work_area, sample) + OPTIONS
                        bsub_cmd = "bsub -J {}_{}_{}_ARC -o {}_ARC.out{}".format(sequencer_and_run, project, sample, sample, cmd)
                        print(bsub_cmd)
                        subprocess.run(bsub_cmd, shell=True)
                    else:
                        print("Multiome sample set not complete yet")
                elif tag == "spatial":
                    scripts.cellranger_spatial.copy_tiff(project)
                    sample_info = scripts.cellranger_spatial.Spatial_sample(sample)
                    tool = config_dict[tag]["tool"]
                    transcriptome = config_dict[tag]["genome"][sample_genome_dict[sample]]
                    cmd = "{}--id=Sample_{}{}".format(tool, sample, transcriptome) + "--fastqs=" + ",".join(sample_fastqfile_dict[sample]) + " --image={} --slide={} --area={}".format(sample_info.tiff_image, sample_info.chip_id, sample_info.chip_position)
                    if sample_info.preservation == "FFPE":
                        probe = config_dict[tag]["probe"][sample_genome_dict[sample]]
                        cmd = cmd + "--probe-set={}".format(probe)
                    bsub_cmd = "bsub -J {}_{}_{}_SPATIAL -o {}_SPATIAL.out{}{}".format(sequencer_and_run, project, sample, sample, cmd, OPTIONS)
                    print(bsub_cmd)
                    # subprocess.run(bsub_cmd, shell=True)
                
                elif tag != "Skip":
                    cmd = generate_cellranger_cmd(sample, tag, sample_genome_dict[sample], sample_fastqfile_dict[sample], sequencer_and_run)
                    print(cmd)
                    subprocess.run(cmd, shell=True)
                    send_json["samples"].append({"sample":"Sample_" + sample, "type":tag, "project":project, "run":sequencer_and_run})
            if send_json["samples"]:
                create_json(send_json, sequencer_and_run, project, tag, work_area)
        else:
            sample_list = project_sample_dict[project]
            # call cellranger for each sample
            for sample in sample_list:
                tag, genome = get_SCRI_tag(sample)
                # if recipe within the tool being set up, lanuch cellranger
                if tag != "Skip" and genome != "na":
                    cmd = generate_cellranger_cmd(sample, tag, genome, sample_fastqfile_dict[sample], sequencer_and_run)
                    cmd = cmd + " --include-introns=true"  # SCRI samples always have include-introns true
                    print(cmd)
                    subprocess.run(cmd, shell=True)

# lanuch cellranger by given project_directory eg: /igo/staging/FASTQ/RUTH_0141_AH27NGDSX5/Project_13586_B
def lanuch_by_project(project_directory, recipe, species):
    # get sample_ID list
    sample_list_ori = os.listdir(project_directory)
    sample_list = []
    for sample in sample_list_ori:
        # remove Sample_ prefix
        sample_list.append(sample[7:])
    # get project and run info from project_directory
    project = project_directory.split("/")[5]
    sequencer_and_run = project_directory.split("/")[4]
    sample_fastqfile_dict = find_fastq_file(sample_list)
    tag = get_tag(recipe)
    send_json = {}
    send_json["samples"] = []
    # CREATE RUN FOLDER AND PROJECT FOLDER IF NOT ALREADY THERE
    os.chdir(STATS_AREA)
    runs = next(os.walk("."))[1]
    if sequencer_and_run not in runs:
        os.mkdir(sequencer_and_run, ACCESS)
                    
    stats_and_run = STATS_AREA + sequencer_and_run
    os.chdir(stats_and_run)
    projects = next(os.walk("."))[1]
    if project not in projects:
        os.mkdir(project, ACCESS)
    work_area = stats_and_run + "/" + project + "/" 
    # GO TO project ID LOCATION to start cellranger command
    os.chdir(work_area)

    # call cellranger for each sample and append info to json dict
    for sample in sample_list:
        # if recipe within the tool being set up, lanuch cellranger
        if tag == "arc":
            validation = multiome_valid(sample_fastqfile_dict[sample])
            if validation[0] == "YES":
                create_library_csv_file(validation[1], validation[2], sample)
                tool = config_dict[tag]["tool"]
                transcriptome = config_dict[tag]["genome"][species]
                cmd = "{}--id=Sample_{}{}".format(tool, sample, transcriptome) + "--libraries={}/Sample_{}.csv".format(work_area, sample) + OPTIONS
                bsub_cmd = "bsub -J {}_{}_{}_ARC -o {}_ARC.out{}".format(sequencer_and_run, project, sample, sample, cmd)
                print(bsub_cmd)
                subprocess.run(cmd, shell=True)
            else:
                print("Multiome sample not finished yet")
                print(validation)
        elif tag != "Skip":
            cmd = generate_cellranger_cmd(sample, tag, species, sample_fastqfile_dict[sample], sequencer_and_run)
            print(cmd)
            subprocess.run(cmd, shell=True)
            send_json["samples"].append({"sample":"Sample_" + sample, "type":tag, "project":project, "run":sequencer_and_run})
    if send_json["samples"]:
        create_json(send_json, sequencer_and_run, project, tag, work_area)


if __name__ == '__main__':
    # launch cellranger commands by project
    # Usage: python cellranger.py [project_directory] [recipe] [species]
    # example: python cellranger.py /igo/staging/FASTQ/RUTH_0141_AH27NGDSX5/Project_13586_B 10X_Genomics_GeneExpression-3 Human
    project_directory = sys.argv[1]
    recipe = sys.argv[2]
    species = sys.argv[3]
    lanuch_by_project(project_directory, recipe, species)
