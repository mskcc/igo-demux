# launch cell ranger pipeline (GE, VDJ) for 10X samples
# launch pipeline by recipe TODO how to decide for SCRI
# put result in /igo/stats/CELLRANGER/<run_ID>

import re
import sys
import os
import glob
import argparse
import json

from os.path import join
from os.path import basename
from os.path import abspath
from os.path import isdir
from subprocess import call

"""
input: sample_sheet object(for sample list), sequencer_and_run(for stats folder and fastq file location)
output: running cmd for cellranger by sample
"""
# config info for cellranger (reference, tools location)

# cellranger applications
CELLRANGER_COUNT = ' /igo/work/nabors/tools/cellranger-6.1.2/cellranger count '
CELLRANGER_VDJ = ' /igo/work/nabors/tools/cellranger-4.0.0/cellranger vdj '
CELLRANGER_ATAC_COUNT = ' /igo/work/nabors/tools/cellranger-atac-1.2.0/cellranger-atac count '
CELLRANGER_CNV = ' /igo/work/nabors/tools/cellranger-dna-1.1.0/cellranger-dna cnv '
ACCESS = 0o775
CELLRANGER_PROJECT_DATA = 'CELLRANGER_JSON.json'

# work folder
STATS_AREA = '/igo/stats/CELLRANGER/'

# cellranger command line options
OPTIONS = ' --nopreflight --jobmode=lsf --mempercore=64 --disable-ui --maxjobs=200'

# References and Transcriptomes
MOUSE_TRANSCRIPTOME = ' --transcriptome=/igo/work/nabors/genomes/10X_Genomics/GEX/refdata-gex-mm10-2020-A '
HUMAN_TRANSCRIPTOME = ' --transcriptome=/igo/work/nabors/genomes/10X_Genomics/GEX/refdata-gex-GRCh38-2020-A '
MOUSE_REFERENCE = ' --reference=/igo/work/nabors/genomes/10X_Genomics/VDJ/refdata-cellranger-vdj-GRCm38-alts-ensembl-2.2.0 '
HUMAN_REFERENCE = ' --reference=/igo/work/nabors/genomes/10X_Genomics/VDJ/refdata-cellranger-vdj-GRCh38-alts-ensembl-2.0.0 '
HUMAN_REFERENCE_ATAC = ' --reference=/igo/work/nabors/genomes/10X_Genomics/ATAC/refdata-cellranger-atac-GRCh38-1.0.1 '
MOUSE_REFERENCE_ATAC = ' --reference=/igo/work/nabors/genomes/10X_Genomics/ATAC/refdata-cellranger-atac-mm10-1.1.0 '
HUMAN_REFERENCE_CNV = ' --reference=/igo/work/nabors/10X_Genomics_references/CNV/refdata-GRCh38-1.0.0 '
MOUSE_REFERENCE_CNV = ' --reference=/igo/work/nabors/10X_Genomics_references/CNV/refdata-GRCm38-1.0.0 '

# 10X recipe list for different pipelines TODO 10X_Genomics_Visium, 10X_Genomics_Multiome
COUNT_FLAVORS = ['10X_Genomics_RNA', '10X_Genomics_GeneExpression', '10X_Genomics_GeneExpression-3', '10X_Genomics_GeneExpression-5']
VDJ_FLAVORS = ['10X_Genomics_VDJ']
ATAC_FLAVORS = ['10X_Genomics_ATAC']
RNA_PLUS_VDJ_FLAVORS = ['10X_Genomics_VDJ_GeneExpression', '10X_Genomics_Expression_VDJ', '10X_Genomics-Expression+VDJ']
CNV_FLAVORS = ['10X_Genomics_CNV']

"""
steps:
1. check which pipeline to run by recipe, if recipe not in the list above, skip for now
2. check whether there is previous fastq existing under /igo/staging/FASTQ (find_fastq_file)
3. generate corresponding commands based on genome
4. run stats
5. create josn file for the run and call endpoint to push to qc website
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

def gene_expression(project, run, samples, genome):
        tag = 'count'
        # GET TRANSCRIPTIION FILE
        if genome == 'human':
                transcriptome = HUMAN_TRANSCRIPTOME
        else:
                transcriptome = MOUSE_TRANSCRIPTOME
        # GO TO WORK AREA
        # work_area = '/igo/work/GCL/hiseq/Stats/CELLRANGER/MICHELLE_0243_cellranger/combined/' + project + '/'
        work_area = STATS_AREA + run + '/' + project + '/' 
        job_name = 'cellranger_' + tag + GAP + run + GAP + project + GAP
        os.chdir(work_area)
        print('<<<<<<<<<<<<<<<<<<<<      ' + work_area + '        >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>')
        for s in samples:
                # make FASTQ location variable
                fastqs = project_area + s
                print(fastqs)
                id_area = run + GAP + project  + GAP + s
                bsub = 'bsub -J ' + job_name + s + ' -o ' + id_area + GAP + tag + '.log'
                # sub_cellranger_count = bsub + CELLRANGER_COUNT + '--id=' + s + GAP + tag + transcriptome + "--fastqs=/igo/staging/FASTQ/DIANA_0429_AHFG3WDRXY/" + project + "
#/" + s + ',' + "/igo/staging/FASTQ/RUTH_0053_AH5NK7DSX3/" + project + '/' + s + ',' + fastqs + OPTIONS
                sub_cellranger_count = bsub + CELLRANGER_COUNT + '--id=' + s + GAP + tag + transcriptome + '--fastqs=' + fastqs + OPTIONS
                print(sub_cellranger_count)
                # execute cellranger count command
                # call(sub_cellranger_count, shell = True)
        create_json(job_name, run, samples, project, tag, work_area)
def cnv(project, run, samples, genome):
        tag = 'cnv'
        # GET TRANSCRIPTIION FILE
        if genome == 'human':
                reference = HUMAN_REFERENCE_CNV
        else:
                reference = MOUSE_REFERENCE_CNV
        # GO TO WORK AREA
        work_area = STATS_AREA + run + '/' + project + '/' 
        job_name = 'cellranger_' + tag + GAP + run + GAP + project + GAP
        os.chdir(work_area)
        print('<<<<<<<<<<<<<<<<<<<<      ' + work_area + '        >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>')
        for s in samples:
                # make FASTQ location variable
                fastqs = project_area + s
                print(fastqs)
                id_area = run + GAP + project  + GAP + s
                bsub = 'bsub -J ' + job_name + s + ' -o ' + id_area + GAP + tag + '.log'
                sub_cellranger_cnv = bsub + CELLRANGER_CNV + '--id=' + s + GAP + tag + reference + '--fastqs=' + fastqs + OPTIONS
                print(sub_cellranger_cnv)
                # execute cellranger cnv command
                call(sub_cellranger_cnv, shell = True)
        create_json(job_name, run, samples, project, tag, work_area)
        
        
def vdj(project, run, samples, genome):
        tag = 'vdj'
        # GET TRANSCRIPTIION FILE
        if genome == 'human':
                reference = HUMAN_REFERENCE
        else:
                reference = MOUSE_REFERENCE
                # GO TO WORK AREA
        work_area = STATS_AREA + run + '/' + project + '/' 
        job_name = 'cellranger_' + tag + GAP + run + GAP + project + GAP
        os.chdir(work_area)
        print('<<<<<<<<<<<<<<<<<<<<      ' + work_area + '        >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>')
        for s in samples:
                # make FASTQ location variable
                fastqs = project_area + s
                print(fastqs)
                id_area = run + GAP + project  + GAP + s
                bsub = 'bsub -J ' + job_name + s + ' -o ' + id_area + GAP + tag + '.log'
                sub_cellranger_vdj = bsub + CELLRANGER_VDJ + '--id=' + s + GAP + tag + reference + '--fastqs=' + fastqs + OPTIONS
                print(sub_cellranger_vdj)
                # execute cellranger vdj command
                call(sub_cellranger_vdj, shell = True)
        create_json(job_name, run, samples, project, tag, work_area)

def atac(project, run, samples, genome):
        tag = 'atac_count'
        # GET TRANSCRIPTIION FILE
        if genome == 'human':
                reference = HUMAN_REFERENCE_ATAC
        else:
                reference = MOUSE_REFERENCE_ATAC
        # GO TO WORK AREA
        work_area = STATS_AREA + run + '/' + project + '/' 
        job_name = 'cellranger_' + tag + GAP + run + GAP + project + GAP
        os.chdir(work_area)
        print('<<<<<<<<<<<<<<<<<<<<      ' + work_area + '        >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>')
        for s in samples:
                # make FASTQ location variable
                fastqs = project_area + s
                print(fastqs)
                id_area = run + GAP + project  + GAP + s
                bsub = 'bsub -J ' + job_name + s + ' -o ' + id_area + GAP + tag + '.log'
                sub_cellranger_atac = bsub + CELLRANGER_ATAC_COUNT + '--id=' + s + GAP + tag + reference + '--fastqs=' + fastqs + OPTIONS
                print(sub_cellranger_atac)
                # execute cellranger atac command
                call(sub_cellranger_atac, shell = True)
        create_json(job_name, run, samples, project, tag, work_area)
        

def create_json(job_name, run, samples, project, tag, work_area):
        send_json = {}
        send_json['samples'] = []
        print(type(samples))
        print(samples)
        for x in samples:
                send_json['samples'].append({'sample':x, 'type':tag, 'project':project, 'run':run})
                        
        json_data_file = 'cellranger_json___' + run + GAP + project + '.json'
        with open(json_data_file, 'w') as jfile:
                json.dump(send_json, jfile)
        
        bsub_json = 'bsub -J create_json___' + run + GAP + project + ' -o ' + 'create_json___' + run + GAP + project + '.log -w \"done(' + job_name + '*)\" sh /home/igo/Scripts/PicardScripts/send_json_data.sh ' + work_area + ' ' + json_data_file
        print(bsub_json)
        call(bsub_json, shell = True)
'''        
if __name__ == '__main__':
        parser = argparse.ArgumentParser(description='Given Cellranger FASTQs, run Cellranger Pipeline, and generate stats for IGO-QC')
        parser.add_argument('--project-dir', type = str, required = True, help = 'parent directory with subdirs per sample containing FASTQs')
        parser.add_argument('--recipe', type = str, required = False, default = '10X_Genomics_GeneExpression', help = 'recipe for the 10X Genomics Sequenced FASTQs')
        parser.add_argument('--genome', type = str, required = False, default = 'human', help = 'genome of the FASTQ data')
        args = parser.parse_args()

        # Get arguments
        project_area = args.project_dir
        recipe = args.recipe
        genome = args.genome
        

        # grab the IGO project ID and run ID from the project-dir argument
        run = project_area.split('/')[4]
        project = project_area.split('/')[5]
        
        # CREATE RUN FOLDER AND PROJECT FOLDER IF NOT ALREADY THERE
        os.chdir(STATS_AREA)
        runs = next(os.walk('.'))[1]
        if run not in runs:
                os.mkdir(run, ACCESS)
                
        stats_and_run = STATS_AREA + run
        os.chdir(stats_and_run)
        projects = next(os.walk('.'))[1]
        if project not in projects:
                os.mkdir(project, ACCESS)

                
        # GO TO project ID LOCATION
        os.chdir(project_area)
        
        # GRAB SAMPLES
        samples = next(os.walk('.'))[1]
        print(type(samples))
        # print(samples)
        print(recipe)
        # print(genome)
        # GO TO THE CORRECT FUNCTION
        if recipe in COUNT_FLAVORS:
                print('lets go ahead and do count...')
                gene_expression(project, run, samples, genome)
        if recipe in CNV_FLAVORS:
                print('lets do CNV this time...')
                cnv(project, run, samples, genome)      
        if recipe in VDJ_FLAVORS:
                print('we are about to do VDJ...')
                vdj(project, run, samples, genome)
        if recipe in ATAC_FLAVORS:
                print('Need to do ATAC...')
                atac(project, run, samples, genome)
        if recipe in RNA_PLUS_VDJ_FLAVORS:
                print('oh snap!  we need to run both COUNT and VDJ')
                gene_expression(project, run, samples, genome)
                vdj(project, run, samples, genome)
                
        
        print('Fin... 4 now')   
'''
# test = {}
# test1 = os.listdir("/Users/luc/Documents/")
# for i in test1:
#     curret = "/Users/luc/Documents/" + i + "/"
#     if os.path.isdir(curret):
#         test[i] = os.listdir(curret)
# with open("result.txt", 'w') as _file:
#             _file.write(str(test1))


# Sample_ID,Sample_Plate,Sample_Well,I7_Index_ID,index,index2,Sample_Project,Description
# AT8286N_IGO_04540_R_1,Human,HumanWholeGenome,UDI0079,ACAGGCGC,AGGCAGAG,Project_04540_R,weigeltb@mskcc.org



sample_ID_list = ["06265_8869_1_IGO_06265_AG_3","Third-Transcriptome_IGO_11969_E_3", "Second_IGO_11969_E_2"])
fastq_file_list_dict = {'06265_8869_1_IGO_06265_AG_3': ['/igo/staging/FASTQ/DIANA_0453_AHFKJ5DRXY/Project_06265_AG/Sample_06265_8869_1_IGO_06265_AG_3'], 'Third-Transcriptome_IGO_11969_E_3': ['/igo/staging/FASTQ/DIANA_0450_AH3JL3DSX3/Project_11969_E/Sample_Third-Transcriptome_IGO_11969_E_3', '/igo/staging/FASTQ/DIANA_0454_BH555MDMXY/Project_11969_E/Sample_Third-Transcriptome_IGO_11969_E_3'], 'Second_IGO_11969_E_2': ['/igo/staging/FASTQ/DIANA_0453_AHFKJ5DRXY/Project_11969_E/Sample_Second_IGO_11969_E_2', '/igo/staging/FASTQ/DIANA_0450_AH3JL3DSX3/Project_11969_E/Sample_Second_IGO_11969_E_2']}


# /igo/work/bin/cellranger-6.0.1/cellranger count --id=06265_8869_1_IGO_06265_AG_3 --transcriptome=/igo/work/nabors/genomes/10X_Genomics/GEX/refdata-gex-GRCh38-2020-A --fastqs=/igo/staging/FASTQ/DIANA_0453_AHFKJ5DRXY/Project_06265_AG/Sample_06265_8869_1_IGO_06265_AG_3 --nopreflight --jobmode=lsf --mempercore=64 --disable-ui --maxjobs=200
# /igo/work/bin/cellranger-6.0.1/cellranger count --id=Second_IGO_11969_E_2 --transcriptome=/igo/work/nabors/genomes/10X_Genomics/GEX/refdata-gex-mm10-2020-A --fastqs=/igo/staging/FASTQ/DIANA_0453_AHFKJ5DRXY/Project_11969_E/Sample_Second_IGO_11969_E_2,/igo/staging/FASTQ/DIANA_0450_AH3JL3DSX3/Project_11969_E/Sample_Second_IGO_11969_E_2 --nopreflight --jobmode=lsf --mempercore=64 --disable-ui --maxjobs=200
# /igo/work/bin/cellranger-6.0.1/cellranger count --id=Third-Transcriptome_IGO_11969_E_3 --transcriptome=/igo/work/nabors/genomes/10X_Genomics/GEX/refdata-gex-mm10-2020-A --fastqs=/igo/staging/FASTQ/DIANA_0450_AH3JL3DSX3/Project_11969_E/Sample_Third-Transcriptome_IGO_11969_E_3,/igo/staging/FASTQ/DIANA_0454_BH555MDMXY/Project_11969_E/Sample_Third-Transcriptome_IGO_11969_E_3 --nopreflight --jobmode=lsf --mempercore=64 --disable-ui --maxjobs=200