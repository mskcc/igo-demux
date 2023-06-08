import pandas as pd
import re
import sys
import os
import glob
import json
import subprocess
from os.path import join, basename, abspath, isdir
from subprocess import call
from cellranger import config_dict, OPTIONS

CONFIG_AREA = "/igo/stats/Multi_config/"
DRIVE_LOCATION = "/skimcs/mohibullahlab/LIMS/LIMS_cellranger_multi/"

# config file class. It contains all the information needed for config file and can create csv file base on those info
class Multi_Config:
    def __init__(self):
        self.gene_expression = {}
        self.lirbaries = {}
        self.vdj = "EMPTY"
        self.samples = "EMPTY"
        self.features = "EMPTY"
    
    def write_to_csv(self, name_of_file):
        with open(name_of_file,'w') as file:
            file.write("[gene-expression]\n")
            for key, value in self.gene_expression.items():
                file.write("{},{}\n".format(key, value))
            
            file.write("\n[libraries]\n")
            file.write("fastq_id,fastqs,feature_types\n")
            
            for key, value in self.lirbaries.items():
                for i in value[0]:
                    file.write("{},{},{}\n".format(key, i, value[1]))
            if self.vdj != "EMPTY":
                file.write("\n[vdj]\nreference,{}\n".format(self.vdj))

            if self.features != "EMPTY":
                file.write("\n[feature]\nreference,{}\n".format(self.features))

            if self.samples != "EMPTY":
                file.write("\n[samples]\nsample_id,cmo_ids\n")
                for item in self.samples:
                    file.write("{},{}\n".format(item[0], item[1]))

# read ch file from shared drive and generate config file per sample and return sample to sample info also
# default all hash tag are totalseq B from biolegend
def ch_file_generation(project_id):
    in_file_location = DRIVE_LOCATION + project_id + "/" + os.listdir("DRIVE_LOCATION + project_id")[0]
    df = pd.read_excel(in_file_location, engine="openpyxl")
    line_number = df[df[df.columns[0]] == "Your Submission:"].index.values
    df = pd.read_excel(in_file_location, engine="openpyxl", skiprows=line_number + 1, header=line_number + 1)
    sample_tag_dict = pd.Series(df['Hashtag Name'].values,index=df['Sample Name']).to_dict()
    tag_seq_dict = pd.Series(df['Hashtag sequence'].values,index=df['Hashtag Name']).to_dict()
    sample_lst = set(df["Sample Name in IGO"].tolist())
    sample_dict = {}
    for sample in sample_lst:
        sample_dict[sample] = df[df["Sample Name in IGO"] == sample]["Sample Name"].tolist()
        for i in range(len(sample_dict[sample])):
            sample_dict[sample][i] = [sample_dict[sample][i], sample_tag_dict[sample_dict[sample][i]]]

    # write ch config file for this project
    file_name = "{}Project_{}/Project_{}_ch.csv".format(CONFIG_AREA, project_id, project_id)
    if not os.path.exists(os.path.dirname(file_name)):
        os.makedirs(os.path.dirname(file_name))

    with open(file_name,'w') as file:
        file.write("id,name,read,pattern,sequence,feature_type\n")
        for key, value in tag_seq_dict.items():
            file.write("{},{},R2,5PNNNNNNNNNN(BC),{},Multiplexing Capture\n".format(key, key, value))

    return(sample_dict)

def gather_config_info(sample_dict, genome, IGO_ID):
    """
    sample_dict contains all the information about samples, sample name, project_ID and recipe name
    example: {"ge":"LJ01_IGO_14396_1", "vdj": "LJ01_VDJ_IGO_14396_C_1", "fb":"", "ch":""}
    one file per sample named with IGO ID (GEX sample ID)
    """
    # library use fastq id as key, then followed by a list which first item is list of fastq path in case top up and second is feature type
    # fb reference file should have format as following: /igo/stats/Multi_config/Project_12345/Project_12345_fb.csv
    # ch reference file should have format as following: /igo/stats/Multi_config/Project_12345/Project_12345_ch.csv
    # how to record vdj-t and vdj-b?
    project_ID = "_".join(IGO_ID.split("IGO_")[1].split("_")[:-1])
    sample_name = IGO_ID.split("IGO_")[0]
    config = Multi_Config()

    config.gene_expression["reference"] = cellranger.config_dict["count"]["genome"][genome]
    if "vdj" in sample_dict.keys():
        config.vdj = cellranger.config_dict["vdj"]["genome"][genome]
    
    # if feature barcoding invovled, add feature list file path
    if "fb" in sample_dict.keys():
        config.features = CONFIG_AREA + "Project_{}/Project_{}_fb.csv".format(project_ID, project_ID)
    # if cell hashing invovled, add cmo-set file path and get sample info from file, id as sample name and name as hashtag name
    if "ch" in sample_dict.keys():
        config.gene_expression["cmo-set"] = CONFIG_AREA + "Project_{}/Project_{}_ch.csv".format(project_ID, project_ID)
        config.samples = ch_file_generation(project_ID)[sample_name]

    # find fastq files for each sample and append information into config["libraries"]
    sample_list = []
    for i in sample_dict.values():
        sample_list.append(i)
    fastq_list = find_fastq_file(sample_list)
    for key, value in sample_dict.items():
        if key == "ge":
            config.lirbaries[value] = [fastq_list[value], "Gene Expression"]
        elif key == "vdj":
            config.lirbaries[value] = [fastq_list[value], "VDJ"]
        elif key == "fb":
            config.lirbaries[value] = [fastq_list[value], "Antibody Capture"]
        elif key == "ch":
            config.lirbaries[value] = [fastq_list[value], "Multiplexing Capture"]
       
    return config

# example config Class
# test = Multi_Config()
# test.gene_expression["reference"] = "/igo/work/nabors/genomes/10X_Genomics/GEX/refdata-gex-GRCh38-2020-A"
# test.lirbaries = {"WO9112_NoPeptide_IGO_14514_1":[["/igo/staging/FASTQ/MICHELLE_0643_BHKL3GDMXY/Project_14514/Sample_WO9112_NoPeptide_IGO_14514_1"], "Gene Expression"], "WO9112_NoPeptide_VDJ_IGO_14514_C_1":[["/igo/staging/FASTQ/MICHELLE_0643_BHKL3GDMXY/Project_14514_C/Sample_WO9112_NoPeptide_VDJ_IGO_14514_C_1", "fake test"],"VDJ"]}
# test.vdj = "/igo/work/genomes/10X_Genomics/VDJ/refdata-cellranger-vdj-GRCh38-alts-ensembl-7.0.0"
# test.features = "/igo/stats/Multi_config/Project_14514_fb.csv"
# test.samples = {"test1":"B301", "test2":"B302"}

# test.write_to_csv("test.csv")


if __name__ == '__main__':
    # input as sample dict, genome and IGO ID
    # Usage: python cellranger_multi.py [sample_dict] [Human/Mouse] [IGO_ID]
    sample_dict = sys.argv[1]
    genome = sys.argv[2]
    IGO_ID = sys.argv[3]
    config = gather_config_info(sample_dict, genome, IGO_ID)
    project_ID = "_".join(IGO_ID.split("IGO_")[1].split("_")[:-1])
    file_name = "{}/Project_{}/{}.csv".format(CONFIG_AREA, project_ID, IGO_ID)
    config.write_to_csv(file_name)

    cmd = "bsub -J {}_multi -o {}_multi.out {} --id={} --csv={} --nopreflight --jobmode=lsf --mempercore=64 --disable-ui --maxjobs=200".format(IGO_ID, IGO_ID, config_dict["multi"]["tool"], IGO_ID, file_name)
    print(cmd)
