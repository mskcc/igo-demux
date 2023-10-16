import pandas as pd
import os
import subprocess
import glob
from subprocess import call
import argparse
from collections import OrderedDict
import requests
import json

# copy those config info from cellranger.py so this code can run independently.
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
    "multi": {
        "tool": " /igo/work/nabors/tools/cellranger-7.0.0/cellranger multi "
    }
}

# cellranger command line options
OPTIONS = " --nopreflight --jobmode=lsf --mempercore=64 --disable-ui --maxjobs=200"
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

CONFIG_AREA = "/igo/stats/Multi_config/"
DRIVE_LOCATION = "/igo/work/igo/Cellranger_Multi_Config/"
ORIGIN_DRIVE_LOCATION = "/rtssdc/mohibullahlab/LIMS/LIMS_cellranger_multi/"
BAMTOFASTQ = "/igo/work/nabors/tools/cellranger-7.0.0/lib/bin/bamtofastq"
STATS_AREA = "/igo/stats/PIPELINE/"
# endpoint for cellranger multi
ENDPOINT= "https://igolims.mskcc.org:8443/LimsRest/getTenxSampleInfo?requestId="

# config file class. It contains all the information needed for config file and can create csv file base on those info
class Multi_Config:
    def __init__(self):
        self.name = "EMPTY"  # gene expression sample name
        self.gene_expression = OrderedDict()   # genome info for the ge sample
        self.lirbaries = OrderedDict()   # location info for fastq files
        self.vdj = "EMPTY"    # genome info for the vdj sample
        self.samples = "EMPTY"  # sub_samples for cell hashing
        self.features = "EMPTY"   # fb info for fb sample, provided by user
        self.sub_sample_info = {} # cell number and fastq path as list for sub samples if include cell hashing
        self.ge_reads_number = 0  # gene expression fastq reads number if needed
    
    # write all info from the class to csv file
    def write_to_csv(self, name_of_file):
        with open(name_of_file,'w') as file:
            file.write("[gene-expression]\n")
            for key, value in self.gene_expression.items():
                file.write("{},{}\n".format(key, value))
            
            file.write("\n[libraries]\nfastq_id,fastqs,feature_types\n")
            
            for key, value in self.lirbaries.items():
                for i in value[0]:
                    file.write("{},{},{}\n".format(key, i, value[1]))
            if self.vdj != "EMPTY":
                file.write("\n[vdj]\nreference,{}\n".format(self.vdj))

            if self.features != "EMPTY":
                file.write("\n[feature]\nreference,{}\n".format(self.features))

            if self.samples != "EMPTY":
                file.write("\n[samples]\nsample_id,cmo_ids\n")
                for key, value in self.samples.items():
                    file.write("{},{}\n".format(key, value))

    # write to csv for ch and ge only, prepare for bam2fastq
    def write_ch_ge_only_to_csv(self, name_of_file):
        with open(name_of_file,'w') as file:
            file.write("[gene-expression]\n")
            for key, value in self.gene_expression.items():
                file.write("{},{}\n".format(key, value))
            
            file.write("\n[libraries]\nfastq_id,fastqs,feature_types\n")
            
            for key, value in self.lirbaries.items():
                if value[1] == "Gene Expression" or value[1] == "Multiplexing Capture":
                    for i in value[0]:
                        file.write("{},{},{}\n".format(key, i, value[1]))
                
            file.write("\n[samples]\nsample_id,cmo_ids\n")
            for key, value in self.samples.items():
                file.write("{},{}\n".format(key, value))

    # write to csv without ch info after bam2fastq for each sub sample
    def new_config_and_generate_cmd(self):
        for key, value in self.sub_sample_info.items():
            name_of_file = "{}Project_{}/{}_{}.csv".format(CONFIG_AREA, "_".join(self.name.split("IGO_")[1].split("_")[:-1]), self.name, key)
            with open(name_of_file,'w') as file:
                file.write("[gene-expression]\n")
                file.write("{},{}\n".format("reference", self.gene_expression["reference"]))
                file.write("force-cells, {}\ncheck-library-compatibility,FALSE\n".format(value[0]))
                
                file.write("\n[libraries]\nfastq_id,fastqs,feature_types\n")
                
                for key1, value1 in self.lirbaries.items():
                    if value1[1] == "Gene Expression":
                        file.write("bamtofastq,{},{}\n".format(value[1] , value1[1]))
                    elif value1[1] != "Multiplexing Capture":
                        for i in value1[0]:
                            file.write("{},{},{}\n".format(key1, i, value1[1]))
                if self.vdj != "EMPTY":
                    file.write("\n[vdj]\nreference,{}\n".format(self.vdj))

                if self.features != "EMPTY":
                    file.write("\n[feature]\nreference,{}\n".format(self.features))
            
            # generate cmd for final cellranger
            # wait bam2fastq finish before excute bsub command
            # cmd = "bsub -J {}_{}_multi -o {}_{}_multi.out -w \"done(*_bamtofastq)\"{}--id={}_{} --csv={}{}".format(self.name, key, self.name, key, config_dict["multi"]["tool"], self.name, key, name_of_file, OPTIONS)
            # no need to wait bam2fastq finish before excute bsub command for now since we wait bam2fastq finish before updating the fastq path
            cmd = "bsub -J {}_{}_multi -o {}_{}_multi.out{}--id={}_{} --csv={}{}".format(self.name, key, self.name, key, config_dict["multi"]["tool"], self.name, key, name_of_file, OPTIONS)
            print(cmd)
            subprocess.run(cmd, shell=True)

    # get reads number and sub sample cell number
    def update_info_from_step1(self, fb_project_id):
        # get total reads number for gene expression library
        reads_file = "/igo/stats/PIPELINE/Project_{}_step1/{}/outs/per_sample_outs/{}/metrics_summary.csv".format(fb_project_id, self.name, list(self.samples.keys())[0])
        summary_metrix = pd.read_csv(reads_file)
        ind = summary_metrix.index[(summary_metrix["Category"] == "Library") & (summary_metrix["Metric Name"] == "Number of reads") & (summary_metrix["Library Type"] == "Gene Expression") & (summary_metrix["Grouped By"] == "Physical library ID")].tolist()
        reads_number = summary_metrix.iloc[ind[0]]["Metric Value"]
        reads_number = int(reads_number.replace(",", "")) + 10000
        self.ge_reads_number = reads_number

        # update sub sample cell number
        cell_file = "/igo/stats/PIPELINE/Project_{}_step1/{}/outs/multi/multiplexing_analysis/tag_calls_summary.csv".format(fb_project_id, self.name)
        cell_matrix = pd.read_csv(cell_file)
        for key, value in self.samples.items():
            if value in cell_matrix["Category"].values:
                cell_number = int(cell_matrix.iloc[cell_matrix.index[cell_matrix["Category"] == value]]["num_cells"])
                self.sub_sample_info[key] = [cell_number]

    # get fastq location after bam2fastq, need to run for each sub sample
    # doesn't know what is the pattern when have two runs for one sample. leave as it is for now
    def update_fastq_location(self, sub_sample_name, destination):
        file = glob.glob(destination + "/" + self.name + "_0_*")
        self.sub_sample_info[sub_sample_name].append(file[0])

# read ch file from shared drive and generate config/ch file per sample and return sample to sub sample info also
# default all hash tag are totalseq B from biolegend
def ch_file_generation(project_id, sample_name):
    in_file_location = DRIVE_LOCATION + project_id + "/" + os.listdir(DRIVE_LOCATION + project_id)[0]
    with open(in_file_location, "rb") as f:
        df = pd.read_excel(f, engine="openpyxl")
    line_number = df[df[df.columns[0]] == "Your Submission:"].index.values
    with open(in_file_location, "rb") as f:
        df = pd.read_excel(f, engine="openpyxl", skiprows=line_number + 1, header=line_number + 1)
    sample_tag_dict = pd.Series(df['Hashtag Name'].values,index=df['Sample Name']).to_dict()
    tag_seq_dict = pd.Series(df['Hashtag sequence'].values,index=df['Hashtag Name']).to_dict()

    sub_sample_dict = {}
    sub_sample_lst = df[df["Sample Name in IGO"] == sample_name]["Sample Name"].tolist()
    for item in sub_sample_lst:
        sub_sample_dict[item] = sample_tag_dict[item]

    # write ch config file for this sample
    file_name = "{}Project_{}/Project_{}_ch_{}.csv".format(CONFIG_AREA, project_id, project_id, sample_name)
    if not os.path.exists(os.path.dirname(file_name)):
        os.makedirs(os.path.dirname(file_name))

    with open(file_name,'w') as file:
        file.write("id,name,read,pattern,sequence,feature_type\n")
        for tag in sub_sample_dict.values():
            file.write("{},{},R2,5PNNNNNNNNNN(BC),{},Multiplexing Capture\n".format(tag, tag, tag_seq_dict[tag]))

    return(sub_sample_dict)

def gather_config_info(sample_dict, genome, IGO_ID):
    """
    sample_dict contains all the information about samples, sample name, project_ID and recipe name
    example: {"ge":"LJ01_IGO_14396_1", "vdj": "LJ01_VDJ_IGO_14396_C_1", "fb":"", "ch":""}
    one file per sample named with IGO ID (GEX sample ID)
    """
    # library use fastq id as key, then followed by a list which first item is list of fastq path in case top up and second is feature type
    # fb reference file should have format as following: /igo/stats/Multi_config/Project_12345/Project_12345_fb.csv
    # ch reference file should have format as following: /igo/stats/Multi_config/Project_12345/Project_12345_ch_ABC.csv
    # how to record vdj-t and vdj-b?
    project_ID = "_".join(IGO_ID.split("IGO_")[1].split("_")[:-1])
    sample_name = IGO_ID.split("_IGO_")[0]
    config = Multi_Config()
    config.name = IGO_ID
    config.gene_expression["reference"] = config_dict["count"]["genome"][genome][17:]
    if "vdj" in sample_dict.keys():
        config.vdj = config_dict["vdj"]["genome"][genome][13:]
    
    # if feature barcoding invovled, add feature list file path
    if "fb" in sample_dict.keys():
        config.features = CONFIG_AREA + "Project_{}/Project_{}_fb.csv".format(project_ID, project_ID)
        
    # if cell hashing invovled, add cmo-set file path and get sample info from file, id as sample name and name as hashtag name
    if "ch" in sample_dict.keys():
        config.gene_expression["cmo-set"] = CONFIG_AREA + "Project_{}/Project_{}_ch_{}.csv".format(project_ID, project_ID, sample_name)
        config.samples = ch_file_generation(project_ID, sample_name)

    # if both ch and fb are there, change the ch name
    if "ch" in sample_dict.keys() and "fb" in sample_dict.keys():
        sample_dict["ch"] = sample_dict["ch"].replace("FB_IGO", "CH_IGO")

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

# for case of ch + vdj +/- fb
def cellragner_ch_vdj(config, file_name, ch_project_ID, project_ID, ge):
    # run ch + ge first under PIPELINE folder with name of Project_fb_step1
    config.write_ch_ge_only_to_csv(file_name)
    cmd = "bsub -K -J {}_multi -o {}_multi.out{}--id={} --csv={}{}".format(ge, ge, config_dict["multi"]["tool"], ge, file_name, OPTIONS)
    os.chdir(STATS_AREA)
    projects = next(os.walk("."))[1]
    project = "Project_{}_step1".format(ch_project_ID)
    if project not in projects:
        os.mkdir(project, ACCESS)
    work_area = STATS_AREA + project + "/" 
    # GO TO project ID LOCATION to start cellranger command
    os.chdir(work_area)
    print(cmd)
    subprocess.run(cmd, shell=True)
    
    # update cell number and ge reads number after ge + ch finish
    config.update_info_from_step1(ch_project_ID)
    # create bamtofastq folder if not exists
    try:
        os.mkdir("{}Project_{}/bamtofastq".format(CONFIG_AREA, project_ID))
    except OSError as error:
        print(error)  
        
    # create bam2fastq cmd per sub sample
    for key in config.sub_sample_info.keys():
        name2 = ge + "_" + key
        source_bam = "/igo/stats/PIPELINE/Project_{}_step1/{}/outs/per_sample_outs/{}/count/sample_alignments.bam".format(ch_project_ID, ge, key)
        destination_bam = "{}Project_{}/bamtofastq/{}".format(CONFIG_AREA, project_ID, name2)
        cmd = "bsub -K -J {}_bamtofastq -o {}_bamtofastq.out -n 8 -M 8 {} --reads-per-fastq={} {} {}".format(name2, name2, BAMTOFASTQ, config.ge_reads_number, source_bam, destination_bam)
        print(cmd)
        subprocess.run(cmd, shell=True)
        # update new fastq file path after bam2fastq for each sub sample
        config.update_fastq_location(key, destination_bam)
    
    os.chdir(STATS_AREA)
    projects = next(os.walk("."))[1]
    project = "Project_{}".format(ch_project_ID)
    if project not in projects:
        os.mkdir(project, ACCESS)
    work_area = STATS_AREA + project + "/" 
    # GO TO project ID LOCATION to start cellranger command
    os.chdir(work_area)

    config.new_config_and_generate_cmd()

# for case of ch + fb - vdj
def cellranger_ch_fb(config, file_name, ch_project_ID, ge, ch, fb):
    # Process: modify fb fastq files and store in new location then proceed
    # 1. copy fb fastq to a separate location
    # 2. modify the fastq using modify_fastq_for_fb.sh
    # 3. same process as normal cases
    new_ch_sample_name = ch.replace("FB_IGO", "CH_IGO")
    DESTINATION_CH_FASTQ_prefix = "/igo/stats/Multi_config/"
    os.chdir(DESTINATION_CH_FASTQ_prefix)
    runs = next(os.walk("."))[1]
    # copy fb fastq file to /igo/stats/Multi_config/<RUN>/<Sample>
    for i in config.lirbaries[fb][0]:
        print(i)
        # /igo/staging/FASTQ/RUTH_0233_AHHYVKDSX5/Project_13422_F/Sample_19288_66_IGO_13422_F_1
        run = i.split("/")[4]
        if run not in runs:
            os.mkdir(run, ACCESS)
        cmd = "cp -R {}/ {}/".format(i, run)
        print(cmd)
        # subprocess.run(cmd, shell=True)
        sample = i.split("/")[6]
        # go to the fastq folder and modify fastq files. The name will be CH replaced FB
        DESTINATION_CH_FASTQ = "{}{}/{}".format(DESTINATION_CH_FASTQ_prefix, run, sample)
        print(DESTINATION_CH_FASTQ)
        os.chdir(DESTINATION_CH_FASTQ)
        job_name = "modify_fb_{}_{}".format(run, sample)
        cmd = "bsub -J {} -o '/igo/stats/Multi_config/{}/{}.out' -n 8 -M 8 sh /igo/work/igo/igo-demux/scripts/modify_fastq_for_fb.sh".format(job_name, run, job_name)
        print(cmd)
        # subprocess.run(cmd, shell=True)

        # append new fastq location for ch sample under config class
        config.lirbaries[new_ch_sample_name][0].append(DESTINATION_CH_FASTQ)
        # create config and submit job to wait for modify fastq finish before excute
        config.write_to_csv(file_name)
        cmd = "bsub -J {}_multi -o {}_multi.out -w \"done(*{})\"{}--id={} --csv={}{}".format(ge, ge, sample, config_dict["multi"]["tool"], ge, file_name, OPTIONS)
        # create project folder if not exists
        os.chdir(STATS_AREA)
        projects = next(os.walk("."))[1]
        project = "Project_" + ch_project_ID
        if project not in projects:
            os.mkdir(project, ACCESS)
        work_area = STATS_AREA + project + "/" 
        # GO TO project ID LOCATION to start cellranger command
        os.chdir(work_area)
        print("Start cellranger from {}".format(work_area))
        print(cmd)
        # subprocess.run(cmd, shell=True)

# for other simple cases
def cellranger_general(config, file_name, ch_project_ID, ge):
    config.write_to_csv(file_name)
    cmd = "bsub -J {}_multi -o {}_multi.out{}--id={} --csv={}{}".format(ge, ge, config_dict["multi"]["tool"], ge, file_name, OPTIONS)
    # create project folder if not exists
    os.chdir(STATS_AREA)
    projects = next(os.walk("."))[1]
    project = "Project_" + ch_project_ID
    if project not in projects:
        os.mkdir(project, ACCESS)
    work_area = STATS_AREA + project + "/" 
    # GO TO project ID LOCATION to start cellranger command
    os.chdir(work_area)
    print("Start cellranger from {}".format(work_area))
    print(cmd)
    subprocess.run(cmd, shell=True)

# parsing LIMS endpoint to get sample set for cellranger multi using gene expression sample
# example input 190121-TSR1_IGO_15041_1, expected return{"ge": "190121-TSR1_IGO_15041_1", "ch": "190121-TSR1_FB_IGO_15041_B_1"}
def gather_sample_set_info(sample_name):
    sample_set = {"ge": sample_name, "ch": None, "vdj": None, "fb": None}
    IGO_ID = sample_name.split("_IGO_")[1]
    sample_number = IGO_ID.split("_")[-1]
    project_ID = "_".join(IGO_ID.split("_")[0:-1])
    response = requests.get(ENDPOINT + project_ID , auth = ("pms", "tiagostarbuckslightbike"), verify = False)
    response_data = json.loads(response.text.encode("utf8"))
    # get ilab request info to get full set sampels later
    for sample in response_data:
        for key, value in sample.items():
            if key == IGO_ID:
                ilab_request = value[0]
                # TODO not sure what is the return value if both are there
                fb_type = value[2].split(",")[1].strip()
                vdj_type = value[2].split(",")[2].strip()
                print(fb_type, vdj_type)
                break

    # using ilab and sample name for this 
    for sample in response_data:
        for key, value in sample.items():
            if value[0].startswith(ilab_request) and key.endswith(sample_number):
                value[2] = value[2].split(",")
                if "10X_Genomics_FeatureBarcoding" in value[2][0]:
                    if fb_type == "Feature Barcoding":
                        sample_set["fb"] = "_IGO_".join([value[1], key])
                    elif fb_type == "Cell Hashing":
                        sample_set["ch"] = "_IGO_".join([value[1], key])
                if "10X_Genomics_VDJ" in value[2][0]:
                    sample_set["vdj"] = "_IGO_".join([value[1], key])

    return sample_set

# TODO check whether a project set is complete to launch pipeline

# TODO fb file generation from user form

if __name__ == '__main__':
    # input as name for each library type plus genome
    # Usage: python cellranger_multi.py -ge=AT3_C1-hashtag_IGO_14767_1 -ch=AT3_C1-hashtag_FB_IGO_14767_B_1 -genome=Mouse
    parser = argparse.ArgumentParser(prog = 'Cellranger Multi', usage = 'Run the pipeline for cellranger multi')
    parser.add_argument('-ge', required = True)
    parser.add_argument('-vdj')
    parser.add_argument('-ch')
    parser.add_argument('-fb')
    parser.add_argument('-genome', help = 'Human or Mouse', required = True)
    args = parser.parse_args()
    sample_dict = {"ge":args.ge}
    if args.vdj:
        sample_dict["vdj"] = args.vdj
    if args.ch:
        sample_dict["ch"] = args.ch
        ch_project_ID = "_".join(args.ch.split("IGO_")[1].split("_")[:-1])

    if args.fb:
        sample_dict["fb"] = args.fb
        ch_project_ID = "_".join(args.fb.split("IGO_")[1].split("_")[:-1])
    
    genome = args.genome
    config = gather_config_info(sample_dict, genome, args.ge)
    project_ID = "_".join(args.ge.split("IGO_")[1].split("_")[:-1])
    file_name = "{}Project_{}/{}.csv".format(CONFIG_AREA, project_ID, args.ge)

    # condition for ch + vdj +/- fb
    if args.ch and args.vdj:
        cellragner_ch_vdj(config, file_name, ch_project_ID, project_ID, args.ge)
    # condition for ch + fb - vdj
    elif args.ch and args.fb:
        cellranger_ch_fb(config, file_name, ch_project_ID, args.ge, args.ch, args.fb)
    
    # other normal cases
    else:
        cellranger_general(config, file_name, ch_project_ID, args.ge)
