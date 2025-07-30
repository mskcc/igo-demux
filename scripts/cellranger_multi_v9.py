import pandas as pd
import os
import subprocess
import glob
import argparse
from collections import OrderedDict
import requests
import json
import openpyxl
import scripts.cellranger_config as CONFIG
from scripts.cellranger import find_fastq_file

# config file class. It contains all the information needed for config file and can create csv file base on those info
class Multi_Config:
    def __init__(self, sample_set, genome, archive):
        self.name = "EMPTY"  # gene expression sample name from fastq file
        self.username = ""  # sample name user submitted
        self.genome = ""   # genome info for the ge sample
        self.lirbaries = OrderedDict()   # location info for fastq files
        self.vdj = "EMPTY"    # genome info for the vdj sample
        self.samples = "EMPTY"  # sub_samples for cell hashing
        self.features = "EMPTY"   # reference file created based on user provided information
        self.ch_project_id = "" # fb project id such as 12345, used for pipeline folder name
        self.ge_project_id = "" # ge project id such as 12345
        self.ch = False
        self.fb = False
        self.gather_config_info(sample_set, genome, archive)

    def gather_config_info(self, sample_set, genome, archive):
        """
        sample_set contains all the information about samples, sample name, project_ID and recipe name
        example: {"ge":"LJ01_IGO_14396_1", "vdj_t": "LJ01_VDJ_IGO_14396_C_1", "fb":"", "ch":""}
        one file per sample named with IGO ID (GEX sample ID)
        """
        # library use fastq id as key, then followed by a list which first item is list of fastq path in case top up and second is feature type
        # fb/ch reference file should have format as following: /igo/stats/Multi_config/Project_12345/Project_12345_reference_sample_name.csv
        self.ge_project_id = "_".join(sample_set["ge"].split("IGO_")[1].split("_")[:-1])
        self.username = sample_set["ge"].split("_IGO_")[0]
        self.features = "{}Project_{}/Project{}_reference_{}.csv".format(CONFIG.CONFIG_AREA, self.ge_project_id, self.ge_project_id, self.username)
        self.name = sample_set["ge"]
        self.genome = CONFIG.config_dict["count"]["genome"][genome][17:]
        if "vdj-t" in sample_set or "vdj-b" in sample_set:
            self.vdj = CONFIG.config_dict["vdj"]["genome"][genome][13:]
        
        # find fastq files for each sample and append information into config["libraries"]
        sample_list = []
        for i in sample_set.values():
            sample_list.append(i)
        fastq_list = find_fastq_file(sample_list, archive)
        for key, value in sample_set.items():
            print("key: {}, value: {}".format(key, value))
            if key == "ge":
                self.lirbaries[value] = [fastq_list[value], "Gene Expression"]
            elif key == "vdj_t":
                self.lirbaries[value] = [fastq_list[value], "VDJ-T"]
            elif key == "vdj_b":
                self.lirbaries[value] = [fastq_list[value], "VDJ-B"]
            elif key == "fb":
                self.lirbaries[value] = [fastq_list[value], "Antibody Capture"]
                self.ch_project_id = "_".join(value.split("IGO_")[1].split("_")[:-1])
                self.fb = True
            elif key =="ch":
                self.lirbaries[value] = [fastq_list[value], "Antibody Capture"]
                self.ch_project_id = "_".join(value.split("IGO_")[1].split("_")[:-1])
                self.ch = True

    # read ch/fb file from shared drive and generate config/ch file per sample and return sample to sub sample info
    # default all hash tag are totalseq B from biolegend
    def reference_file_generation(self):
        # ch file
        if self.ch:
            in_file_location = glob.glob("{}{}/*cell_hash.xlsx".format(CONFIG.DRIVE_LOCATION, self.ge_project_id))[0]
            with open(in_file_location, "rb") as f:
                df_all = pd.read_excel(f, engine="openpyxl")
            line_number = df_all[df_all[df_all.columns[0]] == "Your Submission:"].index.values[0]
            print("starting from line {}".format(line_number))
            new_header = df_all.iloc[line_number + 1]
            df = df_all.iloc[line_number + 2:].copy()
            df.columns = new_header
            df = df.reset_index(drop=True)

            # update column name for new template
            df.columns = ['Sample Name' if col == 'Sample Name Pre-Hashing' else col for col in df.columns]
            df.columns = ['Sample Name in IGO' if col == 'Sample ID from Sample Submission' else col for col in df.columns]

            sample_tag_dict = pd.Series(df['Hashtag Name'].values,index=df['Sample Name']).to_dict()
            tag_seq_dict = pd.Series(df['Hashtag sequence'].values,index=df['Hashtag Name']).to_dict()

            self.samples = {}
            sub_sample_lst = df[df["Sample Name in IGO"].astype(str) == str(self.username)]["Sample Name"].tolist()
            for item in sub_sample_lst:
                self.samples[item] = sample_tag_dict[item]
        
        # fb file
        if self.fb:      
            in_file_location = glob.glob("{}{}/*feature_barcoding.xlsx".format(CONFIG.DRIVE_LOCATION, self.ge_project_id))[0]  
            with open(in_file_location, "rb") as f:
                df_all = pd.read_excel(f, engine="openpyxl")
            line_number = df_all[df_all[df_all.columns[0]] == "Your Submission:"].index.values[0]
            print("starting from line {}".format(line_number))
            new_header = df_all.iloc[line_number + 1]
            df = df_all.iloc[line_number + 2:].copy()
            df.columns = new_header
            df = df.reset_index(drop=True)
            antibody_seq_dict = pd.Series(df['Sequence'].values,index=df['Target']).to_dict()

        # write reference file for this sample
        try:
            os.makedirs(os.path.dirname(self.features))
        except OSError as error:
            print(error) 
        with open(self.features,'w') as file:
            file.write("id,name,read,pattern,sequence,feature_type\n")
            if self.fb:
                for antibody in antibody_seq_dict.keys():
                    file.write("{},{},R2,5PNNNNNNNNNN(BC),{},Antibody Capture\n".format(antibody, antibody, antibody_seq_dict[antibody]))
            if self.ch:
                for tag in self.samples.values():
                    file.write("{},{},R2,5PNNNNNNNNNN(BC),{},Antibody Capture\n".format(tag, tag, tag_seq_dict[tag]))                    

    # write all info from the class to csv file
    def write_to_csv(self):
        name_of_file = "{}Project_{}/{}.csv".format(CONFIG.CONFIG_AREA, self.ge_project_id, self.name)
        with open(name_of_file,'w') as file:
            file.write("[gene-expression]\n")
            file.write("reference,{}\n".format(self.genome))
            file.write("create-bam,true\n")

            file.write("\n[libraries]\nfastq_id,fastqs,feature_types\n")
            
            for key, value in self.lirbaries.items():
                for i in value[0]:
                    file.write("{},{},{}\n".format(key, i, value[1]))
            if self.vdj != "EMPTY":
                file.write("\n[vdj]\nreference,{}\n".format(self.vdj))

            if self.features != "EMPTY":
                file.write("\n[feature]\nreference,{}\n".format(self.features))

            if self.samples != "EMPTY":
                file.write("\n[samples]\nsample_id,hashtag_ids\n")
                for key, value in self.samples.items():
                    file.write("{},{}\n".format(key, value))
        return(name_of_file)


# parsing LIMS endpoint to get sample set for cellranger multi using gene expression sample
# example input 190121-TSR1_IGO_15041_1, expected return{"ge": "190121-TSR1_IGO_15041_1", "ch": "190121-TSR1_FB_IGO_15041_B_1"}
def gather_sample_set_info(sample_name):
    sample_set = {"ge": sample_name, "ch": None, "vdj_t": None, "vdj_b": None, "fb": None}
    IGO_ID = sample_name.split("_IGO_")[1]
    sample_number = IGO_ID.split("_")[-1]
    project_ID = "_".join(IGO_ID.split("_")[0:-1])
    response = requests.get(CONFIG.MULTI_ENDPOINT + project_ID , auth = ("pms", "tiagostarbuckslightbike"), verify = False)
    response_data = json.loads(response.text.encode("utf8"))
    # get ilab request info to get full set sampels later
    for sample in response_data:
        for key, value in sample.items():
            if key == IGO_ID:
                ilab_request = value[0]
                # The return value is whole list, order of Cell Hashing, Feature Barcoding, T Cells, B Cells
                tag_lst = [x.strip() for x in value[2].split(',')]
                print(tag_lst)
                fb_type = []
                if "Cell Hashing" in tag_lst:
                    fb_type.append("Cell Hashing")
                if "Feature Barcoding" in tag_lst:
                    fb_type.append("Feature Barcoding")
                print(fb_type)
                break

    # using ilab and sample name for this 
    for sample in response_data:
        for key, value in sample.items():
            if value[0].startswith(ilab_request[:-1]) and key.endswith(sample_number):
                value[2] = value[2].split(",")
                if "SC_Chromium-FB-5" in value[2][0] or "SC_Chromium-FB-3" in value[2][0]:
                    if "Feature Barcoding" in fb_type:
                        sample_set["fb"] = "_IGO_".join([value[1], key])
                    if "Cell Hashing" in fb_type:
                        sample_set["ch"] = "_IGO_".join([value[1], key])
                if "SC_Chromium-BCR" in value[2][0]:
                    sample_set["vdj_b"] = "_IGO_".join([value[1], key])
                if "SC_Chromium-TCR" in value[2][0]:
                    sample_set["vdj_t"] = "_IGO_".join([value[1], key])

    return sample_set

 # check if the template is fb or ch by reading the first line   
def check_file_type(file_path):
    try:
        workbook = openpyxl.load_workbook(file_path, read_only=True)
        sheet = workbook.active
        first_row = [str(cell.value).strip() for cell in next(sheet.iter_rows(max_row=1)) if cell.value is not None]
        # Check for "cell hashing" or "feature barcoding"
        if any("Cell Hashing" in cell for cell in first_row):
            return "ch"
        elif any("Feature Barcoding" in cell for cell in first_row):
            return "fb"
        else:
            return None                
    except Exception as e:
        print(f"An error occurred: {e} while processing file {file_path}")

# TODO check whether a project set is complete to launch pipeline

# main process to start the pipeline 
def launch_pipeline_by_sample(sample_set, genome, archive):
    sample = Multi_Config(sample_set, genome, archive)
    os.chdir(CONFIG.MULTI_STATS_AREA)
    # create project folder if not exists
    project = "Project_" + sample.ch_project_id
    try:
        os.mkdir(project, CONFIG.ACCESS)
    except OSError as error:
        print(error)
    sample.reference_file_generation()
    file_name = sample.write_to_csv()

    cmd = "bsub -J {}_multi -o {}_multi.out{}--id={} --csv={}{}".format(sample.name, sample.name, CONFIG.config_dict["multi"]["tool"], sample.name, file_name, CONFIG.MULTI_OPTIONS)

    work_area = CONFIG.MULTI_STATS_AREA + project + "/" 
    # GO TO project ID LOCATION to start cellranger command
    os.chdir(work_area)
    print("Start cellranger from {}".format(work_area))
    print(cmd)
    # subprocess.run(cmd, shell=True)

def launch_multi_by_project_location(project_directory, genome):
    project_id = project_directory.split("/")[-1]
    # copy the multi config from shared drive to cluster
    cmd = "cp -R {}{} {}".format(CONFIG.ORIGIN_DRIVE_LOCATION, project_id[8:], CONFIG.DRIVE_LOCATION)
    print(cmd)
    subprocess.run(cmd, shell=True)
    # add file checking for ch and fb because lims request info not accurate at this moment
    file_lst = []
    file_prefix = CONFIG.DRIVE_LOCATION + project_id[8:]
    file_lst = os.listdir(file_prefix)
    print(file_lst)
    ch = False
    fb = False
    for i in file_lst:
        file_path = file_prefix + "/" + i
        file_type = check_file_type(file_path)
        if file_type == "ch":
            ch = True
            ch_file_name = "{}/{}_cell_hash.xlsx".format(file_prefix, project_id)
            os.rename(file_path, ch_file_name)
            print(f"File renamed from {file_path} to {ch_file_name}")
        elif file_type == "fb":
            fb = True
            fb_file_name = "{}/{}_feature_barcoding.xlsx".format(file_prefix, project_id)
            os.rename(file_path, fb_file_name)
            print(f"File renamed from {file_path} to {fb_file_name}")

    # gather sample set info from LIMS for each sample
    archive = False
    if "delivery" in project_directory:
        archive = True

    sample_list_ori = os.listdir(project_directory)
    sample_list = []
    # remove Sample_ prefix
    for sample in sample_list_ori:
        sample_list.append(sample[7:])
    for sample in sample_list:
        sample_set = gather_sample_set_info(sample)
        # update sample_set based on file checking result
        if sample_set["ch"] is not None:
            sample_name = sample_set["ch"]
        elif sample_set["fb"] is not None:
            sample_name = sample_set["fb"]
        del sample_set["fb"]
        del sample_set["ch"]
        if ch:
            sample_set["ch"] = sample_name
        if fb:
            sample_set["fb"] = sample_name
        
        # remove None value
        sample_set = {k: v for k, v in sample_set.items() if v is not None}

        print(sample_set)
        launch_pipeline_by_sample(sample_set, genome, archive)

if __name__ == '__main__':
    # input as name for each library type plus genome
    # Usage: python cellranger_multi.py -ge_path=/igo/staging/FASTQ/RUTH_0141_AH27NGDSX5/Project_13586_B -genome=Mouse
    parser = argparse.ArgumentParser(prog = 'Cellranger Multi', usage = 'Run the pipeline for cellranger multi')
    parser.add_argument('-ge_path', required = True)
    parser.add_argument('-genome', help = 'Human or Mouse', required = True)
    args = parser.parse_args()
    
    launch_multi_by_project_location(args.ge_path, args.genome)
