import pandas as pd
import sys
import numpy
import json
import re

# get total reads number from Demultiplex_Stats.csv file or json file and generate txt files for each sample
# add DLP type function. For DLP, only total reads for each project is needed

"""
input: samplesheet object, demux directory
ouput: txt files containing total reads and ready for upload to website
"""

STATS_DONE_DIR_PREFIX = "/igo/stats/DONE/"

# return dictiona of sampleID -> total reads from demux report and sample_ID_list
def get_total_reads(sample_ID_list, demux_report_file):
    demux_df = pd.read_csv(demux_report_file, index_col=1)
    total_reads_dict = {}
    for sample_ID in sample_ID_list:
        if isinstance(demux_df.loc[sample_ID]["# Reads"], numpy.int64):
            total_reads_dict[sample_ID] = demux_df.loc[sample_ID]["# Reads"] * 2
        else:
            total_reads_dict[sample_ID] = sum(demux_df.loc[sample_ID]["# Reads"]) * 2
    
    return total_reads_dict

# get total reads for DLP projects. return total reads for each project from demux report and sample sheet
def get_total_reads_DLP(sample_sheet, demux_report_file):
    total_reads_dict = {}

    sample_project_dict = pd.Series(sample_sheet.df_ss_data['Sample_Project'].values,index=sample_sheet.df_ss_data['Sample_ID']).to_dict()
    # dictionary of project->sample_ID for each project, separate the samples into three groups: sample, neg control and pos control
    project_sample_dict = {}
    for sample_ID, project_ID in sample_project_dict.items():
        # if project_ID key not exist, create key first
        if project_ID not in project_sample_dict.keys():
            project_sample_dict[project_ID] = {"samples":[], "pos_control":[], "neg_control":[]}
        
        # check which category the sample belongs to and put in corresponding group
        if "DLPNegativeCONTROL" in sample_ID:
            project_sample_dict[project_ID]["neg_control"].append(sample_ID)
        elif "DLPGmCONTROL" in sample_ID:
            project_sample_dict[project_ID]["pos_control"].append(sample_ID)
        else:
            project_sample_dict[project_ID]["samples"].append(sample_ID)

    demux_df = pd.read_csv(demux_report_file, index_col=1)

    # for each project, sum up all the reads for each category as total reads and assign random sample name for each category
    for project in project_sample_dict.keys():
        total_reads_dict[project] = {"samples":["empty",0], "pos_control":["empty",0], "neg_control":["empty",0]}
        
        for category in project_sample_dict[project].keys():
            total_reads_dict[project][category][0] = project_sample_dict[project][category][0]
            for sample_ID in project_sample_dict[project][category]:
                if isinstance(demux_df.loc[sample_ID]["# Reads"], numpy.int64):
                    total_reads_dict[project][category][1] += demux_df.loc[sample_ID]["# Reads"] * 2
                else:
                    total_reads_dict[project][category][1] = sum(demux_df.loc[sample_ID]["# Reads"]) * 2
            
    return total_reads_dict

# generete AM txt file containing total reads info
def write_to_am_txt(run_ID, sample_ID, total_reads, output_path):
    data_list_to_write = [0] * 24
    data_list_to_write[0] = "PAIR"
    data_list_to_write[1] = total_reads
    data_list_to_write[5] = total_reads
    
    lst = sample_ID.split("_")
    write_to_file = output_path + run_ID + "___P" + "_".join(lst[lst.index("IGO")+1:-1]) + "___" + sample_ID + "___grch38___AM.txt"
    data_line = ""
    for i in data_list_to_write:
        data_line = data_line + str(i) + "\t"
    
    with open(write_to_file, 'w') as _file:
        _file.write("#" + sample_ID + "\n")
        for i in range(6):
            _file.write("#\n")
        _file.write(data_line)
            
def run(sample_sheet, sequencer_and_run):
    # remove postfix if existing, eg: DIANA_0502_BHM3VKDSX3_10X will convert to DIANA_0502_BHM3VKDSX3
    sequencer_and_run_prefix = "_".join(sequencer_and_run.split("_")[0:3])
    sequencer = sequencer_and_run.split("_")[0]
    stats_done_dir = STATS_DONE_DIR_PREFIX + sequencer + "/"
    demux_report_file = "/igo/staging/FASTQ/" + sequencer_and_run + "/Reports/Demultiplex_Stats.csv"
    # dictionary of Sample_ID->Project
    sample_project_dict = pd.Series(sample_sheet.df_ss_data['Sample_Project'].values,index=sample_sheet.df_ss_data['Sample_ID']).to_dict()
    sample_ID_list = list(sample_project_dict.keys())
    total_reads_dict = get_total_reads(sample_ID_list, demux_report_file)
    for sample in sample_ID_list:
        write_to_am_txt(sequencer_and_run_prefix, sample, total_reads_dict[sample], stats_done_dir)
    
    print("generate AM txt files to folder: {}".format(stats_done_dir))

# generate AM txt files containing total reads by project ID such as "Project_12754_E"
def by_project(sample_sheet, project_id, sequencer_and_run):
    sequencer_and_run_prefix = "_".join(sequencer_and_run.split("_")[0:3])
    sequencer = sequencer_and_run.split("_")[0]
    stats_done_dir = STATS_DONE_DIR_PREFIX + sequencer + "/"
    demux_report_file = "/igo/staging/FASTQ/" + sequencer_and_run + "/Reports/Demultiplex_Stats.csv"
    # dictionary of Sample_ID->Project
    sample_project_dict = pd.Series(sample_sheet.df_ss_data['Sample_Project'].values,index=sample_sheet.df_ss_data['Sample_ID']).to_dict()
    
    sample_ID_list = []
    # filter sample_ID by projectID and append to sample_ID_list
    for sample, project in sample_project_dict.items():
        if project == project_id:
            sample_ID_list.append(sample)

    total_reads_dict = get_total_reads(sample_ID_list, demux_report_file)
    for sample in sample_ID_list:
        write_to_am_txt(sequencer_and_run_prefix, sample, total_reads_dict[sample], stats_done_dir)

    print("generate AM txt files to folder: {}".format(stats_done_dir))

def by_json(sequencer_and_run):
    # remove postfix if existing, eg: DIANA_0502_BHM3VKDSX3_10X will convert to DIANA_0502_BHM3VKDSX3
    sequencer_and_run_prefix = "_".join(sequencer_and_run.split("_")[0:3])
    sequencer = sequencer_and_run.split("_")[0]
    stats_done_dir = STATS_DONE_DIR_PREFIX + sequencer + "/"
    demux_json_file = "/igo/staging/FASTQ/" + sequencer_and_run + "/Stats/Stats.json"
    # Opening JSON file
    f = open(demux_json_file)
    data = json.load(f)
    f.close()
    sample_reads_dict = {}

    for item in data["ConversionResults"]:
        for sample in item["DemuxResults"]:
            if sample["SampleName"] not in sample_reads_dict.keys():
                sample_reads_dict[sample["SampleName"]] = sample["NumberReads"] * 2
            else:
                sample_reads_dict[sample["SampleName"]] += sample["NumberReads"] * 2

    # generate AM txt files
    for sample in sample_reads_dict.keys():
        if "Sample_" in sample:
            write_to_am_txt(sequencer_and_run_prefix, sample[7:], sample_reads_dict[sample], stats_done_dir)
        else:
            write_to_am_txt(sequencer_and_run_prefix, sample, sample_reads_dict[sample], stats_done_dir)

# crete txt files of summarize total reads for DLP projects
def run_DLP(sample_sheet, sequencer_and_run):
    # remove postfix if existing, eg: DIANA_0502_BHM3VKDSX3_DLP will convert to DIANA_0502_BHM3VKDSX3
    sequencer_and_run_prefix = "_".join(sequencer_and_run.split("_")[0:3])
    sequencer = sequencer_and_run.split("_")[0]
    stats_done_dir = STATS_DONE_DIR_PREFIX + sequencer + "/"
    demux_report_file = "/igo/staging/FASTQ/" + sequencer_and_run + "/Reports/Demultiplex_Stats.csv"

    total_reads_dict = get_total_reads(sample_sheet, demux_report_file)
    
    for project in total_reads_dict.keys():
        for category in total_reads_dict[project].keys():
            data_list_to_write = [0] * 24
            data_list_to_write[0] = "PAIR"
            data_list_to_write[1] = total_reads_dict[project][category][1]
            data_list_to_write[5] = total_reads_dict[project][category][1]
                
            write_to_file = stats_done_dir + run_ID + "___P" + project[8:] + "___" + total_reads_dict[project][category][0] + "___grch38___AM.txt"
            data_line = ""
            for i in data_list_to_write:
                data_line = data_line + str(i) + "\t"
                
            with open(write_to_file, 'w') as _file:
                _file.write("#" + total_reads_dict[project][category][0] + "\n")
                for i in range(6):
                    _file.write("#\n")
                _file.write(data_line)
  
    print("generate AM txt files to folder: {}".format(stats_done_dir))


if __name__ == '__main__':
    # generate txt files with total reads info from bclconvert(CSV)/bcl2fastq(JSON) demux
    # Usage: python get_total_reads_from_demux.py [JSON/CSV][sequencer_and_run]
    run_type = sys.argv[1]
    sequencer_and_run = sys.argv[2]
    if run_type == "CSV":
        sequencer_and_run_prefix = "_".join(sequencer_and_run.split("_")[0:3])
        sequencer = sequencer_and_run.split("_")[0]
        stats_done_dir = STATS_DONE_DIR_PREFIX + sequencer + "/"
        demux_report_file = "/igo/staging/FASTQ/" + sequencer_and_run + "/Reports/Demultiplex_Stats.csv"
        sample_info_file = "/igo/staging/FASTQ/" + sequencer_and_run + "/Reports/fastq_list.csv"
            
        # get sample_list from sample_info_file
        sample_info_full = pd.read_csv(sample_info_file)
        sample_ID_list = list(set(sample_info_full['RGSM'].tolist()))
        total_reads_dict = get_total_reads(sample_ID_list, demux_report_file)
        for sample in sample_ID_list:
            write_to_am_txt(sequencer_and_run_prefix, sample, total_reads_dict[sample], stats_done_dir)
    elif run_type == "JSON":
        by_json(sequencer_and_run)
    else:
        print("run type has to be either CSV or JSON")
