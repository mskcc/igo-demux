import pandas as pd
import sys
import numpy

# get total reads number from Dragen Demultiplex_Stats.csv file and generate txt files for all the samples on samplesheet

"""
main function: run
input: samplesheet object, demux directory
ouput: txt files containing total reads and ready for upload to website
"""

def get_total_reads(sample_ID_list, demux_report_file):
    demux_df = pd.read_csv(demux_report_file, index_col=1)
    total_reads_dict = {}
    for sample_ID in sample_ID_list:
        if isinstance(demux_df.loc[sample_ID]["# Reads"], numpy.int64):
            total_reads_dict[sample_ID] = demux_df.loc[sample_ID]["# Reads"] * 2
        else:
            total_reads_dict[sample_ID] = sum(demux_df.loc[sample_ID]["# Reads"]) * 2
    
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
    sequencer_and_run_prefix = "_".join(sequencer_and_run.split("_")[0:3])
    sequencer = sequencer_and_run.split("_")[0]
    stats_done_dir = "/igo/stats/DONE/" + sequencer + "/"
    demux_report_file = "/igo/staging/FASTQ/" + sequencer_and_run + "/Reports/Demultiplex_Stats.csv"
    # dictionary of Sample_ID->Project
    sample_project_dict = pd.Series(sample_sheet.df_ss_data['Sample_Project'].values,index=sample_sheet.df_ss_data['Sample_ID']).to_dict()
    sample_ID_list = list(sample_project_dict.keys())
    total_reads_dict = get_total_reads(sample_ID_list, demux_report_file)
    for sample in sample_ID_list:
        write_to_am_txt(sequencer_and_run_prefix, sample, total_reads_dict[sample], stats_done_dir)

def test_run(sample_sheet, sequencer_and_run):
    sequencer_and_run_prefix = "_".join(sequencer_and_run.split("_")[0:3])
    sequencer = sequencer_and_run.split("_")[0]
    stats_done_dir = "/Users/luc/Documents/GitHub/igo-demux/test/result_test/"
    demux_report_file = "/Users/luc/Documents/GitHub/igo-demux/test/result_test/Demultiplex_Stats.csv"
    # dictionary of Sample_ID->Project
    sample_project_dict = pd.Series(sample_sheet.df_ss_data['Sample_Project'].values,index=sample_sheet.df_ss_data['Sample_ID']).to_dict()
    sample_ID_list = list(sample_project_dict.keys())
    total_reads_dict = get_total_reads(sample_ID_list, demux_report_file)
    for sample in sample_ID_list:
        write_to_am_txt(sequencer_and_run_prefix, sample, total_reads_dict[sample], stats_done_dir)