# read in csv stats files generated by dragen alignment and create txt files that can be parsed into NGS database
# for WGS now
import pandas as pd
import sys
import os

# input will be the path to all the csv files and output path for txt files

# read in all the names of csv files, get the name before first dot as sampleID, save in a list called sample_list
def get_sample_list(folder_path):
    sample_list = []
    # get all the file name in the stats folder
    file_list = os.listdir(folder_path)
    # if the format of the file name is <sample_name>.<category>.csv, add uniqle sample name to sample_list
    for i in file_list:
        ls = i.split('.')
        if len(ls) == 3 and (ls[-1] == "csv") and (ls[0] not in sample_list):
            sample_list.append(ls[0])
    
    return sample_list


# for each sample in the list, get info needed from csv files and generate txt files
# txt file needed: ___AM.txt, ___WGS.txt, ___MD.txt
class DragenStats:
    recipe = "WGS"
    sep = "___"

    def __init__(self, sample_name):
        self.sample_name = sample_name

    def read_info_from_csv(self, folder_path):
        # open mapping_metrics csv and store info in dataframe
        mapping_file_name = folder_path + self.sample_name + ".mapping_metrics.csv"
        df_mapping_metrics = pd.read_csv(mapping_file_name, index_col=2, nrows=53, names=[0,1,2,3,4])
        self.READ_PAIRS_EXAMINED = int (df_mapping_metrics.loc["Properly paired reads"][3] / 2)
        self.UNMAPPED_READS = int (df_mapping_metrics.loc["Unmapped reads"][3])
        self.TOTAL_READS = int (df_mapping_metrics.loc["Total input reads"][3])
        self.READ_PAIR_DUPLICATES = int (df_mapping_metrics.loc["Number of duplicate marked reads"][3] / 2)
        self.PERCENT_DUPLICATION = df_mapping_metrics.loc["Number of duplicate marked reads"][4] / 100
        self.READS_ALIGNED_IN_PAIRS = int (df_mapping_metrics.loc["Properly paired reads"][3])
        # open coverage_metrics*.csv and store info in dataframe
        coverage_file_name = folder_path + self.sample_name + ".wgs_coverage_metrics.csv"
        df_coverage_metrics = pd.read_csv(coverage_file_name, index_col=2, names=[0,1,2,3,4])
        self.MEAN_TARGET_COVERAGE = df_coverage_metrics.loc["Average alignment coverage over genome"][3]
        self.PF_READS_ALIGNED = int(df_coverage_metrics.loc["Aligned reads"][3])


    # creat ___AM.txt picard stats format file
    # requirments: 7 lines minimun, line with data need to start with Category PAIR, file name ends with ___AM.txt
    def write_to_am_txt(self, output_path):
        data_list_to_write = [0] * 24
        data_list_to_write[0] = "PAIR"
        data_list_to_write[1] = self.TOTAL_READS
        data_list_to_write[5] = self.PF_READS_ALIGNED
        data_list_to_write[16] = self.READS_ALIGNED_IN_PAIRS

        write_to_file = output_path + self.sample_name + "___DRAGEN3_10_8___AM.txt"
        data_line = ""
        for i in data_list_to_write:
            data_line = data_line + str(i) + "\t"
        
        with open(write_to_file, 'w') as _file:
            _file.write("#" + self.sample_name + "\n")
            for i in range(6):
                _file.write("#\n")
            _file.write(data_line)

    # creat ___MD.txt picard stats format file

    def write_to_md_txt(self, output_path):
        data_list_to_write = [0] * 10
        data_list_to_write[0] = self.sample_name
        data_list_to_write[2] = self.READ_PAIRS_EXAMINED
        data_list_to_write[4] = self.UNMAPPED_READS
        data_list_to_write[6] = self.READ_PAIR_DUPLICATES
        data_list_to_write[8] = self.PERCENT_DUPLICATION
        
        header = "LIBRARY	UNPAIRED_READS_EXAMINED	READ_PAIRS_EXAMINED	SECONDARY_OR_SUPPLEMENTARY_RDS	UNMAPPED_READS	UNPAIRED_READ_DUPLICATES	READ_PAIR_DUPLICATES	READ_PAIR_OPTICAL_DUPLICATES	PERCENT_DUPLICATION	ESTIMATED_LIBRARY_SIZE"
        write_to_file = output_path + self.sample_name + "___DRAGEN3_10_8___MD.txt"
        data_line = ""
        for i in data_list_to_write:
            data_line = data_line + str(i) + "\t"
        
        with open(write_to_file, 'w') as _file:
            _file.write("#" + self.sample_name + "\n" + "#\n" + "\n")
            _file.write(header + "\n")
            _file.write(data_line)
           
    # creat ___WGS.txt picard stats format file
    def write_to_wgs_txt(self, output_path):
        data_list_to_write = [0] * 28
        data_list_to_write[1] = self.MEAN_TARGET_COVERAGE
        
        header = "GENOME_TERRITORY	MEAN_COVERAGE	SD_COVERAGE	MEDIAN_COVERAGE	MAD_COVERAGE	PCT_EXC_ADAPTER	PCT_EXC_MAPQ	PCT_EXC_DUPE	PCT_EXC_UNPAIRED	PCT_EXC_BASEQ	PCT_EXC_OVERLAP	PCT_EXC_CAPPED	PCT_EXC_TOTAL	PCT_1X	PCT_5X	PCT_10X	PCT_15X	PCT_20X	PCT_25X	PCT_30X	PCT_40X	PCT_50X	PCT_60X	PCT_70X	PCT_80X	PCT_90X	PCT_100X	FOLD_80_BASE_PENALTY	FOLD_90_BASE_PENALTY	FOLD_95_BASE_PENALTY	HET_SNP_SENSITIVITY	HET_SNP_Q"
        write_to_file = output_path + self.sample_name + "___DRAGEN3_10_8___WGS.txt"
        data_line = ""
        for i in data_list_to_write:
            data_line = data_line + str(i) + "\t"
        
        with open(write_to_file, 'w') as _file:
            _file.write("#" + self.sample_name + "\n")
            for i in range(4):
                _file.write("#\n")
            _file.write("\n" + header + "\n")
            _file.write(data_line)

if __name__ == "__main__":
    # Usage: python dragenstats_csv_to_txt.py [dragen_stats_dir] [output_file_dir]
    # example: python3 /Users/luc/Documents/GitHub/igo-demux/scripts/dragenstats_csv_to_txt.py /Users/luc/Documents/GitHub/igo-demux/test/ /Users/luc/Documents/GitHub/igo-demux/test/result_test/
    dragen_stats_folder = sys.argv[1]
    output_folder_path = sys.argv[2]
    sample_list = get_sample_list(dragen_stats_folder)
    for i in sample_list:
        dragen_stats = DragenStats(i)
        dragen_stats.read_info_from_csv(dragen_stats_folder)
        dragen_stats.write_to_am_txt(output_folder_path)
        dragen_stats.write_to_md_txt(output_folder_path)
        dragen_stats.write_to_wgs_txt(output_folder_path)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
