# read in DRAGEN .csv stats files generated and create txt files that can be parsed into NGS database based on the Picard stats format
import pandas as pd
import sys
import os
import glob

# read in all the names of csv files, get the name before first dot as sampleID, save in a list called sample_list
def get_sample_list(folder_path):
    sample_list = []
    # get all the file name in the stats folder
    file_list = os.listdir(folder_path)
    # if the format of the file name is <sample_name>.<category>.csv, add unique sample name to sample_list
    for i in file_list:
        ls = i.split('.')
        if len(ls) == 3 and (ls[-1] == "csv") and (ls[0] not in sample_list):
            sample_list.append(ls[0])
    return sample_list


# for each sample in the list, get the info needed from .csv files and generate the .txt files
# txt files needed: ___AM.txt, ___WGS.txt, ___MD.txt
class DragenStats:

    def read_info_from_csv(self, sample_id, dragen_metrics_directory):
        print("Reading (" + sample_id + ").mapping_metrics.csv .wgs_coverage_metrics.csv")
        # open mapping_metrics csv and store info in dataframe
        mapping_metrics_file = "{}/{}.mapping_metrics.csv".format(dragen_metrics_directory, sample_id)
        df_mapping_metrics = pd.read_csv(mapping_metrics_file, index_col = 2, nrows = 53, names = [0,1,2,3,4])
        self.READ_PAIRS_EXAMINED = int (df_mapping_metrics.loc["Properly paired reads"][3] / 2)
        self.UNMAPPED_READS = int (df_mapping_metrics.loc["Unmapped reads"][3])
        self.TOTAL_READS = int (df_mapping_metrics.loc["Total input reads"][3])
        self.READ_PAIR_DUPLICATES = int (df_mapping_metrics.loc["Number of duplicate marked reads"][3] / 2)
        self.PERCENT_DUPLICATION = df_mapping_metrics.loc["Number of duplicate marked reads"][4] / 100
        self.READS_ALIGNED_IN_PAIRS = int (df_mapping_metrics.loc["Properly paired reads"][3])

        # open wgs_coverage_metrics.csv and store info in dataframe
        wgs_coverage_metrics_file = "{}/{}.wgs_coverage_metrics.csv".format(dragen_metrics_directory, sample_id)
        df_coverage_metrics = pd.read_csv(wgs_coverage_metrics_file, index_col = 2, names = [0,1,2,3,4])
        self.MEAN_TARGET_COVERAGE = df_coverage_metrics.loc["Average alignment coverage over genome"][3]
        self.PCT_1X = df_coverage_metrics.loc["PCT of genome with coverage [   1x: inf)"][3]
        self.PCT_10X = df_coverage_metrics.loc["PCT of genome with coverage [  10x: inf)"][3]
        self.PCT_15X = df_coverage_metrics.loc["PCT of genome with coverage [  15x: inf)"][3]
        self.PCT_20X = df_coverage_metrics.loc["PCT of genome with coverage [  20x: inf)"][3]
        self.PCT_50X = df_coverage_metrics.loc["PCT of genome with coverage [  50x: inf)"][3]
        self.PCT_100X = df_coverage_metrics.loc["PCT of genome with coverage [ 100x: inf)"][3]
        self.MEDIAN_COVERAGE = df_coverage_metrics.loc["Median autosomal coverage over genome"][3]
        self.PF_READS_ALIGNED = int(df_coverage_metrics.loc["Aligned reads"][3])


    # create ___AM.txt picard stats format file
    # requirements: 7 lines minimun, line with data need to start with Category PAIR, file name ends with ___AM.txt
    def write_to_am_txt(self, output_folder, metrics_file_prefix):
        data_list_to_write = [0] * 24
        data_list_to_write[0] = "PAIR"
        data_list_to_write[1] = self.TOTAL_READS
        data_list_to_write[5] = self.PF_READS_ALIGNED
        data_list_to_write[16] = self.READS_ALIGNED_IN_PAIRS
        
        tab = "\t"
        newline = "\n"
        header = "CATEGORY{0}TOTAL_READS{0}PF_READS{0}PCT_PF_READS{0}PF_NOISE_READS{0}PF_READS_ALIGNED{0}PCT_PF_READS_ALIGNED{0}PF_ALIGNED_BASES{0}PF_HQ_ALIGNED_READS{0}PF_HQ_ALIGNED_BASES{0}PF_HQ_ALIGNED_Q20_BASES{0}PF_HQ_MEDIAN_MISMATCHES{0}PF_MISMATCH_RATE{0}PF_HQ_ERROR_RATE{0}PF_INDEL_RATE{0}MEAN_READ_LENGTH{0}READS_ALIGNED_IN_PAIRS{0}PCT_READS_ALIGNED_IN_PAIRS{0}PF_READS_IMPROPER_PAIRS{0}PCT_PF_READS_IMPROPER_PAIRS{0}BAD_CYCLES{0}STRAND_BALANCE{0}PCT_CHIMERAS{0}PCT_ADAPTER{0}SAMPLE{0}LIBRARY{0}READ_GROUP".format(tab)

        write_to_file ="{}/{}___DRAGEN3_10_8___AM.txt".format(output_folder, metrics_file_prefix)
        print("Writing: " + write_to_file)
        data_line = ""
        for i in data_list_to_write:
            data_line = "{}{}{}".format(data_line, str(i), tab)
        
        with open(write_to_file, 'w') as _file:
            _file.write("#{}{}".format(metrics_file_prefix, newline))
            for i in range(5):
                _file.write("#{}".format(newline))
            _file.write("{}{}".format(header, newline))
            _file.write(data_line)

    # creat ___MD.txt picard stats format file
    def write_to_md_txt(self, output_folder, metrics_file_prefix):
        data_list_to_write = [0] * 10
        data_list_to_write[0] = metrics_file_prefix
        data_list_to_write[2] = self.READ_PAIRS_EXAMINED
        data_list_to_write[4] = self.UNMAPPED_READS
        data_list_to_write[6] = self.READ_PAIR_DUPLICATES
        data_list_to_write[8] = self.PERCENT_DUPLICATION
        
        tab = "\t"
        newline = "\n"
        header = "LIBRARY{0}UNPAIRED_READS_EXAMINED{0}READ_PAIRS_EXAMINED{0}SECONDARY_OR_SUPPLEMENTARY_RDS{0}UNMAPPED_READS{0}UNPAIRED_READ_DUPLICATES{0}READ_PAIR_DUPLICATES{0}READ_PAIR_OPTICAL_DUPLICATES{0}PERCENT_DUPLICATION{0}ESTIMATED_LIBRARY_SIZE".format(tab)
        
        write_to_file = "{}/{}___DRAGEN3_10_8___MD.txt".format(output_folder, metrics_file_prefix)
        print("Writing: " + write_to_file)
        data_line = ""
        for i in data_list_to_write:
            data_line = "{}{}{}".format(data_line, str(i), tab)
        
        with open(write_to_file, 'w') as _file:
            _file.write("#{0}{1}{1}{1}".format(metrics_file_prefix, newline))
            _file.write("{}{}".format(header, newline))
            _file.write(data_line)
            
    def write_to_wgs_txt(self, output_folder, metrics_file_prefix):
        data_list_to_write = [0] * 28
        data_list_to_write[1] = self.MEAN_TARGET_COVERAGE
        data_list_to_write[3] = self.MEDIAN_COVERAGE
        data_list_to_write[13] = self.PCT_1X
        data_list_to_write[15] = self.PCT_10X
        data_list_to_write[16] = self.PCT_15X
        data_list_to_write[17] = self.PCT_20X
        data_list_to_write[21] = self.PCT_50X
        data_list_to_write[26] = self.PCT_100X
        
        tab = "\t"
        newline = "\n"
        header = "GENOME_TERRITORY{0}MEAN_COVERAGE{0}SD_COVERAGE{0}MEDIAN_COVERAGE{0}MAD_COVERAGE{0}PCT_EXC_ADAPTER{0}PCT_EXC_MAPQ{0}PCT_EXC_DUPE{0}PCT_EXC_UNPAIRED{0}PCT_EXC_BASEQ{0}PCT_EXC_OVERLAP{0}PCT_EXC_CAPPED{0}PCT_EXC_TOTAL{0}PCT_1X{0}PCT_5X{0}PCT_10X{0}PCT_15X{0}PCT_20X{0}PCT_25X{0}PCT_30X{0}PCT_40X{0}PCT_50X{0}PCT_60X{0}PCT_70X{0}PCT_80X{0}PCT_90X{0}PCT_100X{0}FOLD_80_BASE_PENALTY{0}FOLD_90_BASE_PENALTY{0}FOLD_95_BASE_PENALTY{0}HET_SNP_SENSITIVITY{0}HET_SNP_Q".format(tab)
        
        write_to_file = "{}/{}___DRAGEN3_10_8___WGS.txt".format(output_folder, metrics_file_prefix)
        print("Writing: " + write_to_file)
        data_line = ""
        for i in data_list_to_write:
            data_line = "{}{}{}".format(data_line, str(i), tab)
                
        with open(write_to_file, 'w') as _file:
            _file.write("#{}{}".format(metrics_file_prefix, newline))
            for i in range(4):
                _file.write("#{}".format(newline))
            _file.write("{0}{1}{0}".format(newline, header))
            _file.write(data_line)


def process_one_sample(dragen_metrics_directory, work_directory, sample_id, metrics_file_prefix, sample_type):
    dragen_stats = DragenStats()
    dragen_stats.read_info_from_csv(sample_id, dragen_metrics_directory)
    dragen_stats.write_to_am_txt(work_directory, metrics_file_prefix)
    dragen_stats.write_to_md_txt(work_directory, metrics_file_prefix)
    if sample_type == "WGS":
        dragen_stats.write_to_wgs_txt(work_directory, metrics_file_prefix)
    

if __name__ == "__main__":
    # Usage: python dragenstats_csv_to_txt.py [dragen_stats_dir] [output_file_dir] [] [TYPE]
    
    # if there are 4 arguments process the one sample named
    # example: python scripts/dragen_csv_to_picard.py /Users/mcmanamd/Downloads/DRAGENRNA/RNA /Users/mcmanamd/Downloads/test FAUCI_0050_A22C3WKLT3___P08822___XPRO_1057_T_RNA_IGO_08822_VD_2 RNA
    if len(sys.argv) == 5:
        dragen_metrics_directory = sys.argv[1]
        work_directory = sys.argv[2]
        metrics_file_prefix = sys.argv[3]
        sample_type = sys.argv[4]  # WGS or RNA

        sample_id = metrics_file_prefix.split("___")[2]

        process_one_sample(dragen_metrics_directory, work_directory, sample_id, metrics_file_prefix, sample_type)

    # if there are 2 arguments to main() then Process all DRAGEN WGS stats in the entire directory
    # example: python3 /Users/luc/Documents/GitHub/igo-demux/scripts/dragen_csv_to_picard.py /Users/luc/Documents/GitHub/igo-demux/test/ /Users/luc/Documents/GitHub/igo-demux/test/result_test/
    if len(sys.argv) == 3:
        dragen_stats_folder = sys.argv[1]
        output_folder_path = sys.argv[2]
        print("Processing all DRAGEN stats .csv files in folder: " + dragen_stats_folder)

        sample_list = get_sample_list(dragen_stats_folder)
        for i in sample_list:
            print("Processing " + i)
            dragen_stats = DragenStats()
            dragen_stats.read_info_from_csv(i, dragen_stats_folder)
            dragen_stats.write_to_am_txt(output_folder_path, i)
            dragen_stats.write_to_md_txt(output_folder_path, i)
            dragen_stats.write_to_wgs_txt(output_folder_path, i)
