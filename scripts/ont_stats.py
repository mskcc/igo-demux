import pandas as pd
import statistics
import sys
import glob
import os
from collections import OrderedDict

# TODO check for multiple run
def get_read_length_and_summary(file_path):
    summary_metrix = pd.read_csv(file_path, delimiter = "\t")
    read_length = summary_metrix[summary_metrix["passes_filtering"]]["sequence_length_template"].tolist()
    read_length.sort(reverse = True)
    median = statistics.median(read_length)
    N50_value = sum(read_length) / 2
    total = 0
    for item in read_length:
        total += item
        if total >= N50_value:
            N50 = item
            break
    
    return(len(read_length), N50_value * 2 / 1000000000, N50, median)

def write_to_csv(sample_dict):
    file_name = "summary.csv"
    with open(file_name,'w') as file:
        file.write("sample_id, Reads, Bases, N50, Meidan Read Length\n")
        for key, value in sample_dict.items():
            file.write("{}, {}, {}, {}, {}\n".format(key, value[0], value[1], value[2], value[3]))

if __name__ == '__main__':
    # Usage: python ont_stats.py [project_directory]
    # example: python ont_stats.py /igo/staging/promethion/Project_14607
    
    project_directory = sys.argv[1]
    os.chdir(project_directory)
    sample_list = next(os.walk("."))[1]
    sample_dict = OrderedDict()
    sample_list.sort()
    for sample in sample_list:
        destination = project_directory + "/" + sample
        file = glob.glob(destination + "/*/sequencing_summary_*")
        if len(file) != 0:
            sample_dict[sample] = get_read_length_and_summary(file[0])
    
    write_to_csv(sample_dict)
