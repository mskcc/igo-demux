import os
import shutil

def isBarcoded(experiment_directory):
    for root, dirs, files in os.walk(experiment_directory):
        if 'barcode01' in dirs:
            return True
    return False

def renameBarcodedSamples(experiment_directory):
    #search lims db for pool name if pool name is an IGO ID (data -> experiemnt -> pool name)
    #fetch parent samples of pool, find their barcodes
    #match parent sample IGO IDs to barcoded directories
    #rename

    # for if sample name is [pool ID] + ["_barcode_..."]
    pass


def reorganizeExperiment(experiment_path):
        sample_list = os.listdir(experiment_path)
        for sample_dir in sample_list:
            sample_path = os.path.join(experiment_path, sample_dir)
            if not os.path.isdir(sample_path):
                continue
            for alias_dir in os.listdir(sample_path):
                alias_path = os.path.join(sample_path, alias_dir)
                if not os.path.isdir(alias_path):
                    continue
                for fastq_pod5_dir in os.listdir(alias_path):
                    fastq_pod5_path = os.path.join(alias_path,fastq_pod5_dir)
                    if not os.path.isdir(fastq_pod5_path):
                        continue
                    for barcode_dir in os.listdir(fastq_pod5_path):
                        barcode_path = os.path.join(fastq_pod5_path,barcode_dir)
                        if not os.path.isdir(barcode_path):
                            continue
                        files_in_barcode = os.listdir(barcode_path)
                        if len(files_in_barcode) > 1:
                            #make new directories
                            new_sample_dir = sample_dir + "_" + barcode_dir
                            new_sample_path = os.path.join(experiment_path, new_sample_dir)
                            new_alias_path = os.path.join(new_sample_path, alias_dir)
                            new_fastq_pod5_path = os.path.join(new_alias_path, fastq_pod5_dir)
                            os.makedirs(new_fastq_pod5_path, exist_ok=True)
                            #copy
                            for file in files_in_barcode:
                                file_path = os.path.join(barcode_path, file)
                                print("File Path: ",file_path)
                                print("Dest Path: ", new_fastq_pod5_path)
                                shutil.copy2(file_path, new_fastq_pod5_path)

def renameNonBarcodeSamples():
    # for if sample name is IGO ID
    # make sample name full deliverable name
    pass


def OrganizeData(data_path):
    for experiment in os.listdir(data_path):
        experiment_path = os.path.join(data_path, experiment)
        if os.path.isdir(experiment_path):
            if isBarcoded(experiment_path):
                reorganizeExperiment(experiment_path)
                renameBarcodedSamples(experiment_path)
    renameNonBarcodeSamples()


if __name__ == "__main__":
    data_path = "/Users/desmondlambe/data"
    OrganizeData(data_path)


