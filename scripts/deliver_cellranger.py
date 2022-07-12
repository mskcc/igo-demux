import os
import sys
import glob

# given project ID, look through cellranger folder and return a list of path of folders need to copy

CELLRANGER_DIR = '/igo/stats/CELLRANGER/'
# structure '/igo/stats/CELLRANGER/RUNNAME/PROJECTID/SAMPLEFOLDER

# find all the cellranger result given project ID, return a list of address
def find_cellranger(project):
    cellranger_list = []
    project_ID = "Project_" + project

    # get whole list of all cellranger folders that available with project folder as tag
    cellranger_list_full = os.listdir(CELLRANGER_DIR)
    # dictionary of run_ID->project_list
    run_project_dict = {}
    for run_ID in cellranger_list_full:
        current_path = CELLRANGER_DIR + run_ID + "/"
        if os.path.isdir(current_path):
            run_project_dict[run_ID] = []
            file_list = os.listdir(current_path)
            for item in file_list:
                if "Project_" in item:
                    run_project_dict[run_ID].append(item)

    # dictionary of sample->cellranger file list
    cellranger_file_list_dict = {}
    for run_ID, project_list in run_project_dict.items():
        if project_ID in project_list:
            cellranger_file_path_prefix = CELLRANGER_DIR + run_ID + "/" + project_ID + "/"
            sample_folder_list = os.listdir(cellranger_file_path_prefix)
            for sample_folder in sample_folder_list:
                sample_folder_fullpath = cellranger_file_path_prefix + sample_folder
                if os.path.isdir(sample_folder_fullpath):
                    print("cellranger output folder found: {}".format(sample_folder_fullpath))
                    sample = sample_folder.split("/")[-1]
                    if sample in cellranger_file_list_dict.keys():
                        cellranger_file_list_dict[sample].append(sample_folder_fullpath)
                    else:
                        cellranger_file_list_dict[sample] = [sample_folder_fullpath]

    # for multi cellranger files for same sample, keep only the newest one
    for sample, cellrangers in cellranger_file_list_dict.items():
        if len(cellrangers) == 1:
            cellranger_list.append(cellrangers[0])
        else:
            latest = cellrangers[0]
            for cellranger in cellrangers:
                if os.path.getmtime(cellranger) > os.path.getmtime(latest):
                    latest = cellranger
            cellranger_list.append(latest)

    return cellranger_list

