# launch cell ranger pipeline (GE, VDJ, ATAC....) for 10X samples by recipe
import pandas as pd
import re
import sys
import os
import json
import subprocess
import os.path
import shutil
import scripts.get_sequencing_read_data
import scripts.cellranger_spatial
import scripts.cellranger_config as CONFIG

"""
input: sample_sheet object(for sample list and essential info), sequencer_and_run(for stats folder and fastq file location)
output: running cmd for cellranger by sample

steps:
1. check whether there is previous fastq existing under /igo/staging/FASTQ (find_fastq_file)
2. get tag by recipe, if recipe not in the list above, skip for now (get_tag)
3. generate corresponding commands based on tag and genome (generate_cellranger_cmd)
4. run stats and create josn file for each project and call endpoint to push to qc website (launch_cellranger)
"""
# return a dictionary sample_ID -> list of fastq file path by given sample name list
def find_fastq_file(sample_ID_list, archive = False):
    # get whole list of all fastq files that available with project folder as tag
    if archive:
        path_prefix = "/igo/delivery/FASTQ/"
    else:
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

def get_tag(recipe):
    tag = "Skip"
    if recipe in CONFIG.COUNT_FLAVORS:
        tag = "count"
    if recipe in CONFIG.VDJ_T_FLAVORS:
        tag = "vdj_t"
    if recipe in CONFIG.VDJ_B_FLAVORS:
        tag = "vdj_b"
    if recipe in CONFIG.ARC_FLAVORS:
        tag = "arc"
    if recipe in CONFIG.SPATIAL_FLAVORS:
        tag = "spaceranger"
    return tag

def generate_cellranger_cmd(sample_ID, tag, genome, fastq_file_path, sequencer_and_run):
    vdj = False
    if tag == "vdj_t" or tag == "vdj_b":
        vdj = True
        if tag == "vdj_t":
            additional_option = " --chain TR"
        if tag == "vdj_b":
            additional_option = " --chain IG"
        tag = "vdj"

    tool = CONFIG.config_dict[tag]["tool"]
    transcriptome = CONFIG.config_dict[tag]["genome"][genome]
    project_ID = "Project_" + "_".join(sample_ID.split("_")[sample_ID.split("_").index("IGO") + 1:-1])
    cellranger_cmd = "{}--id=Sample_{}__{}".format(tool, sample_ID, tag) + transcriptome + "--fastqs=" + ",".join(fastq_file_path) + CONFIG.OPTIONS
    if vdj:
        cellranger_cmd = cellranger_cmd.replace(" --create-bam=true", "")
        cellranger_cmd = cellranger_cmd + additional_option
    job_name = "{}_{}_{}_{}_cellranger".format(sequencer_and_run, project_ID, sample_ID, tag)
    bsub_cmd = "bsub -J {} -o {}.out{}".format(job_name, job_name, cellranger_cmd) 
    return bsub_cmd
        
def create_json(send_json, sequencer_and_run, project, work_area):  
    job_id = sequencer_and_run + "_" + project            
    json_data_file = "cellranger_json___" + sequencer_and_run + "__" + project + ".json"
    with open(json_data_file, "w") as jfile:
        json.dump(send_json, jfile)
        
    bsub_json = "bsub -J create_json___{} -o create_json___{}.log -w \"done({}*)\" sh /igo/work/igo/igo-demux/scripts/send_json_data.sh {} {}".format(job_id, job_id, job_id, work_area, json_data_file)
    print(bsub_json)
    subprocess.run(bsub_json, shell = True)

def create_library_csv_file(ge_sample_path, atac_sample_path, sample_ID):
    """
    ge_sample_path and atac_sample_path will be list just in case if top up is performed
    """
    with open("Sample_{}.csv".format(sample_ID),"a") as file:
        file.write("fastqs,sample,library_type\n")
        for ge in ge_sample_path:
            file.write("{},{},Gene Expression\n".format(ge, sample_ID))
        for atac in atac_sample_path:
            file.write("{},{},Chromatin Accessibility\n".format(atac, sample_ID))

def get_sequencer_runID(fastq_path):
    runID = fastq_path.split("/")[4]
    runID = "_".join(runID.split("_")[0:3])
    sequencer = runID.split("_")[0].lower()
    return sequencer, runID

def multiome_valid(fastq_list):
    """
    check whether the list of fastq files contain both GE and ATAC
    return a list which first item will be YES or NO for validation result
    if yes, second item will be GE fastq list and third item will be ATAC fastq list
    """
    sequencer_prefix = "/igo/sequencers/"
    is_valid = "NO"
    ge_list = []
    atac_list = []
    if len(fastq_list) == 1:
        return [is_valid, ge_list, atac_list]
    else:
        for fastq in fastq_list:
            # find corresponding run folder and check if is atac run
            # append the fastq file to corresponding list
            sequencer, runID = get_sequencer_runID(fastq)
            # find sequenceing path
            if sequencer == "pepe" or sequencer == "amelie":
                sequencer_path_pre = sequencer_prefix + sequencer + "/output"
            else:
                sequencer_path_pre = sequencer_prefix + sequencer
            run_list = os.listdir(sequencer_path_pre)
            for run in run_list:
                if re.match(".*" + runID, run):
                    sequencer_path = sequencer_path_pre + "/" + run
                    atac = scripts.get_sequencing_read_data.get_sequencing_read_data(sequencer_path)[0]
                    if atac:
                        atac_list.append(fastq)
                    else:
                        ge_list.append(fastq)
                    break

    if len(ge_list) > 0 and len(atac_list) > 0:
        is_valid = "YES"
    
    return [is_valid, ge_list, atac_list]

# lanuch cellranger per project
def lanuch_by_project(sequencer_and_run, project, sample_id_list, sample_genome_dict, sample_recipe_dict, archive = False):
    sample_fastqfile_dict = find_fastq_file(sample_id_list, archive)
    send_json = {}
    send_json["samples"] = []
    # CREATE RUN FOLDER AND PROJECT FOLDER IF NOT ALREADY THERE    
    work_area = CONFIG.STATS_AREA + sequencer_and_run + "/" + project + "/" 
    if not os.path.exists(work_area):
        os.makedirs(work_area, CONFIG.ACCESS)

    # GO TO project ID LOCATION to start cellranger command
    os.chdir(work_area)

    # call cellranger for each sample and append info to json dict
    for sample in sample_id_list:
        if sample_genome_dict[sample] != "Human" and sample_genome_dict[sample] != "Mouse":
            sample_genome_dict[sample] = "Mouse"
        tag = get_tag(sample_recipe_dict[sample])
        # if recipe within the tool being set up, lanuch cellranger
        if tag == "arc":
            validation = multiome_valid(sample_fastqfile_dict[sample])
            if validation[0] == "YES":
                create_library_csv_file(validation[1], validation[2], sample)
                tool = CONFIG.config_dict[tag]["tool"]
                transcriptome = CONFIG.config_dict[tag]["genome"][sample_genome_dict[sample]]
                cmd = "{}--id=Sample_{}{}".format(tool, sample, transcriptome) + "--libraries={}Sample_{}.csv".format(work_area, sample) + CONFIG.ARC_OPTIONS
                bsub_cmd = "bsub -J {}_{}_{}_ARC -o {}_ARC.out{}".format(sequencer_and_run, project, sample, sample, cmd)
                print(bsub_cmd)
                subprocess.run(bsub_cmd, shell=True)
            else:
                print("Multiome sample set not complete yet")
        elif tag == "spaceranger":
            sample_info = scripts.cellranger_spatial.Spatial_sample(sample, project)
            if sample_info.tiff_image == "EMPTY":
                print("check tif image for sample {}".format(sample))
            else:
                tool = CONFIG.config_dict[tag]["tool"]
                transcriptome = CONFIG.config_dict[tag]["genome"][sample_genome_dict[sample]]
                cmd = "{}--id=Sample_{}{}".format(tool, sample, transcriptome) + "--fastqs=" + ",".join(sample_fastqfile_dict[sample]) + " --image={} --slide={} --area={}".format(sample_info.tiff_image, sample_info.chip_id, sample_info.chip_position)
                
                if sample_info.cytAssist:
                    cmd = "{}--id=Sample_{}{}".format(tool, sample, transcriptome) + "--fastqs=" + ",".join(sample_fastqfile_dict[sample]) + " --cytaimage={} --slide={} --area={}".format(sample_info.tiff_image, sample_info.chip_id, sample_info.chip_position)
                    if sample_genome_dict[sample] == "Human":
                        probe = CONFIG.config_dict[tag]["probe"]["Human_CytAssist"]
                    elif sample_genome_dict[sample] == "Mouse":
                        if sample_info.chip_id.startswith("H1"):
                            probe = CONFIG.config_dict[tag]["probe"]["Mouse_HD"]
                        else:
                            probe = CONFIG.config_dict[tag]["probe"]["Mouse"]
                    cmd = cmd + " --probe-set={}".format(probe)
                        
                elif sample_info.preservation == "FFPE":
                    probe = CONFIG.config_dict[tag]["probe"][sample_genome_dict[sample]]
                    cmd = cmd + " --probe-set={}".format(probe)
                
                # Eventhough HE image is required internal, the pipeline doesn't need it. Add it if exists
                if sample_info.HE_tiff_image != "EMPTY":
                    cmd = cmd + " --image={}".format(sample_info.HE_tiff_image)
                    # copy microsope image here in sub folder for delivery 
                    HE_folder_loc = work_area + "Microscope/"
                    if not os.path.exists(HE_folder_loc):
                        os.makedirs(HE_folder_loc)
                    shutil.copy(sample_info.HE_tiff_image , HE_folder_loc)
                
                # if there is manual alignment json file availabe, add that to the cmd
                if sample_info.json != "EMPTY":
                    cmd = cmd + " --loupe-alignment={}".format(sample_info.json)

                bsub_cmd = "bsub -J {}_{}_{}_SPATIAL -o {}_SPATIAL.out{}{}".format(sequencer_and_run, project, sample, sample, cmd, CONFIG.OPTIONS)
                print(bsub_cmd)
                subprocess.run(bsub_cmd, shell=True)

        elif tag != "Skip":
            cmd = generate_cellranger_cmd(sample, tag, sample_genome_dict[sample], sample_fastqfile_dict[sample], sequencer_and_run)
            print(cmd)
            subprocess.run(cmd, shell=True)
            send_json["samples"].append({"sample":"Sample_" + sample, "type":tag, "project":project, "run":sequencer_and_run})
    
    if send_json["samples"]:
        create_json(send_json, sequencer_and_run, project, work_area)

# Main function: launch cellranger cmd by given samplesheet object and sequencer_and_run
def launch_cellranger_by_sample_sheet(sample_sheet, sequencer_and_run):
    # get parameters from sample_sheet
    # dictionary of Sample_ID->Project
    sample_project_dict = pd.Series(sample_sheet.df_ss_data["Sample_Project"].values,index=sample_sheet.df_ss_data["Sample_ID"]).to_dict()
    # dictionary of project->sample_ID
    project_sample_dict = {}
    for sample_ID, project_ID in sample_project_dict.items():
        if project_ID in project_sample_dict.keys():
            project_sample_dict[project_ID].append(sample_ID)
        else:
            project_sample_dict[project_ID] = [sample_ID]
    # dictionary of sample_ID->recipe
    sample_recipe_dict = pd.Series(sample_sheet.df_ss_data["Sample_Well"].values,index=sample_sheet.df_ss_data["Sample_ID"]).to_dict()
    # dictionary of sample_ID->genome
    sample_genome_dict = pd.Series(sample_sheet.df_ss_data["Sample_Plate"].values,index=sample_sheet.df_ss_data["Sample_ID"]).to_dict()
    # launch cellranger cmd for each project
    for project in project_sample_dict.keys():
        # SCRI or SAIL samples don't need to run cellranger
        if (not any(prj in project for prj in CONFIG.DO_NOT_PROCESS)):
            sample_list = project_sample_dict[project]
            lanuch_by_project(sequencer_and_run, project, sample_list, sample_genome_dict, sample_recipe_dict)

def launch_cellranger_by_project_location(project_directory, recipe, species):
    # get sample_ID list
    sample_list_ori = os.listdir(project_directory)
    sample_list = []
    for sample in sample_list_ori:
        # remove Sample_ prefix
        sample_list.append(sample[7:])
    # get project and run info from project_directory
    project = project_directory.split("/")[5]
    sequencer_and_run = project_directory.split("/")[4]
    sample_genome_dict = {}
    sample_recipe_dict = {}
    for sample in sample_list:
        sample_genome_dict[sample] = species
        sample_recipe_dict[sample] = recipe
    # add checker for if fastq file should come from archived location
    archive = False
    if "delivery" in project_directory:
        archive = True
    lanuch_by_project(sequencer_and_run, project, sample_list, sample_genome_dict, sample_recipe_dict, archive)


if __name__ == '__main__':
    # launch cellranger commands by project
    # Usage: python cellranger.py [project_directory] [recipe] [species]
    # example: python cellranger.py /igo/staging/FASTQ/RUTH_0141_AH27NGDSX5/Project_13586_B 10X_Genomics_GeneExpression-3 Human
    project_directory = sys.argv[1]
    recipe = sys.argv[2]
    species = sys.argv[3]
    launch_cellranger_by_project_location(project_directory, recipe, species)
