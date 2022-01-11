import os
from re import sub
import subprocess
from datetime import datetime, timedelta
from SampleSheet import SampleSheet
import scripts.organise_fastq_split_by_lane
import pandas

from airflow import DAG
from airflow.operators.python import PythonOperator


"""
    Runs the demux such as: 
    bsub -n48 -q dragen /opt/edico/bin/dragen --bcl-conversion-only true --bcl-sampleproject-subdirectories true --force
      --bcl-input-directory /igo/sequencers/johnsawyers/211108_JOHNSAWYERS_0312 
      --output-directory /igo/staging/FASTQ/JOHNSAWYERS_0312 
      --sample-sheet /home/igo/DividedSampleSheets/SampleSheet_211108_JOHNSAWYERS_0312.csv 
"""
with DAG(
    dag_id="demux_run",
    schedule_interval=None,
    start_date=datetime(2022, 1, 1),
    catchup=False,
    tags=["demux_run"],
) as dag:

    """ 
    Read the input arguments such as:

    'params': {'samplesheet': '/igo/work/igo/SampleSheetCopies/SampleSheet_211206_JOHNSAWYERS_0317_000000000-K3LFK.csv',
             'sequencer_path': '/igo/sequencers/johnsawyers/211206_JOHNSAWYERS_0317_000000000-K3LFK'},
    """
    def demux(ds, **kwargs):
        sequencer_path = kwargs["params"]["sequencer_path"]
        samplesheet_path = kwargs["params"]["samplesheet"]

        samplesheet = os.path.basename(samplesheet_path)
        samplesheet_no_ext = os.path.splitext(samplesheet)[0]  # SampleSheet_210331_MICHELLE_0360_BH5KFYDRXY
        sequencer_and_run = samplesheet_no_ext[19:]            # remove 'SampleSheet_210331_'

        sample_sheet = SampleSheet(samplesheet_path)
        
        is_DLP = False
        if "DLP" in sample_sheet.recipe_set:
            is_DLP = True
        is_10X = False
        if len(sample_sheet.barcode_list_10X) > 0:
            is_10X = True

        if is_DLP:
            output_directory = "/igo/staging/FASTQ/" + sequencer_and_run + "_DLPDGN"
        if is_10X:
            output_directory = "/igo/staging/FASTQ/" + sequencer_and_run + "_10XDGN"
        else:
            output_directory = "/igo/staging/FASTQ/" + sequencer_and_run + "_DGN"
        
        is_DRAGEN_demux = True
        
        if is_10X:
            is_DRAGEN_demux = False
            # TODO for 10X build correct mkfastq command, special 10X barcodes can't go to dragen
            print("Building mkfastq command")
        else:
            # DLP can demux with the default command as long as the [Settings] have 'NoLaneSplitting,true'
            # -K - wait for the job to complete
            bsub_command = "bsub -K -n48 -q dragen -e /igo/work/igo/igo-demux/logs/demux.log -o /igo/work/igo/igo-demux/logs/demux.log "
            command = bsub_command + "/opt/edico/bin/dragen --bcl-conversion-only true --bcl-only-matched-reads true --force --bcl-sampleproject-subdirectories true --bcl-input-directory \'{}\' --output-directory \'{}\' --sample-sheet \'{}\'".format(
            sequencer_path, output_directory, samplesheet_path)
        print("Running demux command: " + command)
        subprocess.run(command, shell=True, check=True)
        # if the demux was successful:
        if is_DRAGEN_demux and not is_DLP:
            print("Adding sample sub-folders to the DRAGEN demux.")
            scripts.organise_fastq_split_by_lane.create_fastq_folders(output_directory)
        
        # for DLP projects create the .yaml file
        if is_DLP:
            # example: make create-metadata-yaml ss=/igo/home/igo/DividedSampleSheets/SampleSheet_211022_DIANA_0415_AHMJWMDSX2_DLP.csv prj=Project_09443_CI project_path=/igo/staging/FASTQ/DIANA_0415_AHMJWMDSX2_DLP/Project_09443_CI
            for project in sample_sheet.project_set:
                fastq_project_dir = output_directory + "/" + project
                make_command = "make create-metadata-yaml ss={} prj={} project_path={}".format(samplesheet_path, project, fastq_project_dir)
                print("Calling DLP make command: {}".format(make_command))
                subprocess.check_output(make_command, cwd="/home/igo/shared-single-cell", shell=True)

        # TODO email demux complete starting stats for non "REFERENCE" demuxes
        if "REFERENCE" in samplesheet_path:
            return command

        launch_stats(sample_sheet, output_directory, sequencer_and_run)

        return command

    demux_run = PythonOperator(
        task_id='start_the_demux',
        python_callable=demux,
    )

    demux_run

    """
    Process dictionary of sample sheet project,recipe and launch stats for each project on the run.
    """
    def launch_stats(sample_sheet, output_directory, sequencer_and_run):
        nf_working_dir = "/igo/staging/working/" + sequencer_and_run
        print("Creating nextflow working directory - {}".format(nf_working_dir))
        if not os.path.exists(nf_working_dir):
            os.mkdir(nf_working_dir)
        os.chdir(nf_working_dir)

        for project, recipe in sample_sheet.project_dict.items():
            cmd_basic = "nohup /home/igo/bin/nextflow /home/igo/nf-fastq-plus/samplesheet_stats_main.nf"
            cmd = "{} --ss {} --dir {}  --filter {}".format(cmd_basic, sample_sheet.path, output_directory, project.replace("Project_", ""))
            print(project, recipe, cmd)
            subprocess.run(cmd, shell=True)