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
            output_directory = "/igo/work/FASTQ/" + sequencer_and_run + "_DLP"
        if is_10X:
            output_directory = "/igo/work/FASTQ/" + sequencer_and_run + "_10X"
        if "REFERENCE" in samplesheet_path:
            output_directory = "/igo/work/FASTQ/" + sequencer_and_run + "_REFERENCE"
        else:
            output_directory = "/igo/work/FASTQ/" + sequencer_and_run
        
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

        return command


    def stats(ds, **kwargs):
        sequencer_path = kwargs["params"]["sequencer_path"]
        samplesheet_path = kwargs["params"]["samplesheet"]
        samplesheet = os.path.basename(samplesheet_path)
        samplesheet_no_ext = os.path.splitext(samplesheet)[0]  # SampleSheet_210331_MICHELLE_0360_BH5KFYDRXY
        sequencer_and_run = samplesheet_no_ext[19:]            # remove 'SampleSheet_210331_'

        sample_sheet = SampleSheet(samplesheet_path)
        if "REFERENCE" in samplesheet_path:
            return "No stats for reference "  + samplesheet_path

        if "DLP" in sample_sheet.recipe_set:
           return "No DLP stats"
        
        launch_stats(sample_sheet, sequencer_and_run)

        return "Completed"

    demux_run = PythonOperator(
        task_id='start_the_demux',
        python_callable=demux,
    )

    launch_stats = PythonOperator(
        task_id='launch_stats',
        python_callable=stats,
    )

    demux_run >> launch_stats

    """
    Process dictionary of sample sheet project,recipe and launch stats for each project on the run.
    """
    def launch_stats(sample_sheet, sequencer_and_run):
        working_dir = "/igo/staging/stats/" + sequencer_and_run
        if not os.path.exists(working_dir):
            os.mkdir(working_dir)
        os.chdir(working_dir)

        # TODO email demux complete starting stats for non "REFERENCE" demuxes

        # sh /home/igo/Scripts/Automate-Stats/Launch_Stats_GRCh37.sh 191029_JOHNSAWYERS_0218_000000000-G4BYY 2>&1 >> /home/igo/log/statsTest.log
        cmd = "sh /home/igo/Scripts/Automate-Stats/Launch_Stats_Airflow.sh {} 2>&1 >> /igo/staging/stats.log".format(sequencer_and_run)
        print(cmd)
        subprocess.run(cmd, shell=True)
        #TODO Launch DetectStatsCompletion.sh
        # sh $SCRIPTPATH/../Automate-Stats/DetectStatsCompletion.sh $RUNNAME &