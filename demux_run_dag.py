import os
import subprocess
from datetime import datetime, timedelta
from SampleSheet import SampleSheet
import scripts.organise_fastq_split_by_lane
import pandas

from airflow import DAG
from airflow.operators.python import PythonOperator


with DAG(
    dag_id="demux_run",
    schedule_interval=None,
    start_date=datetime(2021, 1, 1),
    catchup=False,
    tags=["demux_run"],
) as dag:
    """
    Runs the demux such as: 
    bsub -n48 -q dragen /opt/edico/bin/dragen --bcl-conversion-only true --bcl-sampleproject-subdirectories true --force
      --bcl-input-directory /igo/sequencers/johnsawyers/211108_JOHNSAWYERS_0312 
      --output-directory /igo/staging/FASTQ/JOHNSAWYERS_0312 
      --sample-sheet /home/igo/DividedSampleSheets/SampleSheet_211108_JOHNSAWYERS_0312.csv 
    """

    def demux(ds, **kwargs):
        sequencer_path = kwargs["params"]["sequencer_path"]
        samplesheet_path = kwargs["params"]["samplesheet"]

        samplesheet = os.path.basename(samplesheet_path)
        samplesheet_no_ext = os.path.splitext(samplesheet)[0]  # SampleSheet_210331_MICHELLE_0360_BH5KFYDRXY
        sequencer_and_run = samplesheet_no_ext[19:]            # remove 'SampleSheet_210331_'

        sample_sheet = SampleSheet(pandas.DataFrame(),pandas.DataFrame(),"").read_csv(samplesheet_path)
        
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
        if is_DLP:
            is_DRAGEN_demux = False
            # demux with --no-lane-splitting
            # TODO setup bcl2fastq command for DLP
            command = "bsub -J {job_name} -o {output_log_file} -n 36 -M 6 /opt/common/CentOS_6/bcl2fastq/bcl2fastq2-v2.20.0.422/bin/bcl2fastq --minimum-trimmed-read-length 0 --mask-short-adapter-reads 0 --ignore-missing-bcl --runfolder-dir /igo/sequencers/{sequencer}/{run_name} --sample-sheet {sample_sheet} --output-dir /igo/work/FASTQ/{fastq_folder} --ignore-missing-filter --ignore-missing-positions --ignore-missing-control --barcode-mismatches 0 --no-lane-splitting --loading-threads 12 --processing-threads 24"
        elif is_10X:
            is_DRAGEN_demux = False
            # TODO for 10X build correct mkfastq command, special 10X barcodes can't go to dragen
            print("Building mkfastq command")
        else:
            # -K - wait for the job to complete
            bsub_command = "bsub -K -n48 -q dragen -e /igo/work/igo/igo-demux/logs/demux.log -o /igo/work/igo/igo-demux/logs/demux.log "
            command = bsub_command + "/opt/edico/bin/dragen --bcl-conversion-only true --force --bcl-sampleproject-subdirectories true --bcl-input-directory \'{}\' --output-directory \'{}\' --sample-sheet \'{}\'".format(
            sequencer_path, output_directory, samplesheet_path)
        print("Running demux command: " + command)
        subprocess.run(command, shell=True, check=True)
        # if the demux was successful:
        if is_DRAGEN_demux and not is_DLP:
            print("Adding sample sub-folders to the DRAGEN demux.")
            scripts.organise_fastq_split_by_lane.create_fastq_folders(output_directory)
        # TODO launch stats and/or pipeline for all projects on the run which needs stats/pipeline
        return command

    demux_run = PythonOperator(
        task_id='start_the_demux',
        python_callable=demux,
    )

    demux_run
""" 
Read the input arguments such as:

    'params': {'samplesheet': '/igo/work/igo/SampleSheetCopies/SampleSheet_211206_JOHNSAWYERS_0317_000000000-K3LFK.csv',
             'sequencer_path': '/igo/sequencers/johnsawyers/211206_JOHNSAWYERS_0317_000000000-K3LFK'},
"""