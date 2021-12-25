import os
import subprocess
from datetime import datetime, timedelta

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
        samplesheet_path = kwargs["params"]["samplesheet"]
        samplesheet = os.path.basename(samplesheet_path)
        samplesheet_no_ext = os.path.splitext(samplesheet)[0]  # SampleSheet_210331_MICHELLE_0360_BH5KFYDRXY
        sequencer_and_run = samplesheet_no_ext[19:]            # remove 'SampleSheet_210331_'
        sequencer_path = kwargs["params"]["sequencer_path"]

        # TODO for some 10X build correct mkfastq command, special 10X barcodes can't go to dragen

        output_directory = "/igo/staging/FASTQ/" + sequencer_and_run + "_DGN"
        # -K - wait for the job to complete
        bsub_command = "bsub -K -n48 -q dragen -e /igo/work/igo/igo-demux/logs/demux.log -o /igo/work/igo/igo-demux/logs/demux.log "
        command = bsub_command + "/opt/edico/bin/dragen --bcl-conversion-only true --force --bcl-sampleproject-subdirectories true --bcl-input-directory \'{}\' --output-directory \'{}\' --sample-sheet \'{}\'".format(
            sequencer_path, output_directory, samplesheet_path)
        print("Running demux: " + command)
        subprocess.run(command, shell=True, check=True)
        # if the demux was successful:
        # TODO for non DLP call organise_fastq_split_by_lane.py to create Sample sub-dirs like bcl2fastq
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