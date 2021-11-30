import os

from datetime import datetime, timedelta

from airflow import DAG
from airflow.decorators import task
from airflow.operators.bash import BashOperator
from airflow.models import Variable

with DAG(
    dag_id='demux_run',
    tags=['demux_run'],
) as dag:
    """
    Runs the demux such as: 
    bsub -n48 -q dragen $JOB /opt/edico/bin/dragen --bcl-conversion-only true --bcl-sampleproject-subdirectories --force
      --bcl-input-directory /igo/sequencers/johnsawyers/211108_JOHNSAWYERS_0312 
      --output-directory /igo/staging/FASTQ/JOHNSAWYERS_0312 
      --sample-sheet /home/igo/DividedSampleSheets/SampleSheet_211108_JOHNSAWYERS_0312.csv 
    """
    # get the sample sheet path and bcl directory
    samplesheets_list = Variable.get("ready_to_demux")
    sample_sheet_path = samplesheets_list.pop(0)
    sequencer_run_dir = samplesheets_list.pop(0)
    Variable.set("ready_to_demux", samplesheets_list)
    
    sample_sheet = os.path.basename(sample_sheet_path)
    sample_sheet_no_ext = os.path.splitext(sample_sheet)[0]

    output_directory = "/igo/staging/FASTQ/" + str.replace(sample_sheet_no_ext("SampleSheet_","")) + "_DGN"
    job_name = "demux_" + sample_sheet
    
    # build the demux command
    command = "bsub -n48 -q dragen /opt/edico/bin/dragen --bcl-conversion-only true --force --bcl-sampleproject-subdirectories --bcl-input-directory {} --output-directory {} --sample-sheet {}"
    print("Running demux: " + command)
    
    demux_command = BashOperator(
        task_id='demux_command',
        bash_command=command,
    )
    
    demux_command

if __name__ == "__main__":
    dag.cli()