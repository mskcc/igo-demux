
import demux
import os
from datetime import datetime

from airflow.operators.bash import BashOperator
from airflow.decorators import dag, task

@dag(dag_id='demultiplex', schedule_interval=None, start_date=datetime(2021, 1, 1), catchup=False, tags=['demultiplex'])
def demultiplex_runs():
    """
    ### TaskFlow API Documentation
    [here](https://airflow.apache.org/docs/apache-airflow/stable/tutorial_taskflow_api.html)
    """
    @task()
    def find_runs(config_file):
        """
        #### Find Runs task
        Finds the list of runs which are recently completed.
        """
        sequencers = demux.read_config(config_file)
        runs_to_demux_list = demux.find_completed_runs(sequencers, 30)
        return runs_to_demux_list

    @task()
    def demux_runs(runs_to_demux_list):
        """
        #### Demux Runs task
        Run the DRAGEN command to generate the fastq.gz files
        """
        for run_value in runs_to_demux_list:
            demux_run(run_value)
        return

    def demux_run(run_path):
        print("Demultiplexing " + run_path)
        demux_dir = str(os.path.split(run_path).tail)
        samplesheet = "SampleSheet_" + demux_dir + ".csv"
        #TODO just leave these paths right here for now
        orig_samplesheet = "/pskis34/LIMS/LIMS_SampleSheets/" + samplesheet
        dest_samplesheet = "/home/igo/SampleSheetCopies/" + samplesheet
        
        cp_command = "cp {} {}".format(orig_samplesheet,dest_samplesheet)
        print("Copying sample sheet " + cp_command)
        run_this = BashOperator(
            task_id='run_after_loop',
            bash_command='echo 1',
        )

        # bash_copy_samplesheet_task = BashOperator(...
        # Split sample sheet for DLP, PED-PEG & 10X
        # for DLP run make create-metadata-yaml ss=${SAMPLESHEET} prj=${prj} project_path=${project_path}
        # execute DRAGEN demux command(s)

    runs_to_demux = find_runs("dag_runs.conf")
    demux_runs(runs_to_demux)

demultiplex_runs = demultiplex_runs()