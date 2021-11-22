
import demux
from datetime import datetime

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
            print("Demultiplexing " + run_value)
            # TODO copy sample sheet from /pskis to /home/igo/SampleSheetCopies
            # bash_copy_samplesheet_task = BashOperator(...
            # Split sample sheet for DLP, PED-PEG & 10X
            # for DLP run make create-metadata-yaml ss=${SAMPLESHEET} prj=${prj} project_path=${project_path}
            # execute DRAGEN demux command(s)
        return

    runs_to_demux = find_runs("dag_runs.conf")
    demux_runs(runs_to_demux)

demultiplex_runs = demultiplex_runs()