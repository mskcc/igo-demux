import sequencer
import os
from datetime import datetime

from airflow.operators.bash import BashOperator
from airflow.decorators import dag, task
from airflow.models import Variable

@dag(dag_id='find_completed_runs', schedule_interval=None, start_date=datetime(2021, 1, 1), catchup=False, tags=['find_completed_runs'])
def find_completed_runs():
    """
    #### Find Completed Runs
    Find recently completed runs by looking for the last file written by the sequencers then copy the sample sheet for the completed run to the list of completed runs variable
    so a later task can start the demux for the sample sheet.
    """
    @task()
    def find_runs(config_file):
        """
        #### Find Runs task
        Finds the list of runs which are recently completed.
        """
        sequencers = sequencer.read_config(config_file)
        time_to_search = Variable.get("completed_run_search_interval_mins", defalt_var=30)
        completed_runs_path = sequencer.find_completed_runs(sequencers, time_to_search)

        for run_path in completed_runs_path:
            copy_samplesheet(run_path)

    def copy_samplesheet(completed_run_path):
        print("Preparing sample sheet(s) for completed run:" + completed_run_path)
        completed_run = str(os.path.split(completed_run_path).tail)
        samplesheet = "SampleSheet_" + completed_run + ".csv"
        orig_samplesheet = "/pskis34/LIMS/LIMS_SampleSheets/" + samplesheet
        dest_samplesheet = "/igo/work/igo/SampleSheetCopies/" + samplesheet
        
        cp_command = "cp {} {}".format(orig_samplesheet, dest_samplesheet)
        print("Copying sample sheet:" + cp_command)
        run_cp_task = BashOperator(
            task_id='copy_samplesheet',
            bash_command=cp_command,
        )

        # TODO Split sample sheet for DLP, PED-PEG & 10X

        # append new sample sheet to sample sheets variable list
        samplesheets_list = Variable.get("samplesheets")
        samplesheets_list.append(dest_samplesheet)
        Variable.set("samplesheets", samplesheets_list)
        print("Current sample sheet list ready for demux: " + samplesheets_list)
        
    find_runs("dag_runs.conf")

find_completed_runs = find_completed_runs()