import sequencer
import os
from datetime import datetime

from airflow.operators.trigger_dagrun import TriggerDagRunOperator
from airflow.operators.bash import BashOperator
from airflow.decorators import dag, task
from airflow.models import Variable

# defines the list of all sequencers and for each sequencer 1) name 2) location it writes runs to and 3) the last file the sequencer writes when a run is completed to signal demux can begin
sequencers = {"sequencers":[{"name":"ayyan","path":"/igo/sequencers/ayyan","last_file":"RTAComplete.txt"},{"name":"diana","path":"/igo/sequencers/diana","last_file":"CopyComplete.txt"},{"name":"michelle","path":"/igo/sequencers/michelle","last_file":"CopyComplete.txt"},{"name":"ruth","path":"/igo/sequencers/ruth","last_file":"CopyComplete.txt"},{"name":"johnsawyers","path":"/igo/sequencers/johnsawyers","last_file":"RTAComplete.txt"},{"name":"pepe","path":"/igo/sequencers/pepe/output","last_file":"CopyComplete.txt"},{"name":"scott","path":"/igo/sequencers/scott","last_file":"RunCompletionStatus.xml"}]}

@dag(dag_id='find_completed_runs', schedule_interval=None, start_date=datetime(2021, 1, 1), catchup=False, tags=['find_completed_runs'])
def find_completed_runs():
    """
    #### Find Completed Runs
    Find recently completed runs by looking for the last file written by the sequencers then copy the sample sheet for the completed run to the list of completed runs variable
    so a later task can start the demux for the sample sheet.
    """
    @task()
    def find_completed_runs():
        """
        #### Find Runs task
        Finds the list of runs which are recently completed.
        """
        print("Processing sequencer list: {}".format(sequencers))
        
        time_to_search = Variable.get("completed_run_search_interval_mins", default_var=30)
        print("Searching for runs completed in the last {} minutes, variable completed_run_search_interval_mins".format(time_to_search))
        completed_runs_path = sequencer.find_completed_runs(sequencers, time_to_search)

        for run_path in completed_runs_path:
            copy_samplesheet(os.path.dirname(run_path)) # remove file name from example: /igo/sequencers/michelle/211129_MICHELLE_0461_AHMJFJDSX2/CopyComplete.txt

    def copy_samplesheet(completed_run_path):
        print("Copying sample sheet(s) for completed run:" + completed_run_path)

        completed_run = str(os.path.basename(completed_run_path))
        samplesheet = "SampleSheet_" + completed_run + ".csv"

        orig_samplesheet_dir = Variable.get("original_samplesheet_dir", default_var="/pskis34/LIMS/LIMS_SampleSheets/")
        orig_samplesheet =  orig_samplesheet_dir + samplesheet
        dest_samplesheet_dir = Variable.get("destination_samplesheet_dir", default_var="/igo/work/igo/SampleSheetCopies/")
        dest_samplesheet =  dest_samplesheet_dir + samplesheet
        
        cp_command = "cp {} {}".format(orig_samplesheet, dest_samplesheet)
        print("Copying sample sheet:" + cp_command)
        run_cp_task = BashOperator(
            task_id='copy_samplesheet',
            bash_command=cp_command,
        )

        # TODO Split sample sheet for DLP, PED-PEG & 10X

        samplesheets_list = Variable.get("ready_to_demux") # TODO, add default if list is empty
        samplesheets_list.append(dest_samplesheet)
        samplesheets_list.append(completed_run_path)
        Variable.set("ready_to_demux", samplesheets_list)

        trigger_dag_demux = TriggerDagRunOperator(
            task_id='demux_run',
            trigger_dag_id='demux_run',
        )

        trigger_dag_demux
        
    find_completed_runs()

find_completed_runs = find_completed_runs()