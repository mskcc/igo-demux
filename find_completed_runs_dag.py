import sequencer
import os
import datetime
import json

from airflow.providers.http.operators.http import SimpleHttpOperator
from airflow.operators.bash import BashOperator
from airflow import DAG
from airflow.models import Variable

# defines the list of all sequencers and for each sequencer 1) name 2) location it writes runs to and 3) the last file the sequencer writes when a run is completed to signal demux can begin
sequencers = {"sequencers":[{"name":"ayyan","path":"/igo/sequencers/ayyan","last_file":"RTAComplete.txt"},{"name":"diana","path":"/igo/sequencers/diana","last_file":"CopyComplete.txt"},{"name":"michelle","path":"/igo/sequencers/michelle","last_file":"CopyComplete.txt"},{"name":"ruth","path":"/igo/sequencers/ruth","last_file":"CopyComplete.txt"},{"name":"johnsawyers","path":"/igo/sequencers/johnsawyers","last_file":"RTAComplete.txt"},{"name":"pepe","path":"/igo/sequencers/pepe/output","last_file":"CopyComplete.txt"},{"name":"scott","path":"/igo/sequencers/scott","last_file":"RunCompletionStatus.xml"}]}

with DAG(
    dag_id='find_completed_runs', 
    schedule_interval=None, 
    start_date=datetime.datetime(2021, 1, 1), 
    catchup=False,
    tags=["find_completed_runs"],
) as dag:
    """
    #### Find Completed Runs
    Find recently completed runs by looking for the last file written by the sequencers then copy the sample sheet for the completed run to the list of completed runs variable
    so a later task can start the demux for the sample sheet.
    """
    print("Processing sequencer list: {}".format(sequencers))
        
    time_to_search = Variable.get("completed_run_search_interval_mins", default_var=30)
    print("Searching for runs completed in the last {} minutes, variable completed_run_search_interval_mins".format(time_to_search))
    completed_runs_path = sequencer.find_completed_runs(sequencers, time_to_search)


    for run_path in completed_runs_path:
        # remove file name from path, for example: /igo/sequencers/michelle/211129_MICHELLE_0461_AHMJFJDSX2/CopyComplete.txt
        completed_run_path = os.path.dirname(run_path)
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
            task_id='copy_samplesheet_'+completed_run, # make the task_id unique for each run
            bash_command=cp_command,
        )

        # TODO Split sample sheet for DLP, PED-PEG & 10X
        
        demux_dict = {}
        demux_dict['samplesheet'] = dest_samplesheet
        demux_dict['sequencer_path'] = completed_run_path
        demux_args_json = json.dumps(demux_dict)

        future = '"'+(datetime.datetime.now() + datetime.timedelta(seconds=60)).strftime("%Y-%m-%dT%H:%M:%SZ")+'"'
        # Airflow required arguments to trigger a dag - execution date and conf arguments
        dag_json = '{"execution_date": '+future +',"conf": '+demux_args_json+'}'
        print("Calling demux with execution time and args:" + dag_json)

        trigger_dag_demux = SimpleHttpOperator(
            task_id="start_demux_"+completed_run,
            http_conn_id='airflow-api',
            endpoint='api/v1/dags/demux_run/dagRuns',
            method='POST',
            headers={'Content-Type': 'application/json'},
            data=dag_json,
        )

        run_cp_task >> trigger_dag_demux