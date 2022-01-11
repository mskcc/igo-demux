import sequencer
import os
import datetime
import json
import shutil

from SampleSheet import SampleSheet

from airflow.providers.http.operators.http import SimpleHttpOperator
from airflow.operators.bash import BashOperator
from airflow import DAG
from airflow.models import Variable
from airflow.decorators import task
from airflow.operators.email_operator import EmailOperator

# defines the list of all sequencers and for each sequencer 1) name 2) location it writes runs to and 3) the last file the sequencer writes when a run is completed to signal demux can begin
sequencers = {"sequencers":[{"name":"ayyan","path":"/igo/sequencers/ayyan","last_file":"RTAComplete.txt"},{"name":"diana","path":"/igo/sequencers/diana","last_file":"CopyComplete.txt"},{"name":"michelle","path":"/igo/sequencers/michelle","last_file":"CopyComplete.txt"},{"name":"ruth","path":"/igo/sequencers/ruth","last_file":"CopyComplete.txt"},{"name":"johnsawyers","path":"/igo/sequencers/johnsawyers","last_file":"RTAComplete.txt"},{"name":"pepe","path":"/igo/sequencers/pepe/output","last_file":"CopyComplete.txt"},{"name":"scott","path":"/igo/sequencers/scott","last_file":"RunCompletionStatus.xml"}]}

with DAG(
    dag_id='find_completed_runs', 
    schedule_interval='@hourly', 
    start_date=datetime.datetime(2022, 1, 1), 
    catchup=False,
    tags=["find_completed_runs"],
) as dag:
    """
    #### Find Completed Runs
    Find recently completed runs by looking for the last file written by the sequencers then copy the sample sheet for the completed run to the list of completed runs variable
    so a later task can start the demux for the sample sheet.
    """
    
    demux_run_path = Variable.get("demux_run_path", default_var="")
    if len(demux_run_path) != 0 :
        completed_runs_path = [demux_run_path + "/RTAComplete.txt"]
        print("Demuxing single run only {}".format(completed_runs_path))
    else:
        print("Processing sequencer list: {}".format(sequencers))
        time_to_search = Variable.get("completed_run_search_interval_mins", default_var=60) # should match schedule_interval above
        print("Searching for runs completed in the last {} minutes, variable completed_run_search_interval_mins".format(time_to_search))
        completed_runs_path = sequencer.find_completed_runs(sequencers, time_to_search)

    for run_path in completed_runs_path:
        # remove file name from path, for example: /igo/sequencers/michelle/211129_MICHELLE_0461_AHMJFJDSX2/CopyComplete.txt
        completed_run_path = os.path.dirname(run_path)
        print("Copying sample sheet(s) for completed run:" + completed_run_path)
        completed_run = str(os.path.basename(completed_run_path)) # ex: 211129_MICHELLE_0461_AHMJFJDSX2
        run_name_only = completed_run[7:]  # ex: MICHELLE_0461_AHMJFJDSX2
        samplesheet = "SampleSheet_" + completed_run + ".csv"

        orig_samplesheet_dir = Variable.get("original_samplesheet_dir", default_var="/pskis34/LIMS/LIMS_SampleSheets/")
        orig_samplesheet =  orig_samplesheet_dir + samplesheet
        dest_samplesheet_dir = Variable.get("destination_samplesheet_dir", default_var="/igo/work/igo/SampleSheetCopies/")
        dest_samplesheet =  dest_samplesheet_dir + samplesheet

        ss_orig = SampleSheet(orig_samplesheet)
        ss_orig.remove_sample_prefix()
        ss_orig.path = dest_samplesheet
        
        ss_list = ss_orig.split_sample_sheet()

        email_to = Variable.get("email_to", default_var="skigodata@mskcc.org")
        send_email = EmailOperator(
            task_id='send_email',
            to=email_to,
            subject='IGO Cluster New Run Sent for Demuxing',
            html_content=" <h3>{}</h3> sent to DRAGEN split into {} sample sheets".format(run_name_only, len(ss_list)),
            dag=dag
        )

        counter = 0
        for samplesheet in ss_list:
            counter = counter + 1
            print("Saving sample sheet {}".format(samplesheet.path))
            samplesheet.write_csv()

            demux_dict = {}
            demux_dict['samplesheet'] = samplesheet.path
            demux_dict['sequencer_path'] = completed_run_path
            demux_args_json = json.dumps(demux_dict)

            future = '"'+(datetime.datetime.now() + datetime.timedelta(seconds=(counter*10))).strftime("%Y-%m-%dT%H:%M:%SZ")+'"'
            # Airflow required arguments to trigger a dag - execution date and conf arguments
            dag_json = '{"execution_date": '+future +',"conf": '+demux_args_json+'}'
            print("Calling demux with execution time and args:" + dag_json)

            trigger_dag_demux = SimpleHttpOperator(
                task_id="demux_"+str(os.path.basename(samplesheet.path)).replace(".csv", "").replace("SampleSheet_",""),
                http_conn_id='airflow-api',
                endpoint='api/v1/dags/demux_run/dagRuns',
                method='POST',
                headers={'Content-Type': 'application/json'},
                data=dag_json,
            )

            # first steps of pipeline are in pure Python that copy the sample sheet from /pskis34
            send_email >> trigger_dag_demux