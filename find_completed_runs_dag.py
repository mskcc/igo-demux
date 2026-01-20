import sequencer
import os
import datetime
import json
import time

from SampleSheet import SampleSheet

from airflow.providers.http.operators.http import HttpOperator
from airflow.operators.bash import BashOperator
from airflow import DAG
from airflow.models import Variable
from airflow.decorators import task
from airflow.providers.smtp.operators.email import EmailOperator
from airflow.sdk.dag import dag; from airflow.sdk.task import task

# defines the list of all sequencers and for each sequencer 1) name 2) location it writes runs to and 3) the last file the sequencer writes when a run is completed to signal demux can begin
sequencers = {
   "sequencers": [
      {
         "name": "ayyan",
         "path": "/igo/sequencers/ayyan",
         "last_file": "CopyComplete.txt"
      },
      {
         "name": "diana",
         "path": "/igo/sequencers/diana",
         "last_file": "CopyComplete.txt"
      },
      {
         "name": "fauci",
         "path": "/igo/sequencers/fauci",
         "last_file": "CopyComplete.txt"
      },
      {
         "name": "bono",
         "path": "/igo/sequencers/bono",
         "last_file": "CopyComplete.txt"
      },
      {
         "name": "fauci2",
         "path": "/igo/sequencers/fauci2",
         "last_file": "CopyComplete.txt"
      },
      {
         "name": "ruth",
         "path": "/igo/sequencers/ruth",
         "last_file": "CopyComplete.txt"
      },
      {
         "name": "johnsawyers",
         "path": "/igo/sequencers/johnsawyers",
         "last_file": "CopyComplete.txt"
      },
      {
         "name": "pepe",
         "path": "/igo/sequencers/pepe/output",
         "last_file": "CopyComplete.txt"
      },
      {
         "name": "amelie",
         "path": "/igo/sequencers/amelie/output",
         "last_file": "CopyComplete.txt"
      },
   ]
}

"""
Find recently completed runs by looking for the last file written by the sequencers,
then split and copy the sample sheet for the completed run and launch the demux task
"""
with DAG(
    dag_id='find_completed_runs',
    schedule_interval='@hourly',
    start_date=datetime.datetime(2022, 1, 1),
    catchup=False,
    tags=["find_completed_runs"],
) as dag:

   completed_runs_path = list()

   # TODO - Consider making separate DAG to trigger 1 specific demux?
   demux_special = "/igo/sequencers/run_to_demux.txt"
   if os.path.exists(demux_special):
      run_to_demux_file = open("/igo/sequencers/run_to_demux.txt", "r")
      demux = run_to_demux_file.readline().strip()
      run_to_demux_file.close()
      if len(demux) > 0:
         completed_runs_path.append(demux + "/RTAComplete.txt")
      # TODO check if path exists, if not then done

   if len(completed_runs_path) == 0:
      print("Processing sequencer list: {}".format(sequencers))
      # should match schedule_interval above
      time_to_search = Variable.get("completed_run_search_interval_mins", default_var=60)
      print("Searching for runs completed in the last {} minutes".format(time_to_search))
      completed_runs_path = sequencer.find_completed_runs(sequencers, time_to_search)

   for run_path in completed_runs_path:
      # remove file name from path, for example: /igo/sequencers/michelle/211129_MICHELLE_0461_AHMJFJDSX2/CopyComplete.txt
      completed_run_path = os.path.dirname(run_path)
      print("Copying sample sheet(s) for completed run:" + completed_run_path)
      # ex: 211129_MICHELLE_0461_AHMJFJDSX2
      completed_run = str(os.path.basename(completed_run_path))
      run_name_only = completed_run[7:]  # ex: MICHELLE_0461_AHMJFJDSX2
      samplesheet = "SampleSheet_" + completed_run + ".csv"

      orig_samplesheet_dir = Variable.get("original_samplesheet_dir", default_var="/rtssdc/mohibullahlab/LIMS/LIMS_SampleSheets/")
      orig_samplesheet = orig_samplesheet_dir + samplesheet
      dest_samplesheet_dir = Variable.get("destination_samplesheet_dir", default_var="/igo/work/igo/SampleSheetCopies/")
      dest_samplesheet = dest_samplesheet_dir + samplesheet

      print("Reading the LIMS sample sheet {}".format(orig_samplesheet))
      ss_orig = SampleSheet(orig_samplesheet)
      ss_orig.path = dest_samplesheet

      ss_list = ss_orig.split_sample_sheet()
      ss_list_str = "\n" # format string to be readable in Data Team emails
      for sheet in ss_list:
         ss_list_str += sheet.path + "\n"

      email_to = Variable.get("email_to", default_var="skigodata@mskcc.org")
      send_demux_email = EmailOperator(
         task_id='send_demux_email'+run_name_only,
            to=email_to,
            subject='IGO Cluster New Run Sent for Demuxing',
            html_content="<h3>{}</h3> sent to DRAGEN split into: {} ".format(run_name_only, ss_list_str),
            dag=dag
        )

      counter = 0
      for samplesheet in ss_list:
         counter = counter + 1
         print("Saving sample sheet {}".format(samplesheet.path))
         samplesheet.write_csv()

         demux_dict = {}
         demux_dict['dragen_demux'] = 'False'
         demux_dict['samplesheet'] = samplesheet.path
         demux_dict['sequencer_path'] = completed_run_path
         demux_args_json = json.dumps(demux_dict)

         exec_time = (datetime.datetime.now() + datetime.timedelta(seconds=(counter))).strftime("%Y-%m-%dT%H:%M:%SZ")
         future = '"'+exec_time+'"'
         # Airflow required arguments to trigger a dag - execution date and conf arguments
         dag_json = '{"execution_date": ' + future + ',"conf": '+demux_args_json+'}'
         print("Calling demux with execution time and args:" + dag_json)

         trigger_dag_demux = SimpleHttpOperator(
            task_id="demux_" +
            str(os.path.basename(samplesheet.path)).replace(".csv", "").replace("SampleSheet_", ""),
               http_conn_id='airflow-api',
               endpoint='api/v1/dags/demux_run/dagRuns',
               method='POST',
               headers={'Content-Type': 'application/json'},
               data=dag_json,
            )

         send_demux_email >> trigger_dag_demux

      # Airflow can only call the demux endpoint once at a specific time, do not allow demux endpoint calls to overlap
      time.sleep(5)