from datetime import datetime, timedelta
from airflow import DAG
from airflow.operators.python import PythonOperator
from airflow.models import DagRun
from airflow.utils.state import State
import subprocess

"""
This DAG triggers automatically after find_completed_runs_dag identifies completed runs.
It uses the existing find_completed_runs_dag output instead of duplicating its logic.
"""

default_args = {
    "owner": "igo",
    "depends_on_past": False,
    "email": ["skigodata@mskcc.org"],
    "email_on_failure": True,
    "email_on_retry": False,
    "retries": 1,
    "retry_delay": timedelta(minutes=10),
}

dag = DAG(
    dag_id="auto_copy_after_find_runs",
    description="Automatically copy sequencing analysis files after find_completed_runs_dag completes.",
    schedule_interval='@hourly',
    start_date=datetime(2025, 1, 1),
    catchup=False,
    tags=["igo", "sequencer", "copy"],
)

copy_script = "/igo/work/igo/igo-demux/scripts/move_novaseqx_analysis_files.py"
completed_runs_file = "/igo/sequencers/completed_runs.json"


def get_completed_runs_from_file(**context):
    """
    Reads completed runs from a JSON file written by the find_completed_runs DAG.
    Example file content:
    {
        "timestamp": "2025-11-06T12:00:00",
        "completed_runs": ["FAUCI2_0073_B232GYLLT3", "BONO_0123_AHJKLXYZ"]
    }
    """
    if not os.path.exists(completed_runs_file):
        print(f"Completed runs file not found: {completed_runs_file}")
        return []

    try:
        with open(completed_runs_file, "r") as f:
            data = json.load(f)
        completed_runs = data.get("completed_runs", [])
        print(f"✅ Found {len(completed_runs)} completed runs: {completed_runs}")
        return completed_runs
    except Exception as e:
        print(f"❌ Error reading completed runs file: {e}")
        return []


def copy_runs_to_staging(**context):
    """
    Runs the copy script for each completed run ID found in the file.
    """
    ti = context["ti"]
    completed_runs = ti.xcom_pull(task_ids="get_completed_runs_from_file")

    if not completed_runs:
        print("No completed runs found.")
        return

    for run_id in completed_runs:
        print(f"Copying run {run_id} to staging...")
        try:
            subprocess.run(
                ["/home/igo/miniconda_airflow/bin/python3.9", copy_script, run_id],
                check=True,
            )
            print(f"✅ Successfully copied {run_id}")
        except subprocess.CalledProcessError as e:
            print(f"❌ Copy failed for {run_id}: {e}")

def clear_completed_runs_file(**context):
    """
    Clears the contents of /igo/sequencers/completed_runs.json after successful copy.
    """
    if not os.path.exists(completed_runs_file):
        print("No completed_runs.json file found to clear.")
        return

    try:
        with open(completed_runs_file, "w") as f:
            json.dump({"timestamp": datetime.now().isoformat(), "completed_runs": []}, f, indent=2)
        print(f"Cleared completed runs in {completed_runs_file}")
    except Exception as e:
        print(f"Failed to clear completed_runs.json: {e}")


get_runs_task = PythonOperator(
    task_id="get_completed_runs_from_file",
    python_callable=get_completed_runs_from_file,
    provide_context=True,
    dag=dag,
)

copy_runs_task = PythonOperator(
    task_id="copy_runs_to_staging",
    python_callable=copy_runs_to_staging,
    provide_context=True,
    dag=dag,
)

clear_file_task = PythonOperator(
    task_id="clear_completed_runs_file",
    python_callable=clear_completed_runs_file,
    provide_context=True,
    dag=dag,
)

get_runs_task >> copy_runs_task >> clear_file_task