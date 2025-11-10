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

copy_script = "/igo-demux/scripts/move_novaseqx_analysis_files.py"


def get_latest_completed_runs(**context):
    """
    Fetch list of completed runs from the latest successful
    find_completed_runs_dag execution via XCom.
    """
    session = context["session"]
    dag_runs = (
        session.query(DagRun)
        .filter(DagRun.dag_id == "find_completed_runs_dag", DagRun.state == State.SUCCESS)
        .order_by(DagRun.execution_date.desc())
        .all()
    )

    if not dag_runs:
        print("⚠️ No successful find_completed_runs_dag runs found.")
        return []

    latest = dag_runs[0]
    ti = latest.get_task_instances()
    completed_runs = None

    for t in ti:
        if t.task_id == "find_completed_runs":
            completed_runs = t.xcom_pull(task_ids="find_completed_runs", key="completed_runs")
            break

    if not completed_runs:
        print("⚠️ No completed_runs XCom found in the last successful run.")
        return []

    print(f"✅ Found completed runs: {completed_runs}")
    return completed_runs


def copy_runs_to_staging(**context):
    """
    Execute move_novaseqx_analysis_files.py for each completed run.
    """
    completed_runs = context["ti"].xcom_pull(task_ids="get_latest_completed_runs")
    if not completed_runs:
        print("No completed runs found.")
        return

    for run_id in completed_runs:
        print(f"Copying run {run_id} to staging...")
        try:
            subprocess.run(["/home/igo/miniconda_airflow/bin/python3.9 ", copy_script, run_id], check=True)
            print(f"✅ Successfully copied {run_id}")
        except subprocess.CalledProcessError as e:
            print(f"❌ Copy failed for {run_id}: {e}")


get_runs_task = PythonOperator(
    task_id="get_latest_completed_runs",
    python_callable=get_latest_completed_runs,
    provide_context=True,
    dag=dag,
)

copy_runs_task = PythonOperator(
    task_id="copy_runs_to_staging",
    python_callable=copy_runs_to_staging,
    provide_context=True,
    dag=dag,
)

get_runs_task >> copy_runs_task