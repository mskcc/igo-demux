from datetime import datetime, timedelta
from airflow import DAG
from airflow.operators.python import PythonOperator
import os
import subprocess
from pathlib import Path
import re

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
    dag_id="auto_copy_after_copycomplete",
    description="Automatically copy analysis and FASTQ files after CopyComplete.txt appears.",
    schedule_interval='@hourly',
    start_date=datetime(2025, 1, 1),
    catchup=False,
    tags=["igo", "sequencer", "copy"],
)

copy_script = "/igo/work/igo/igo-demux/scripts/move_novaseqx_analysis_files.py"

SEQUENCERS = ["bono", "fauci2"]  # Extend if needed
SEQUENCER_BASE = "/igo/sequencers"


def find_runs_ready_for_copy(**context):
    """
    Scans sequencer directories for run folders and checks if:
        Analysis/1/CopyComplete.txt exists.
    Returns list of full run paths ready for copying.
    """

    ready_runs = []

    run_pattern = re.compile(r"^\d{6}_(?P<seq>[A-Za-z0-9]+)_(?P<runid>.+)$")

    for seq in SEQUENCERS:
        seq_dir = Path(SEQUENCER_BASE) / seq

        if not seq_dir.exists():
            continue

        for run_folder in seq_dir.iterdir():
            if not run_folder.is_dir():
                continue

            m = run_pattern.match(run_folder.name)
            if not m:
                continue

            # CopyComplete checking path
            cc_path = run_folder / "Analysis" / "1" / "CopyComplete.txt"

            if cc_path.exists():
                print(f"✅ Found CopyComplete.txt: {cc_path}")
                ready_runs.append(str(run_folder))
            else:
                print(f"⏳ Not complete: {cc_path}")

    print(f"Total runs ready for copying: {len(ready_runs)}")
    return ready_runs


def run_copy_script(**context):
    """
    Executes the copy script one time.
    The script itself discovers all sequencers and runs independently.
    """
    ready_runs = context["ti"].xcom_pull(task_ids="find_runs_ready_for_copy")

    if not ready_runs:
        print("No completed runs detected. Nothing to copy.")
        return

    print("Starting copy script...")

    try:
        subprocess.run(
            ["/home/igo/miniconda_airflow/bin/python3.9", copy_script],
            check=True,
        )
        print("✅ Copy script executed successfully.")
    except subprocess.CalledProcessError as e:
        print(f"❌ Copy script failed: {e}")


find_runs_task = PythonOperator(
    task_id="find_runs_ready_for_copy",
    python_callable=find_runs_ready_for_copy,
    provide_context=True,
    dag=dag,
)

copy_runs_task = PythonOperator(
    task_id="run_copy_script",
    python_callable=run_copy_script,
    provide_context=True,
    dag=dag,
)

find_runs_task >> copy_runs_task