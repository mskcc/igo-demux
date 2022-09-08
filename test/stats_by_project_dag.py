from datetime import datetime
from airflow import DAG
from airflow.operators.python import PythonOperator


"""
Airflow DAG to run stats by project giving projectID, recipe parameters
"""
with DAG(
    dag_id="stats_by_project",
    schedule_interval=None,
    start_date=datetime(2022, 1, 1),
    catchup=False,
    tags=["stats_by_project"],
) as dag:

    """ 
    Read the input arguments such as:
    {"project":"13097","recipe":"RNASeq-TruSeqPolyA"}
    """
    def run_stats(ds, **kwargs):
        project = kwargs["params"]["project"]
        recipe = kwargs["params"]["recipe"]
        print("running stats for project {}".format(project))

        # main process of calling stats here

        return "Stats done for project {}".format(project)      

    run_stats_by_project = PythonOperator(
        task_id='stats_by_project',
        python_callable=run_stats,
        provide_context=True,
        email_on_failure=True,
        email='skigodata@mskcc.org',
        dag=dag
    )

    run_stats_by_project