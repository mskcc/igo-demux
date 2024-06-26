from datetime import datetime

import scripts.deliver_pipeline

from airflow import DAG
from airflow.operators.python import PythonOperator


"""
Airflow DAG to call the deliver_pipeline.py code.
"""
with DAG(
    dag_id="deliver_pipeline",
    schedule_interval=None,
    start_date=datetime(2022, 1, 1),
    catchup=False,
    tags=["deliver_pipeline"],
) as dag:

    """ 
    Read the input arguments such as:

    {"project":"13097","pi":"abdelwao","recipe":"RNASeq-TruSeqPolyA"}
    """
    def deliver(ds, **kwargs):
        project = kwargs["params"]["project"]
        pi = kwargs["params"]["pi"]
        # recipe here is actually request name
        recipe = kwargs["params"]["recipe"]
        print("Delivering the pipeline output and/or .bams for {} {} {}".format(project, pi, recipe))

        result = scripts.deliver_pipeline.deliver_pipeline_output(project, pi, recipe)
        return result      

    deliver_pipeline_output = PythonOperator(
        task_id='deliver_pipeline',
        python_callable=deliver,
        provide_context=True,
        email_on_failure=True,
        email='skigodata@mskcc.org',
        dag=dag
    )

    deliver_pipeline_output