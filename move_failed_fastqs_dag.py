from datetime import datetime

import scripts.move_failed_fastqs

from airflow import DAG
from airflow.providers.standard.operators.python import PythonOperator

"""
DAG to wraps the Python script which removes failed fastqs and .bams
"""
with DAG(
    dag_id="move_failed_fastqs",
    schedule=None,
    start_date=datetime(2022, 1, 1),
    catchup=False,
    tags=["move_failed_fastqs"],
) as dag:

    """ 
    Read the input arguments such as:

    {"igo_id":"13038_E_10","run":"MICHELLE_0521_BHHWTKDSX3"}
    """
    def move_fastqs(ds, **kwargs):
        igo_id = kwargs["params"]["igo_id"]
        run = kwargs["params"]["run"]
        print("Moving failed fastqs for sample {} and run {}".format(igo_id, run))

        result = scripts.move_failed_fastqs.move_failed_fastqs(igo_id, run)
        return result      

    move_failed_fastqs = PythonOperator(
        task_id='move_failed_fastqs',
        python_callable=move_fastqs,
        provide_context=True,
        email_on_failure=True,
        email='skigodata@mskcc.org',
        dag=dag
    )

    move_failed_fastqs
