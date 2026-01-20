from datetime import datetime
from airflow import DAG
from airflow.providers.standard.operators.python import PythonOperator
import scripts.cellranger

"""
Airflow DAG to generate bams only by project giving projectID, recipe parameters
"""
with DAG(
    dag_id="bams_only_by_project",
    schedule_interval=None,
    start_date=datetime(2023, 1, 1),
    catchup=False,
    tags=["bams_only_by_project"],
) as dag:

    """ 
    Read the input arguments such as:
    {"project_directory":"/igo/staging/FASTQ/RUTH_0141_AH27NGDSX5/Project_13586_B","recipe":"RNASeq_PolyA", "species":"human"}
    """
    def generate_bams(ds, **kwargs):
        project_directory = kwargs["params"]["project_directory"]
        recipe = kwargs["params"]["recipe"]
        species = kwargs["params"]["species"]
        print("creating only the bams for project in this directory {}".format(project_directory))

        # main process of calling stats here
        # let's go ahead and run stats by project
        if "10X_" in recipe:
            scripts.cellranger.lanuch_by_project(project_directory, recipe, species)
        else:
            scripts.alignment_only.main([project_directory, recipe, species])

        return "Bams generated for project in this directory {}".format(project_directory)      

    generate_bams_only_by_project = PythonOperator(
        task_id='bams_only_by_project',
        python_callable=generate_bams,
        provide_context=True,
        email_on_failure=True,
        email='skigodata@mskcc.org',
        dag=dag
    )

    generate_bams_only_by_project