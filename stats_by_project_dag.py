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
    {"project_directory":"/igo/staging/FASTQ/RUTH_0141_AH27NGDSX5/Project_13586_B","recipe":"RNASeq_PolyA", "species":"human"}
    """
    def run_stats(ds, **kwargs):
        import scripts.calculate_stats
        import scripts.cellranger
        import subprocess

        project_directory = kwargs["params"]["project_directory"]
        recipe = kwargs["params"]["recipe"]
        species = kwargs["params"]["species"]
        print("running stats for project in this directory {}".format(project_directory))

        # main process of calling stats here
        # let's go ahead and run stats by project
        if "10X_" in recipe:
            scripts.cellranger.lanuch_by_project(project_directory, recipe, species)
        elif "ONT" in recipe:
            cmd = "bsub -J ont_stats -n 8 -M 8 python /igo/work/igo/igo-demux/scripts/ont_stats.py {}".format(project_directory)
            print(cmd)
            subprocess.run(cmd, shell=True)
        else:
            scripts.calculate_stats.main([project_directory, recipe, species])

        return "Stats done for project in this directory {}".format(project_directory)      

    run_stats_by_project = PythonOperator(
        task_id='stats_by_project',
        python_callable=run_stats,
        provide_context=True,
        email_on_failure=True,
        email='skigodata@mskcc.org',
        dag=dag
    )

    run_stats_by_project