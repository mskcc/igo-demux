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
        import scripts.cellranger_multi_v9
        import os
        import scripts.get_total_reads_from_demux

        project_directory = kwargs["params"]["project_directory"]
        recipe = kwargs["params"]["recipe"]
        species = kwargs["params"]["species"]
        print("running stats for project in this directory {}".format(project_directory))

        # main process of calling stats here
        # let's go ahead and run stats by project
        # add multi process, use recipe as 10X_multi, the project folder has to be gene expression project folder
        project_id = project_directory.split("/")[-1]

        if recipe == "10X_multi":
            scripts.cellranger_multi_v9.main(project_directory, species)

        elif "SC_Chromium" in recipe or "ST_Visium" in recipe:
            scripts.cellranger.launch_cellranger_by_project_location(project_directory, recipe, species)
        elif "Nanopore" in recipe:
            os.chdir("/igo/staging/PIPELINE")
            cmd = "bsub -J ont_stats_{} -o ont_stats_{}.out -n 16 -M 16 /igo/work/nabors/tools/venvpy3/bin/python /igo/work/igo/igo-demux/scripts/ont_stats.py {}".format(project_id, project_id, project_directory)
            print(cmd)
            subprocess.run(cmd, shell=True)
        elif recipe == "demux_stats":
            scripts.get_total_reads_from_demux.by_project_location(project_directory)
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