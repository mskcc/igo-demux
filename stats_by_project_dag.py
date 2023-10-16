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
        import scripts.cellranger_multi
        import os

        project_directory = kwargs["params"]["project_directory"]
        recipe = kwargs["params"]["recipe"]
        species = kwargs["params"]["species"]
        print("running stats for project in this directory {}".format(project_directory))

        # main process of calling stats here
        # let's go ahead and run stats by project
        # add multi process, use recipe as 10X_multi, the project folder has to be gene expression project folder
        if recipe == "10X_multi":
            project_id = project_directory.split("/")[-1]
            # copy the multi config from shared drive to cluster
            cmd = "cp -R {}{} {}".format(scripts.cellranger_multi.ORIGIN_DRIVE_LOCATION, project_id[8:], scripts.cellranger_multi.DRIVE_LOCATION)
            print(cmd)
            subprocess.run(cmd, shell=True)
            os.chdir(scripts.cellranger_multi.STATS_AREA)
            # gather sample set info from LIMS for each sample
            sample_list_ori = os.listdir(project_directory)
            sample_list = []
            for sample in sample_list_ori:
                # remove Sample_ prefix
                sample_list.append(sample[7:])
            for sample in sample_list:
                sample_set = scripts.cellranger_multi.gather_sample_set_info(sample)
                cmd = "bsub -J {}_{}_multi -o {}_{}_multi.out /igo/work/nabors/tools/venvpy3/bin/python /igo/work/igo/igo-demux/scripts/cellranger_multi.py ".format(project_id, sample, project_id, sample)
                for key, value in sample_set.items():
                    if value is not None:
                        cmd = cmd + "-{}={} ".format(key, value)
                cmd = cmd + "-genome={}".format(species)
                print(cmd)
                subprocess.run(cmd, shell=True)

        elif "10X_" in recipe:
            scripts.cellranger.lanuch_by_project(project_directory, recipe, species)
        elif "ONT" in recipe:
            cmd = "bsub -J ont_stats -n 16 -M 16 /igo/work/nabors/tools/venvpy3/bin/python /igo/work/igo/igo-demux/scripts/ont_stats.py {}".format(project_directory)
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