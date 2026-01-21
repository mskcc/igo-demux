from datetime import datetime, timedelta
from airflow import DAG
from airflow.operators.python import PythonOperator
import os
import subprocess
from pathlib import Path
import re
from airflow.utils.email import send_email
from airflow.operators.short_circuit import ShortCircuitOperator

from SampleSheet import SampleSheet
import scripts.cellranger
import scripts.calculate_stats
import scripts.get_sequencing_read_data
import scripts.upload_stats
import Fingerprinting.fingerprinting_dag
import scripts.launch_tcrseq_analysis



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
    dag_id="copy_novaseqx_fastqs_and_analysis",
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


def stats(ds, **kwargs):
    import requests
    from pathlib import Path
    sequencer_path = kwargs["params"]["sequencer_path"]
    samplesheet_path = kwargs["params"]["samplesheet"]
    samplesheet = os.path.basename(samplesheet_path)
    samplesheet_no_ext = os.path.splitext(samplesheet)[0]  # SampleSheet_210331_MICHELLE_0360_BH5KFYDRXY
    sequencer_and_run = samplesheet_no_ext[19:]            # remove 'SampleSheet_210331_'

    done = Path(f"/igo/staging/stats/{sequencer_and_run}/DONE")
    if done.exists():
        return f"Stats already completed for {sequencer_and_run}"

    if ("BONO" in sequencer_and_run) or ("FAUCI2" in sequencer_and_run):
        sample_sheet = SampleSheet(samplesheet_path)
        if "REFERENCE" in samplesheet_path:
            return "No stats for reference "  + samplesheet_path

        if "SC_DLP" in sample_sheet.recipe_set or "SC_SCD-WGS" in sample_sheet.recipe_set:
            scripts.get_total_reads_from_demux.run_DLP(sample_sheet, sequencer_and_run)
            scripts.upload_stats.upload_stats(sequencer_and_run)

            # create the .yaml file for each DLP projects
            output_directory = "/igo/staging/FASTQ/" + sequencer_and_run
            sample_sheet_path = output_directory + "/Reports/SampleSheet.csv"
            stats = output_directory + "/Reports/Demultiplex_Stats.csv"
            run_info = output_directory + "/Reports/RunInfo.xml"
            #python scripts/yaml/generate_metadata.py /igo/delivery/FASTQ/MICHELLE_0480_AH5KTWDSX3_DLP/Project_09443_CT/ \
            #/igo/delivery/FASTQ/MICHELLE_0480_AH5KTWDSX3_DLP/Reports/SampleSheet.csv \
            #/igo/delivery/FASTQ/MICHELLE_0480_AH5KTWDSX3_DLP/Reports/Demultiplex_Stats.csv \
            #/igo/delivery/FASTQ/MICHELLE_0480_AH5KTWDSX3_DLP/Reports/RunInfo.xml \
            #Project_09443_CT \
            #/igo/delivery/FASTQ/MICHELLE_0480_AH5KTWDSX3_DLP/Project_09443_CT/070PP_DLP_UNSORTED_metadata.yaml --revcomp_i5
            for project in sample_sheet.project_set: # such as: Project_09443_CT from the "Sample_Project" column
                fastq_project_dir = output_directory + "/" + project + "/"
                chip_number = get_dlp_chip(sample_sheet, project)
                output_yaml = str(fastq_project_dir) + str(chip_number) + "_metadata.yaml"
                fld_endpoint = "https://igolims.mskcc.org:8443/LimsRest/getDLPFieldMapFile?chipNumber={}".format(chip_number)
                fld_file_path = str(fastq_project_dir) + str(chip_number) + ".fld"
                fld_file = Path(fld_file_path)
                response = requests.get(fld_endpoint, auth = ("pms", "tiagostarbuckslightbike"), verify = False)
                fld_file.write_bytes(response.content)
                python_cmd = "python scripts/yaml/generate_metadata.py " + fastq_project_dir + " " + sample_sheet_path + " " + stats + " " + run_info + " " + fld_file_path + " " + project + " " + output_yaml
                print("Calling DLP generate yaml command: {}".format(python_cmd))
                subprocess.check_output(python_cmd, cwd="/home/igo/shared-single-cell", shell=True)

            return "DLP stats posted and yaml file generated"

        # check if the run is 10X by read length
        atac, use_bases_mask = scripts.get_sequencing_read_data.main(sequencer_path)
        print("read length: {}".format(use_bases_mask))
        if use_bases_mask == [29, 89] or use_bases_mask == [44, 51] or atac:
            # if is atac run, demux is using cellranger mkfastq

            # step 1, generate txt files containing total reads and upload to qc website
            if atac:
                scripts.get_total_reads_from_demux.by_json(sequencer_and_run)
            else:
                scripts.get_total_reads_from_demux.run(sample_sheet, sequencer_and_run)

            scripts.upload_stats.upload_stats(sequencer_and_run)

            # step 2, start cell ranger based on recipe/barcode, check whether multiple fastq files existing
            # trim sequencer_and_run if postfix like _10X exsiting
            sequencer_and_run_prefix = "_".join(sequencer_and_run.split("_")[0:3])
            scripts.cellranger.launch_cellranger_by_sample_sheet(sample_sheet, sequencer_and_run_prefix)

            # add DONE file when all the 10X pipeline finished, -K to wait until finish
            cmd = 'bsub -K -J wait_stats_done_for_{} -w \"ended(create_json___{}*)\" touch /igo/staging/CELLRANGER/{}/DONE'.format(sequencer_and_run_prefix, sequencer_and_run_prefix, sequencer_and_run_prefix)
            print(cmd)
            subprocess.run(cmd, shell=True)

            return "10X Pipeline stats done"

        # if "HumanWholeGenome" in sample_sheet.recipe_set:
        # launch_wgs_stats(sample_sheet, sequencer_and_run)
        # print("DRAGEN WGS stats are running for {}".format(sequencer_and_run))

        # this routine start the DRAGEN and Picard Analysis after the run has demuxed
        scripts.calculate_stats.main(samplesheet_path)

        # use sample sheet to go ahead and start the TCRSeq analysis after fastqs have been created
        scripts.launch_tcrseq_analysis.main(samplesheet_path)

        # add DONE file when all the stats finished, -K to wait until finish
        cmd = 'bsub -K -J wait_stats_done_for_{} -w \"done(uplaodWGSstats{}*)\" touch /igo/staging/stats/{}/DONE'.format(sequencer_and_run, sequencer_and_run, sequencer_and_run)
        print(cmd)
        subprocess.run(cmd, shell=True)

    return "Completed"

def get_dlp_chip(samplesheet, project):
    samplesheet.df_ss_data.reset_index()
    for index, row in samplesheet.df_ss_data.iterrows():
        if (row['Sample_Well'] == 'SC_DLP' or row['Sample_Well'] == 'SC_SCD-WGS') and project == row['Sample_Project']:
            # return chip from 071PP_DLP_UNSORTED_128624A_13_12_IGO_09443_CU_1_1_121
            sample = row['Sample_ID']
            return get_dlp_chip_from_sample_name(sample)

def get_dlp_chip_from_sample_name(sample):
    # given a sample name such as "1_cellcover_4C_128749A_37_42" return the DLP chip ie 128749A
    pattern = r"_[0-9]{6}[A-Z]_"
    match = re.search(pattern, sample)

    if match:
        return match.group()[1:-1]
    else:
        return None

def should_run_stats(**kwargs):
    samplesheet_path = kwargs["params"]["samplesheet"]
    samplesheet = os.path.basename(samplesheet_path)
    run = os.path.splitext(samplesheet)[0][19:]

    stats_done = Path(f"/igo/staging/stats/{run}/DONE")
    cellranger_done = Path(f"/igo/staging/CELLRANGER/{'_'.join(run.split('_')[:3])}/DONE")

    if stats_done.exists() or cellranger_done.exists():
        print(f"⏭ Stats already completed for {run}, skipping.")
        return False

    print(f"▶ Stats not done yet for {run}, continuing.")
    return True

def fingerprinting(ds, **kwargs):
    # read in sample sheet as arguments, filter out projects that need to run fingerprinting
    recipe_list_for_fp = ["PED-PEG", "WGS_Deep", "HC_IMPACT", "HC_IMPACT-Heme", "HC_ACCESS", "WES_Human", "HC_CMOCH"]
    # call fingerprinting_dag.py for each project
    samplesheet_path = kwargs["params"]["samplesheet"]

    if "REFERENCE" in samplesheet_path:
        return "No fingerprinting for reference " + samplesheet_path

    # get project list for running fingerprinting by recipe
    sample_sheet = SampleSheet(samplesheet_path)
    # dictionary of project_ID->genome
    project_genome_dict = pandas.Series(sample_sheet.df_ss_data['Sample_Plate'].values,index=sample_sheet.df_ss_data['Sample_Project']).to_dict()
    project_list_to_run = []
    for project, recipe in sample_sheet.project_dict.items():
        # fingerprinting only support human
        if project_genome_dict[project] == "Human" and recipe in recipe_list_for_fp:
            project_list_to_run.append(project)

    print("Projects need to run fp: {}".format(project_list_to_run))
    if len(project_list_to_run) == 0:
        return "No project need to run fingerprinting"
    else:
        for project in project_list_to_run:
            Fingerprinting.fingerprinting_dag.fingerprint(project[8:])

    return "Completed"

def email_notifier(ds, **kwargs):
    samplesheet_path = kwargs["params"]["samplesheet"]
    samplesheet = os.path.basename(samplesheet_path)
    samplesheet_no_ext = os.path.splitext(samplesheet)[0]  # SampleSheet_210331_MICHELLE_0360_BH5KFYDRXY
    sequencer_and_run = samplesheet_no_ext[19:]            # remove 'SampleSheet_210331_'

    # DLP and reference demux don't have stats, only demux
    if samplesheet_no_ext.endswith("_DLP") or samplesheet_no_ext.endswith("_REFERENCE"):
        content = "{} no stats for DLP | REFERENCE".format(sequencer_and_run)
    else:
        content = "{} stats done".format(sequencer_and_run)

    send_email(
        to=["skigodata@mskcc.org"],
        subject='IGO Cluster Stats Finished',
        html_content=content
    )

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
check_stats_not_done = ShortCircuitOperator(
    task_id="check_stats_not_done",
    python_callable=should_run_stats,
    provide_context=True,
    dag=dag,
)
launch_stats = PythonOperator(
    task_id='launch_stats',
    python_callable=stats,
    provide_context=True,
    email_on_failure=True,
    email='skigodata@mskcc.org',
    dag=dag
)

# step for calling fingerprinting if needed
launch_fingerprinting = PythonOperator(
    task_id='launch_fingerprinting',
    python_callable=fingerprinting,
    provide_context=True,
    email_on_failure=True,
    email='skigodata@mskcc.org',
    dag=dag
)

# step for sending email on stats finish successfully
send_stats_email = PythonOperator(
    task_id='send_stats_email',
    python_callable=email_notifier,
    provide_context=True,
    email_on_failure=True,
    email='skigodata@mskcc.org',
    dag=dag
)

find_runs_task >> copy_runs_task >> check_stats_not_done
check_stats_not_done >> launch_stats >> launch_fingerprinting >> send_stats_email
