import os
from re import sub
import subprocess
from datetime import datetime, timedelta
from SampleSheet import SampleSheet
import scripts.organise_fastq_split_by_lane
import pandas

from airflow import DAG
from airflow.operators.python import PythonOperator
from airflow.operators.email_operator import EmailOperator
from airflow.models import Variable


"""
    Runs the demux such as: 
    bsub -n48 -q dragen /opt/edico/bin/dragen --bcl-conversion-only true --bcl-sampleproject-subdirectories true --force
      --bcl-input-directory /igo/sequencers/johnsawyers/211108_JOHNSAWYERS_0312 
      --output-directory /igo/staging/FASTQ/JOHNSAWYERS_0312 
      --sample-sheet /home/igo/DividedSampleSheets/SampleSheet_211108_JOHNSAWYERS_0312.csv 
"""
with DAG(
    dag_id="demux_run",
    schedule_interval=None,
    start_date=datetime(2022, 1, 1),
    catchup=False,
    tags=["demux_run"],
) as dag:

    """ 
    Read the input arguments such as:

    'params': {'samplesheet': '/igo/work/igo/SampleSheetCopies/SampleSheet_211206_JOHNSAWYERS_0317_000000000-K3LFK.csv',
             'sequencer_path': '/igo/sequencers/johnsawyers/211206_JOHNSAWYERS_0317_000000000-K3LFK'},
    """
    def demux(ds, **kwargs):
        sequencer_path = kwargs["params"]["sequencer_path"]
        samplesheet_path = kwargs["params"]["samplesheet"]
        print("Starting demux {} {}".format(sequencer_path, samplesheet_path))

        samplesheet = os.path.basename(samplesheet_path)
        samplesheet_no_ext = os.path.splitext(samplesheet)[0]  # SampleSheet_210331_MICHELLE_0360_BH5KFYDRXY
        sequencer_and_run = samplesheet_no_ext[19:]            # remove 'SampleSheet_210331_'

        sample_sheet = SampleSheet(samplesheet_path)
        
        is_DLP = False
        if "DLP" in sample_sheet.recipe_set:
            is_DLP = True
        is_WGS = False
        if "HumanWholeGenome" in sample_sheet.recipe_set:
            is_WGS = True
        is_10X = False
        if len(sample_sheet.barcode_list_10X) > 0:
            is_10X = True

        if is_DLP:
            output_directory = "/igo/staging/FASTQ/" + sequencer_and_run + "_DLP"
        if is_WGS:
            output_directory = "/igo/staging/FASTQ/" + sequencer_and_run + "_WGS"
        if is_10X:
            output_directory = "/igo/staging/FASTQ/" + sequencer_and_run + "_10X"
        else:
            output_directory = "/igo/staging/FASTQ/" + sequencer_and_run
        
        is_DRAGEN_demux = True
        
        if is_10X:
            is_DRAGEN_demux = False
            # TODO for 10X build correct mkfastq command, special 10X barcodes can't go to dragen
            print("mkfastq command is not yet supported, this must be launched at the command line by the Data Team")
        else:
            # DLP can demux with the default command as long as the [Settings] have 'NoLaneSplitting,true'
            # -K - wait for the job to complete
            bsub_command = "bsub -K -n48 -q dragen -e /igo/work/igo/igo-demux/logs/demux.log -o /igo/work/igo/igo-demux/logs/demux.log "
            command = bsub_command + "/opt/edico/bin/dragen --bcl-conversion-only true --bcl-only-matched-reads true --force --bcl-sampleproject-subdirectories true --bcl-input-directory \'{}\' --output-directory \'{}\' --sample-sheet \'{}\'".format(
            sequencer_path, output_directory, samplesheet_path)
            print("Running demux command: " + command)
            subprocess.run(command, shell=True, check=True)

        # if the demux was successful:
        if is_DRAGEN_demux and not is_DLP:
            print("Adding sample sub-folders to the DRAGEN demux.")
            scripts.organise_fastq_split_by_lane.create_fastq_folders(output_directory)
            scripts.organise_fastq_split_by_lane.correct_fastq_list_csv(output_directory+"/Reports")

        # Call CopyIlluminaReports.sh /igo/staging/FASTQ/RUTH_0066_BHTJ33DRXY
        copy_reports_cmd = "/igo/work/igo/igo-demux/CopyIlluminaReports.sh /igo/staging/FASTQ/" + sequencer_and_run
        print("Running command to copy demux reports: " + copy_reports_cmd)
        subprocess.run(copy_reports_cmd, shell=True, check=True)
        
        # for DLP projects create the .yaml file
        if is_DLP:
            # example: make create-metadata-yaml ss=/igo/home/igo/DividedSampleSheets/SampleSheet_211022_DIANA_0415_AHMJWMDSX2_DLP.csv prj=Project_09443_CI project_path=/igo/staging/FASTQ/DIANA_0415_AHMJWMDSX2_DLP/Project_09443_CI
            for project in sample_sheet.project_set:
                fastq_project_dir = output_directory + "/" + project
                make_command = "make create-metadata-yaml ss={} prj={} project_path={}".format(samplesheet_path, project, fastq_project_dir)
                print("Calling DLP make command: {}".format(make_command))
                subprocess.check_output(make_command, cwd="/home/igo/shared-single-cell", shell=True)

        return command


    def stats(ds, **kwargs):
        sequencer_path = kwargs["params"]["sequencer_path"]
        samplesheet_path = kwargs["params"]["samplesheet"]
        samplesheet = os.path.basename(samplesheet_path)
        samplesheet_no_ext = os.path.splitext(samplesheet)[0]  # SampleSheet_210331_MICHELLE_0360_BH5KFYDRXY
        sequencer_and_run = samplesheet_no_ext[19:]            # remove 'SampleSheet_210331_'

        sample_sheet = SampleSheet(samplesheet_path)
        if "REFERENCE" in samplesheet_path:
            return "No stats for reference "  + samplesheet_path

        if "DLP" in sample_sheet.recipe_set:
           return "No DLP stats"
        
        if "HumanWholeGenome" in sample_sheet.recipe_set:
            launch_wgs_stats(sample_sheet, sequencer_and_run)
            return "DRAGEN WGS stats are running for " + sequencer_and_run

        if any("10X_" in s for s in sample_sheet.recipe_set):
            # TODO write code to launch 10X pipelines
            return "Not launching 10X Pipeline"

        launch_stats_via_bash_script(sample_sheet, sequencer_and_run)

        return "Completed"

    demux_run = PythonOperator(
        task_id='start_the_demux',
        python_callable=demux,
        provide_context=True,
        email_on_failure=True,
        email='skigodata@mskcc.org',
        dag=dag
    )

    launch_stats = PythonOperator(
        task_id='launch_stats',
        python_callable=stats,
        provide_context=True,
        email_on_failure=True,
        email='skigodata@mskcc.org',
        dag=dag
    )

    demux_run >> launch_stats


    def launch_wgs_stats(sample_sheet, sequencer_and_run):
        stats_path = "/igo/staging/stats/" + sequencer_and_run
        if not os.path.exists(stats_path):
            os.makedirs(stats_path)
            print("Created the stats directory: {}".format(stats_path))
        
        cmds_dragen = build_dragen_cmds(sample_sheet, sequencer_and_run)
        for cmd in cmds_dragen:
            subprocess.run(cmd, shell=True)
        
        cmds_bwamem2 = build_bwamem2_cmds(sample_sheet, sequencer_and_run)
        for cmd in cmds_bwamem2:
            subprocess.run(cmd, shell=True)
        
        email_to = Variable.get("email_to", default_var="skigodata@mskcc.org")
        msg_body = " \n ".join(cmds_dragen)
        msg_subject = "DRAGEN commands launched for " + sequencer_and_run
        # use 'echo -e' to try to preserve newline characters
        mail_cmd = "echo -e {} | mail -s {} {}".format(msg_body, msg_subject, email_to)
        subprocess.run(mail_cmd, shell=True)

    def build_bwamem2_cmds(sample_sheet, sequencer_and_run):
        # For all ped-peg create the BWA-MEM2 GRCh37 .bam (~30x slower than the DRAGEN GRCh38 .bam)
        # python3 /igo/work/nabors/tools/wgs_python/wgs_stats_bwa_mem2.py 
        # --project-dir /igo/staging/FASTQ/MICHELLE_0457_AHGFTGDSX2_WGS/Project_08822_NZ/ 
        # --output-dir /igo/staging/stats/naborsd_workspace/PPG/MICHELLE_0457/08822_NZ
        cmd_list = []
        for project in sample_sheet.project_set:
            if "08822" in project:
                project_dir = "/igo/staging/FASTQ/" + sequencer_and_run + "_WGS/" + project
                if not os.path.exists(project_dir):
                    project_dir = "/igo/staging/FASTQ/" + sequencer_and_run + "/" + project
                output_dir = "/igo/staging/stats/" + sequencer_and_run + "/" + project
                cmd = "python3 /igo/work/nabors/tools/wgs_python/wgs_stats_bwa_mem2.py --project-dir {} --output-dir {}".format(project_dir, output_dir)
                print(cmd)
                cmd_list.append(cmd)
        return cmd_list

    def build_dragen_cmds(sample_sheet, sequencer_and_run):
        print("Creating DRAGEN pipeline command for each sample on " + sequencer_and_run)
        # dictionary of Sample_ID->Project
        sample_dict = pandas.Series(sample_sheet.df_ss_data['Sample_Project'].values,index=sample_sheet.df_ss_data['Sample_ID']).to_dict()
        # Create DRAGEN pipeline command, for example:
        # bsub -J RAD_Pt_20_T_IGO_04540_P_15 -o RAD_Pt_20_T_IGO_04540_P_15.out -q dragen -n 48 -M 4 
        # /opt/edico/bin/dragen --ref-dir /staging/ref/GRCh38_graph --enable-duplicate-marking true --enable-map-align-output true --fastq-list /igo/work/luc/DIANA_0441_fastq_list.csv 
        # --output-directory /igo/staging/stats/DIANA_0441_AH2V3TDSX3 --fastq-list-sample-id RAD_Pt_20_T_IGO_04540_P_15 --output-file-prefix DIANA_0441_AH2V3TDSX3___P04540_P__RAD_Pt_20_T_IGO_04540_P_15
        cmd_list = []
        for sample, project in sample_dict.items():
            #for example: DIANA_0441_AH2V3TDSX3___P04540_P__RAD_Pt_20_T_IGO_04540_P_15
            output_prefix = "{}___P{}___{}".format(sequencer_and_run, project.replace("Project_",""), sample)

            bsub = "bsub -J {} -eo /igo/staging/stats/{}/{}.out -q dragen -n 48 -M 4 ".format(sample, sequencer_and_run, sample)
            dragen_cmd_1 = "/opt/edico/bin/dragen --ref-dir /staging/ref/GRCh38_graph --enable-duplicate-marking true --enable-map-align-output true "
            dragen_cmd_2 = "--fastq-list /igo/staging/FASTQ/{}/Reports/fastq_list.csv --output-directory /igo/staging/stats/{} ".format(sequencer_and_run, sequencer_and_run)
            dragen_cmd_3 = "--fastq-list-sample-id {} --output-file-prefix {}".format(sample, output_prefix)
            cmd = bsub + dragen_cmd_1 + dragen_cmd_2 + dragen_cmd_3
            print(cmd)
            cmd_list.append(cmd)
        return cmd_list

    """
    Process dictionary of sample sheet project,recipe and launch stats for each project on the run.
    """
    def launch_stats_via_bash_script(sample_sheet, sequencer_and_run):
        # Make sure output directory exists or DRAGEN commands will fail
        working_dir = "/igo/staging/stats/" + sequencer_and_run
        if not os.path.exists(working_dir):
            os.mkdir(working_dir)
        os.chdir(working_dir)

        # TODO email demux complete starting stats for non "REFERENCE" demuxes

        # sh /home/igo/Scripts/Automate-Stats/Launch_Stats_GRCh37.sh 191029_JOHNSAWYERS_0218_000000000-G4BYY 2>&1 >> /home/igo/log/statsTest.log
        cmd = "sh /home/igo/Scripts/Automate-Stats/Launch_Stats_Airflow.sh {} 2>&1 >> /igo/staging/stats.log".format(sequencer_and_run)
        print(cmd)
        subprocess.run(cmd, shell=True)
        #TODO Launch DetectStatsCompletion.sh
        # sh $SCRIPTPATH/../Automate-Stats/DetectStatsCompletion.sh $RUNNAME &

def test_build_dragen_cmds():
    sample_sheet = SampleSheet("test/DIANA_0441_WGS.csv")
    cmd_list = build_dragen_cmds(sample_sheet, "DIANA_0441_AH2V3TDSX3")
    assert(cmd_list[0]== "bsub -J PS4268T_IGO_04540_Q_10 -eo /igo/staging/stats/DIANA_0441_AH2V3TDSX3/PS4268T_IGO_04540_Q_10.out -q dragen -n 48 -M 4 /opt/edico/bin/dragen --ref-dir /staging/ref/GRCh38_graph --enable-duplicate-marking true --enable-map-align-output true --fastq-list /igo/staging/FASTQ/DIANA_0441_AH2V3TDSX3/Reports/fastq_list.csv --output-directory /igo/staging/stats/DIANA_0441_AH2V3TDSX3 --fastq-list-sample-id PS4268T_IGO_04540_Q_10 --output-file-prefix DIANA_0441_AH2V3TDSX3___P04540_Q___PS4268T_IGO_04540_Q_10___GRCh38")
    print(*cmd_list)
