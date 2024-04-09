import os
import re
import subprocess
from datetime import datetime, timedelta

import pandas
from SampleSheet import SampleSheet
import scripts.organise_fastq_split_by_lane
import scripts.get_total_reads_from_demux
import scripts.cellranger
import scripts.calculate_stats
import scripts.get_sequencing_read_data
import scripts.upload_stats
import Fingerprinting.fingerprinting_dag

from airflow import DAG
from airflow.operators.python import PythonOperator
from airflow.operators.email_operator import EmailOperator
from airflow.models import Variable
from airflow.utils.email import send_email


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

    'params': {'dragen_demux', 'False',
             'samplesheet': '/igo/work/igo/SampleSheetCopies/SampleSheet_211206_JOHNSAWYERS_0317_000000000-K3LFK.csv',
             'sequencer_path': '/igo/sequencers/johnsawyers/211206_JOHNSAWYERS_0317_000000000-K3LFK'},
    """
    def demux(ds, **kwargs):
        dragen_demux = eval(kwargs["params"]["dragen_demux"])
        sequencer_path = kwargs["params"]["sequencer_path"]
        samplesheet_path = kwargs["params"]["samplesheet"]
        print("Starting demux {} {}".format(sequencer_path, samplesheet_path))

        samplesheet = os.path.basename(samplesheet_path)
        samplesheet_no_ext = os.path.splitext(samplesheet)[0]  # SampleSheet_210331_MICHELLE_0360_BH5KFYDRXY
        sequencer_and_run = samplesheet_no_ext[19:]            # remove 'SampleSheet_210331_'

        sample_sheet = SampleSheet(samplesheet_path)
        output_directory = "/igo/staging/FASTQ/" + sequencer_and_run

        # if output directory alrady exists, delete it before start demux
        if os.path.exists(output_directory):
            remove_cmd = "rm -rf {}".format(output_directory) 
            print(remove_cmd)
            subprocess.run(remove_cmd, shell=True)

        # Let's check to see if this run is an Cellranger ATAC run
        atac, use_bases_mask = scripts.get_sequencing_read_data.main(sequencer_path)
        
        # check if the sample sheet contains DLP project
        is_DLP = False
        if "DLP" in sample_sheet.recipe_set:
            is_DLP = True
            dragen_demux = True
        
        demux_command = ""
        # -K - wait for the job to complete
        if (dragen_demux) or ("AMELIE" in sequencer_path):
            bsub_command = "bsub -K -n48 -q dragen -eo " + output_directory + "/dragen-demux.log "
            # same as bcl-convert arguments except:  "--bcl-conversion-only true --bcl-only-matched-reads true"
            demux_command = bsub_command + "/opt/edico/bin/dragen --bcl-conversion-only true --bcl-only-matched-reads true --force --bcl-sampleproject-subdirectories true --bcl-input-directory \'{}\' --output-directory \'{}\' --sample-sheet \'{}\'".format(sequencer_path, output_directory, samplesheet_path)
        elif (atac):
            # use cellranger atac mkfastq to demux if run is ATAC
            bsub_command = "bsub -K -n72 -M 6 -eo " +  output_directory + "/bcl2fastq--demux.log "
            demux_command =  bsub_command + "/opt/common/CentOS_6/bcl2fastq/bcl2fastq2-v2.20.0.422/bin/bcl2fastq --minimum-trimmed-read-length 0 --mask-short-adapter-reads 0 --ignore-missing-bcl --runfolder-dir \'{}\' --sample-sheet \'{}\' --output-dir \'{}\' --use-bases-mask \'{}\' --create-fastq-for-index-reads --ignore-missing-filter --ignore-missing-positions --ignore-missing-control --barcode-mismatches 1 --processing-threads 72".format(sequencer_path, samplesheet_path, output_directory, use_bases_mask)
            print("Running demux command: " + demux_command)
            subprocess.run(demux_command, shell=True, check=True)
            # cellranger mkfastq/bcl2fastq doesn't need make demux report or fix fastq list csv file
            scripts.organise_fastq_split_by_lane.create_fastq_folders(output_directory)
            return demux_command
        else: # default to bcl-convert
            bsub_command = "bsub -K -n72 -m \"is01 is02 is03 is04 is05 is06 is07 is08\" -eo " + output_directory + "/bcl-convert.log "
            demux_command = bsub_command + "/usr/bin/bcl-convert --force --bcl-sampleproject-subdirectories true --bcl-input-directory \'{}\' --output-directory \'{}\' --sample-sheet \'{}\'".format(sequencer_path, output_directory, samplesheet_path)
        print("Running demux command: " + demux_command)
        subprocess.run(demux_command, shell=True, check=True)

        # if the demux was successful:
        if not is_DLP:
            print("Adding sample sub-folders to the DRAGEN demux.")
            scripts.organise_fastq_split_by_lane.create_fastq_folders(output_directory)
            scripts.organise_fastq_split_by_lane.correct_fastq_list_csv(output_directory+"/Reports")

        # Call CopyIlluminaReports.sh /igo/staging/FASTQ/RUTH_0066_BHTJ33DRXY
        copy_reports_cmd = "/igo/work/igo/igo-demux/scripts/CopyIlluminaReports.sh /igo/staging/FASTQ/" + sequencer_and_run
        print("Running command to copy demux reports: " + copy_reports_cmd)
        subprocess.run(copy_reports_cmd, shell=True)
        
        return demux_command

    def get_dlp_chip(samplesheet, project):
        samplesheet.df_ss_data.reset_index()
        for index, row in samplesheet.df_ss_data.iterrows():
            if row['Sample_Well'] == 'DLP' and project == row['Sample_Project']:
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

    def stats(ds, **kwargs):
        import requests
        from pathlib import Path
        sequencer_path = kwargs["params"]["sequencer_path"]
        samplesheet_path = kwargs["params"]["samplesheet"]
        samplesheet = os.path.basename(samplesheet_path)
        samplesheet_no_ext = os.path.splitext(samplesheet)[0]  # SampleSheet_210331_MICHELLE_0360_BH5KFYDRXY
        sequencer_and_run = samplesheet_no_ext[19:]            # remove 'SampleSheet_210331_'

        sample_sheet = SampleSheet(samplesheet_path)
        if "REFERENCE" in samplesheet_path:
            return "No stats for reference "  + samplesheet_path

        if "DLP" in sample_sheet.recipe_set:
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
                output_yaml = fastq_project_dir + chip_number + "_metadata.yaml"
                fld_endpoint = "https://igolims.mskcc.org:8443/LimsRest/getDLPFieldMapFile?chipNumber={}".format(chip_number)
                fld_file_path = fastq_project_dir + chip_number + ".fld"
                fld_file = Path(fld_file_path)
                response = requests.get(fld_endpoint, auth = ("pms", "tiagostarbuckslightbike"), verify = False)
                fld_file.write_bytes(response.content)
                if "FAUCI" in sequencer_and_run:
                    python_cmd = "python scripts/yaml/generate_metadata.py " + fastq_project_dir + " " + sample_sheet_path + " " + stats + " " + run_info + " " + fld_file_path + " " + project + " " + output_yaml
                else:
                    python_cmd = "python scripts/yaml/generate_metadata.py " + fastq_project_dir + " " + sample_sheet_path + " " + stats + " " + run_info + " " + fld_file_path + " " + project + " " + output_yaml + " --revcomp_i5"

                print("Calling DLP generate yaml command: {}".format(python_cmd))
                subprocess.check_output(python_cmd, cwd="/home/igo/shared-single-cell", shell=True)

            return "DLP stats posted and yaml file generated"

        # check if the run is 10X by read length
        atac, use_bases_mask = scripts.get_sequencing_read_data.main(sequencer_path)
        print("read length: {}".format(use_bases_mask))
        if use_bases_mask == [29, 89] or atac:
            # if is atac run, demux is using cellranger mkfastq
            if atac:
                scripts.get_total_reads_from_demux.by_json(sequencer_and_run)
                scripts.upload_stats.upload_stats(sequencer_and_run)

                # launch cell ranger based on recipe
                sequencer_and_run_prefix = "_".join(sequencer_and_run.split("_")[0:3])
                scripts.cellranger.launch_cellranger_by_sample_sheet(sample_sheet, sequencer_and_run_prefix)

            else:
                # step 1, generate txt files containing total reads and upload to qc website
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

        scripts.calculate_stats.main(samplesheet_path)

        # add DONE file when all the stats finished, -K to wait until finish
        cmd = 'bsub -K -J wait_stats_done_for_{} -w \"done(uplaodWGSstats{}*)\" touch /igo/staging/stats/{}/DONE'.format(sequencer_and_run, sequencer_and_run, sequencer_and_run)
        print(cmd)
        subprocess.run(cmd, shell=True)

        return "Completed"

    def fingerprinting(ds, **kwargs):
        # read in sample sheet as arguments, filter out projects that need to run fingerprinting
        recipe_list_for_fp = [".*IMPACT*", ".*Heme*", "IDT_Exome*", "WholeExomeSequencing", "Twist_Exome", "MSK-ACCESS*", "CMO-CH", "HumanWholeGenome"]
        # call fingerprinting_dag.py for each project
        samplesheet_path = kwargs["params"]["samplesheet"]

        if "REFERENCE" in samplesheet_path:
            return "No fingerprinting for reference "  + samplesheet_path

        # get project list for running fingerprinting by recipe
        sample_sheet = SampleSheet(samplesheet_path)
        # dictionary of project_ID->genome
        project_genome_dict = pandas.Series(sample_sheet.df_ss_data['Sample_Plate'].values,index=sample_sheet.df_ss_data['Sample_Project']).to_dict()
        project_list_to_run = []        
        for project, recipe in sample_sheet.project_dict.items():
            # fingerprinting only support human
            if project_genome_dict[project] == "Human":
                for recipe_list_item in recipe_list_for_fp:
                    print(project, recipe)
                    expr = re.compile(recipe_list_item)
                    if expr.match(recipe):
                        project_list_to_run.append(project)
                        break
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
            content = "{} demux done".format(sequencer_and_run)
        else:
            content = "{} stats done".format(sequencer_and_run)

        send_email(
            to=["skigodata@mskcc.org"],
            subject='IGO Cluster Stats Finished',
            html_content=content
        )

    # add retry number to fix sample_well issue when two task launch at the same time
    demux_run = PythonOperator(
        task_id='start_the_demux',
        python_callable=demux,
        provide_context=True,
        retries=2, 
        retry_delay=timedelta(seconds=60),
        email_on_retry=True, 
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

    demux_run >> launch_stats >> launch_fingerprinting >> send_stats_email


    def launch_wgs_stats(sample_sheet, sequencer_and_run):
        # Make sure DRAGEN commands do not fail for non-existent directory
        stats_path = "/igo/staging/stats/" + sequencer_and_run
        if not os.path.exists(stats_path):
            os.makedirs(stats_path)
            print("Created the stats directory: {}".format(stats_path))
        
        cmds_dragen = build_dragen_cmds(sample_sheet, sequencer_and_run)
        for cmd in cmds_dragen:
            subprocess.run(cmd, shell=True)
        
        # Only for 08822* PED-PEG Projects also create bwamem2 .bam
        cmds_bwamem2 = build_bwamem2_cmds(sample_sheet, sequencer_and_run)
        for cmd in cmds_bwamem2:
            subprocess.run(cmd, shell=True)

        # create txt stats files from dragen result after dragen command finish for one run directly to /igo/stats/DONE/<sequncer> folder
        sequencer = sequencer_and_run.split("_")[0]
        stats_path_for_conversion = stats_path + "/"
        stats_done_dir = "/igo/stats/DONE/" + sequencer + "/"
        cmd_conversion = "python /igo/work/igo/igo-demux/scripts/dragen_csv_to_picard.py {} {}".format(stats_path_for_conversion, stats_done_dir)
        bsub_command_conversion = "bsub -J create_txt_{} -o {}create_txt.out -w \"done({}*)\" {}".format(sequencer_and_run, stats_path_for_conversion, sequencer_and_run, cmd_conversion)
        print(bsub_command_conversion)
        subprocess.run(bsub_command_conversion, shell=True)

        # call endpoint to push data to ngs database and LIMS
        upload_stats = "python /igo/work/igo/igo-demux/scripts/upload_stats.py {}".format(sequencer_and_run)
        upload_stats_cmd = 'bsub -J uplaodWGSstats{} -o {}uplaodWGSstats.out -w \"done(create_txt_{}*)\" \"{}\"'.format(sequencer_and_run, stats_path_for_conversion, sequencer_and_run, upload_stats)
        print(upload_stats_cmd)
        subprocess.run(upload_stats_cmd, shell=True)


    def build_bwamem2_cmds(sample_sheet, sequencer_and_run):
        # For all ped-peg create the BWA-MEM2 GRCh37 .bam (~30x slower than the DRAGEN GRCh38 .bam)
        # python3 /igo/work/nabors/tools/wgs_python/wgs_stats_bwa_mem2.py 
        # --project-dir /igo/staging/FASTQ/MICHELLE_0457_AHGFTGDSX2_WGS/Project_08822_NZ/ 
        # --output-dir /igo/staging/stats/naborsd_workspace/PPG/MICHELLE_0457/08822_NZ
        cmd_list = []
        for project in sample_sheet.project_set:
            if "08822" in project:
                project_dir = "/igo/staging/FASTQ/" + sequencer_and_run + "_PPG/" + project # note special "_PPG" fastq directory
                output_dir = "/igo/staging/stats/" + sequencer_and_run + "/" + project
                cmd = "python3 /igo/work/nabors/tools/wgs_python/bwa_mem2_only.py --project-dir {} --output-dir {}".format(project_dir, output_dir)
                print(cmd)
                cmd_list.append(cmd)
        return cmd_list

    def build_dragen_cmds(sample_sheet, sequencer_and_run):
        print("Creating DRAGEN pipeline command for each sample on " + sequencer_and_run)
        # dictionary of Sample_ID->Project
        sample_dict = pandas.Series(sample_sheet.df_ss_data['Sample_Project'].values,index=sample_sheet.df_ss_data['Sample_ID']).to_dict()
        # Create DRAGEN pipeline command, for example:
        # bsub -J RAD_Pt_20_T_IGO_04540_P_15 -o RAD_Pt_20_T_IGO_04540_P_15.out -q dragen -n 48 -M 4 
        # /opt/edico/bin/dragen --ref-dir /staging/ref/hg38_alt_masked_graph_v2+cnv+graph+rna-8-1644018559 --enable-duplicate-marking true --enable-map-align-output true --fastq-list /igo/work/luc/DIANA_0441_fastq_list.csv 
        # --output-directory /igo/staging/stats/DIANA_0441_AH2V3TDSX3 --fastq-list-sample-id RAD_Pt_20_T_IGO_04540_P_15 --output-file-prefix DIANA_0441_AH2V3TDSX3___P04540_P__RAD_Pt_20_T_IGO_04540_P_15
        cmd_set = set()
        
        # get prefix from the sequencer_and_run with keeping only machineName_runID_flowcellID
        sequencer_and_run_prefix = "_".join(sequencer_and_run.split("_")[0:3])

        for sample, project in sample_dict.items():
            print("PROJECT: {} {}".format(project, sample_sheet.project_dict[project]))
            if sample_sheet.sample_dict[sample] == "HumanWholeGenome":
                #for example: DIANA_0441_AH2V3TDSX3___P04540_P__RAD_Pt_20_T_IGO_04540_P_15
                output_prefix = "{}___P{}___{}".format(sequencer_and_run_prefix, project.replace("Project_",""), sample)
                job_name = sequencer_and_run + "_" + sample
                bsub = "bsub -J {} -eo /igo/staging/stats/{}/{}.out -q dragen -m \"id01 id02 id03\" -n 48 -M 4 ".format(job_name, sequencer_and_run, sample)
                dragen_cmd_1 = "/opt/edico/bin/dragen --ref-dir /staging/ref/hg38_alt_masked_graph_v2+cnv+graph+rna-8-1644018559 --intermediate-results-dir /staging/temp --enable-duplicate-marking true --enable-map-align-output true "
                dragen_cmd_2 = "--fastq-list /igo/staging/FASTQ/{}/Reports/fastq_list.csv --output-directory /igo/staging/stats/{} ".format(sequencer_and_run, sequencer_and_run)
                dragen_cmd_3 = "--fastq-list-sample-id {} --output-file-prefix {}".format(sample, output_prefix)
                cmd = bsub + dragen_cmd_1 + dragen_cmd_2 + dragen_cmd_3
                cmd_set.add(cmd)
        return cmd_set 
