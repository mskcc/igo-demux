import os
from subprocess import call
import sys
import csv
import pickle
from dataclasses import dataclass
from collections import OrderedDict
import glob
import shutil
import scripts.generate_run_params

# setting up the data classes for the sample sheet structure for launching the metrics
@dataclass
class Fastqs:
	r1: str
	r2: str

@dataclass
class Lanes:
	lanes: list()

@dataclass
class Sample:
	sample_id: str
	genome: str
	recipe: str
	project: str
	all_fastqs: str

# Global Variable : we do not want to process these experiments in this script
DO_NOT_PROCESS = ["HumanWholeGenome", "10X_Genomics", "DLP"]
# These recipes will be evaluated using DRAGEN because of their larger size of fastqs
RUN_ON_DRAGEN = ["MissionBio", "SingleCellCNV", "CustomCapture", "MouseWholeGenome"]
# this list contains the headers of the columns.  we will access the data using these listings
data_headers = list()
PICARD_VERSION = "2_23_2"
PICARD_JAR = "/igo/home/igo/resources/picard2.23.2/picard.jar "
PICARD_AND_JAR = "java -Dpicard.useLegacyParser=false -jar " + PICARD_JAR

# this class handles obtaining the data from the sample sheet and storing it in a data class to use later
class GetSampleData:
	#
	def __init__(self):
		self.all_sample_ids = list()
		self.all_samples = list()
		self.duplicate_sample = ""
		self.all_fastqs = list()
	
	# let's grab the sample sheet and run info to store into the data classes
	# the sample sheet is in csv file format
	def get_samples(self, sample_sheet, run):
		#
		global DO_NOT_PROCESS, data_headers
		
		csv_sample_data = list()
		testing = list()
		with open(sample_sheet) as csv_file:
			csv_reader = csv.reader(csv_file, delimiter = ",")
			got_data = False
			# once we find the "Lane" header, let's store the header row and the data
			for row in csv_reader:
				if (row[0] == "Lane"):
					got_data = True
					data_headers = row
					print(data_headers)
				elif (row[0] != "Lane") and got_data:
					# do not process hWGS, DLP, 10X samples or MissionBio - they have their own processes
					if any(s in row[data_headers.index("Sample_Well")] for s in DO_NOT_PROCESS):
						continue
					self.all_sample_ids.append(row[data_headers.index("Sample_ID")])
					# check for duplicate samples in the sample sheet
					self.duplicate_sample = self.check_this_sample(row[data_headers.index("Sample_ID")], self.all_sample_ids)
					if not self.duplicate_sample:
						# go to a routine to pair the reads.  and return them
						self.all_lanes = self.get_fastqs(self, row, sample_sheet, run)
						sr = Sample(row[data_headers.index("Sample_ID")], row[data_headers.index("Sample_Plate")], row[data_headers.index("Sample_Well")], row[data_headers.index("Sample_Project")], self.all_lanes)
						self.all_samples.append(sr)
				else:
					continue
		#
		return(self.all_samples)
	
	# checking for duplicate reads in the sample sheet
	@staticmethod
	def check_this_sample(sample_id, all_sample_ids):
		#
		x = 0
		duplicate_sample = False
		x = all_sample_ids.count(sample_id)
		if (x > 1):
			duplicate_sample = True
		return(duplicate_sample)
		
	# let's start the process of obtaining the fastqs	
	@staticmethod
	def get_fastqs(self, row, sample_sheet, run):
		#
		# get run from the sample sheet
		fastq_directory = "/igo/staging/FASTQ/{}/{}/Sample_{}/".format(run, row[data_headers.index("Sample_Project")], row[data_headers.index("Sample_ID")])
		fastqs  = os.listdir(fastq_directory)
		# right here, let's eliminate all fastqs from the list that isn't R1 or R2
		fastqs = [x for x in fastqs if ("_R1_001.fastq" in x) or ("_R2_001.fastq" in x)]
		# check the run type: PE or SE
		run_type = self.determine_run_type(fastqs)
		# put the fastqs in order by lane
		fastqs = self.pair_fastqs_by_lane(fastqs, run_type)
		# put the fastqs in the Fastq and Lane data classes
		fastqs = self.convert_fastqs_to_class(fastqs)
		return(fastqs)
	
	# pairing the fastqs according to lane
	@staticmethod
	def pair_fastqs_by_lane(fastqs, run_type):
		""" PAIR FASTQS OF THE SAMPLE BY LANE """
		# grab the number of lanes used
		if run_type == "se":
			pairs = list()
			for fq in fastqs:
				se_pair = [fq, None]
				pairs.append(se_pair)
		else:
			lane_size = int(len(fastqs) / 2)
			pairs = list()
			for i in range(0, lane_size, 1):
				fastq1 = fastqs[0][:-15]
				fastq_pair = [i for i in fastqs if fastq1 in i]
				fastq_pair.sort()
				fastqs.remove(fastq_pair[0]), fastqs.remove(fastq_pair[1])
				pairs.append(fastq_pair)
			pairs.sort()
		return(pairs)
	
	# get that run type: PE or SE
	@staticmethod
	def determine_run_type(fastqs):
		r2 = [i for i in fastqs if "_R2_001.fastq" in i]
		if len(r2) == 0:
			return("se")
		else:
			return("pe")
		
	# storing the fastqs in to the data class
	@staticmethod
	def convert_fastqs_to_class(fastqs):
		#
		paired_fastqs = list()
		for sample in fastqs:
			paired = Fastqs(sample[0], sample[1])
			paired_fastqs.append(paired)
		all_fastqs = Lanes(paired_fastqs)
		return(all_fastqs)
		
# class to grab the run information from the sample sheet
class GetRun(object):
	
	def __init__(self):
		self.run = ""
		
	def get_run(self, sample_sheet):
		self.run = sample_sheet.split("/")[5][19:-4]
		return(self.run)
				
# class to lauch Picard and RNA stats
class LaunchMetrics(object):
	#
	def __init__(self):
		self.bam = ""
		self.rna_samples = list()
		
	def launch_metrics(self, all_samples, run):
		#
		global RUN_ON_DRAGEN
		# create output directories
		parent_directory = "/igo/staging/stats/"
		work_directory = parent_directory + run + "/"
		rna_directory = work_directory + "RNA/"
		dragen_directory = work_directory + "DRAGEN/"
		
		# check for the work directory
		if not os.path.isdir(work_directory):
			# create the directory if it is not already there
			os.mkdir(work_directory)
			# shutil.rmtree(work_directory)
			
			
		
		for sample in all_samples:
			# grab the sample parameters (bait set, type, gtag, etc)
			sample_parameters = self.get_parameters(sample.genome, sample.recipe)
			# process the RNA data seperately
			if (sample_parameters["TYPE"] == "RNA"):
				if not os.path.isdir(rna_directory):
					os.mkdir(rna_directory)
				self.rna_alignment_and_metrics(sample, run, sample_parameters, rna_directory)
				continue
			# check to see if we need to run the samples on dragen
			if any(s in sample.recipe for s in RUN_ON_DRAGEN):
				if not os.path.isdir(dragen_directory):
					os.mkdir(dragen_directory)
				self.dragen(sample, run, sample_parameters, dragen_directory)
				continue
			# do the bam alignment
			bams_by_lane = self.alignment_to_genome(self, sample, run, sample_parameters, work_directory)
			# launch the Picard tools
			self.launch_picard(bams_by_lane, run, sample, sample_parameters, work_directory)
			# launch rename of RNA metric files
		self.post_data_files(run, work_directory, rna_directory, dragen_directory)
			
			
	# grab the parameters needed for bwa-mem and picard
	@staticmethod
	def get_parameters(genome, recipe):
		#
		# call outside scripts and return the parameter data
		sample_parameters = scripts.generate_run_params.main(["--recipe", recipe, "--species", genome])
		return(sample_parameters)
		

	# let's align the fastqs to the genome!	
	@staticmethod
	def alignment_to_genome(self, sample, run, sample_parameters, work_directory):
		#
		# BIG_NODES = "-m \"is01 is02 is03 is04 is05 is06 is07 is08\" -n 60 -M 8 "
		# aorrg_bams_by_lane = list()
		bwa_job_name_header = run + "___BWA_MEM___"
		os.chdir(work_directory)
		bams_by_lane = list()
		for fq_pair in sample.all_fastqs.lanes:
			fastq_directory = "/igo/staging/FASTQ/{}/{}/Sample_{}/".format(run, sample.project, sample.sample_id)
			fastq_by_lane = fq_pair.r1[:-16]
			bam_by_lane = fastq_by_lane + ".bam"
			bwa_mem_job_name = bwa_job_name_header + fastq_by_lane
			bams_by_lane.append(bam_by_lane)
			fq_for_bwa = [fq_pair.r1, fq_pair.r2]
			if fq_for_bwa[1] is None:
				fq_for_bwa.pop(1)
			bwa_mem = "\"/igoadmin/opt/common/CentOS_7/bwa/bwa-0.7.17/bwa mem -M -t 40 {} {} | /igoadmin/opt/common/CentOS_7/samtools/samtools-1.9/bin/samtools view -bS - > {} \"".format(sample_parameters["REFERENCE"], " ".join(fastq_directory + fq for fq in fq_for_bwa), bam_by_lane)
			bsub_bwa_mem = "bsub -J {0} -o {0}.out -cwd \"{1}\" -n 40 -M 8 ".format(bwa_mem_job_name, work_directory) + bwa_mem
			print(bsub_bwa_mem)
			call(bsub_bwa_mem, shell = True)
			# do add or replace read group after alignment by bwa-mem
			# aorrg_bams_by_lane.append(self.add_or_replace_read_groups(bam_by_lane, sample, bwa_mem_job_name, run))
		return(bams_by_lane)
		
	# processing the RNA data: Using DRAGEN for the alignment and CollectRNASeqMetrics Picard tool
	@staticmethod
	def rna_alignment_and_metrics(sample, run, sample_parameters, rna_directory):
		# 
		os.chdir(rna_directory)
		
		# prjct = sample.project.split("_")[1]
		prjct = sample.project[8:]
		
		# get the correct path for the reference
		gtag = sample_parameters["GTAG"]
		if (gtag == "GRCh38"):
			rna_path = "/staging/ref/hg38_alt_masked_graph_v2+cnv+graph+rna-8-1644018559"
		else:
			rna_path = "/staging/ref/RNA/grcm39"
			
		rna_dragen_job_name_header = run + "___RNA_DRAGEN___"
		metric_file = "{}___P{}___{}___{}".format(run, prjct, sample.sample_id, sample_parameters["GTAG"])
		fastq_list = "/igo/staging/FASTQ/" + run + "/Reports/fastq_list.csv "
		
		
		launch_dragen_rna = "/opt/edico/bin/dragen -f -r {} --fastq-list {} --fastq-list-sample-id {} -a {} --intermediate-results-dir /staging/temp --enable-map-align true --enable-sort=true --enable-bam-indexing true --enable-map-align-output true --output-format=BAM --enable-rna=true --enable-duplicate-marking true --enable-rna-quantification true --output-file-prefix {} --output-directory {}".format(rna_path, fastq_list, sample_parameters["GTF"], metric_file, rna_directory)
		bsub_launch_dragen_rna = "bsub -J {0}{1} -o {0}{1}.out -cwd \"{2}\" -m \"id01 id02 id03\" -q dragen -n 48 -M 4 ".format(rna_dragen_job_name_header, sample.sample_id, rna_directory) + launch_dragen_rna
		print(bsub_launch_dragen_rna)
		call(bsub_launch_dragen_rna, shell = True)
		
		# run Picard RNA metrics tools
		rna_metrics_job_name_header = run + "___RNA_METRICS___"
		rnaseq = " CollectRnaSeqMetrics --RIBOSOMAL_INTERVALS {0} --STRAND_SPECIFICITY NONE --REF_FLAT {1} --INPUT {2}.bam --OUTPUT {2}___{4}___RNA.txt".format(sample_parameters["RIBOSOMAL_INTERVALS"], sample_parameters["REF_FLAT"], metric_file, PICARD_VERSION)
		bsub_rnaseq = "bsub -J {0}{1} -o {0}{1}.out -w \"done({2}{1})\" -cwd \"{3}\" -n 8 -M 8 ".format(rna_metrics_job_name_header, sample.sample_id, rna_dragen_job_name_header, rna_directory) + rnaseq
		print(bsub_rnaseq)
		call(bsub_rnaseq, shell = True)
		
		
	@staticmethod
	def dragen(sample, run, sample_parameters, dragen_directory):
		#
		os.chdir(dragen_directory)
		
		dragen_job_name_header = run + "___DRAGEN___"
		
		# create metrics file name
		prjct = sample.project[8:]
		
		# get the correct path for the reference
		gtag = sample_parameters["GTAG"]
		if (gtag == "GRCh38"):
			rna_path = "/staging/ref/hg38_alt_masked_graph_v2+cnv+graph+rna-8-1644018559"
		else:
			rna_path = "/staging/ref/grcm39"
			
		metric_file = "{}___P{}___{}___{}".format(run, prjct, sample.sample_id, sample_parameters["GTAG"])
		fastq_list = "/igo/staging/FASTQ/" + run + "/Reports/fastq_list.csv "
		launch_dragen = "/opt/edico/bin/dragen --ref-dir {} --fastq-list {} --fastq-list-sample-id {} --intermediate-results-dir /staging/temp --output-directory {} --output-file-prefix {} --enable-duplicate-marking true".format(rna_path, fastq_list, sample.sample_id, dragen_directory, metric_file)
		bsub_launch_dragen = "bsub -J {0}{1} -o {0}{1}.out -cwd \"{2}\" -m \"id01 id02 id03\" -q dragen -n 48 -M 4 ".format(dragen_job_name_header, sample.sample_id, dragen_directory) + launch_dragen
		print(bsub_launch_dragen)
		call(bsub_launch_dragen, shell = True)
		
		
	# launch the picrd tools to process the bams
	@staticmethod
	def launch_picard(bams_by_lane, run, sample, sample_parameters, work_directory):
		#
		
		# prjct = sample.project.split("_")[1]
		prjct = sample.project[8:]
		metric_file = "{}___P{}___{}___{}".format(run, prjct, sample.sample_id, sample_parameters["GTAG"])
	
		# merge bams
		# special header for add_or_replace_read_groups
		# aorrg_job_name_header = run + "___AddOrReplaceReadGroups___"
		bwa_mem_job_header = run + "___BWA_MEM___"
		merge_bams_job_name_header = run + "___MERGE_BAMS___"
		merge_bams = PICARD_AND_JAR + "MergeSamFiles --SORT_ORDER coordinate --CREATE_INDEX true --OUTPUT {}.merged.bam {}".format(sample.sample_id, " ".join("--INPUT " + i for i in bams_by_lane)) 		
		bsub_merge =  "bsub -w \"ended({0}{1}*)\" -J {2}{1} -o {2}{1}.out -cwd \"{3}\" -n 40 -M 8 ".format(bwa_mem_job_header, sample.sample_id, merge_bams_job_name_header, work_directory)
		bsub_merge_bams = bsub_merge + merge_bams
		print(bsub_merge_bams)
		call(bsub_merge_bams, shell = True)
		
		# mark duplicates
		mark_duplicates_job_name_header = run + "___MARK_DUPLICATES___"
		mark_dup = PICARD_AND_JAR + "MarkDuplicates --CREATE_INDEX true --METRICS_FILE {0}___{1}___MD.txt  " + "--OUTPUT {2}___MD.bam  --INPUT {2}.merged.bam".format(metric_file, PICARD_VERSION, sample.sample_id)
		bsub_mark_dup = "bsub -J {0}{2} -o {0}{2}.out -w \"done({1}{2})\" -cwd \"{3}\" -n 40 -M 8 ".format(mark_duplicates_job_name_header, merge_bams_job_name_header, sample.sample_id, work_directory) + mark_dup
		print(bsub_mark_dup)
		call(bsub_mark_dup, shell = True)
		
		# alignment summary
		alignment_job_name_header = run + "___ALIGNMENT_SUMMARY___"
		alignment = PICARD_AND_JAR + "CollectAlignmentSummaryMetrics --REFERENCE_SEQUENCE {} --INPUT {}___MD.bam --OUTPUT {}___{}___AM.txt".format(sample_parameters["REFERENCE"], sample.sample_id, metric_file, PICARD_VERSION) 
		bsub_alignment = "bsub -J {0}{1} -o {0}{1}.out -w \"done({2}{1})\" -cwd \"3\" -n 8 -M 8 ".format(alignment_job_name_header, sample.sample_id, mark_duplicates_job_name_header, work_directory) + alignment
		print(bsub_alignment)
		call(bsub_alignment, shell = True)
		
		# determining if we need CollectHsMetrics Picard tool
		if ("BAITS" in sample_parameters.keys()):
			hs_metrics_job_name_header = run + "___HS_METRICS___"
			hs_metrics = PICARD_AND_JAR + "CollectHsMetrics --INPUT {}___MD.bam --OUTPUT {}___{}___HS.txt --REFERENCE_SEQUENCE {} --BAIT_INTERVALS {} --TARGET_INTERVALS {}".format(sample.sample_id, metric_file, PICARD_VERSION, sample_parameters["REFERENCE"], sample_parameters["BAITS"], sample_parameters["TARGETS"])
			bsub_hs_metrics = "bsub -J {0}{1} -o {0}{1}.out -w \"done({2}{1})\" -cwd \"{3}\" -n 8 -M 8 ".format(hs_metrics_job_name_header, sample.sample_id, mark_duplicates_job_name_header,work_directory) + hs_metrics
			print(bsub_hs_metrics)
			call(bsub_hs_metrics, shell = True)
			
		# let's determine if we need WGS stats
		if (sample_parameters["TYPE"] == "WGS"):
			wgs_metrics_job_name_header = run + "___WGS_METRICS___"
			bsub_wait_wgs = "bsub -w \"done({0}{1})\" -J {2}{1} -o {2}{1}.out -cwd \"{3}\" -n 8 -M 8 ".format(mark_duplicates_job_name_header, sample.sample_id, wgs_metrics_job_name_header, work_directory)
			collect_wgs = PICARD_AND_JAR + "CollectWgsMetrics --INPUT {}___MD.bam --OUTPUT {}___{}___WGS.txt --REFERENCE_SEQUENCE {}".format(sample.sample_id, metric_file, PICARD_VERSION, sample_parameters["REFERENCE"])
			bsub_collect_wgs = bsub_wait_wgs + collect_wgs
			print(bsub_collect_wgs)
			call(bsub_collect_wgs, shell = True)
			
		
	# let's gather the txt data files and move them
	@staticmethod
	def post_data_files(run, work_directory, rna_directory, dragen_directory):
		#
		sequencer = run.split("_")[0]
		done_directory = "/igo/stats/DONE/" + sequencer + "/"
	
		os.chdir(work_directory)
		
		# check to see if this run had any rna samples
		if os.path.isdir(rna_directory):
			# create the MD, AM and WGS data files, put them back into the directory 
			rna_metrics_job_name = run + "___RNA_METRICS___*"
			csv_2_txt = "/igo/work/nabors/tools/venvpy3/bin/python /igo/work/igo/igo-demux/scripts/dragen_parse_for_am_and_md.py {} {}".format(rna_directory, work_directory)
			bsub_csv_2_txt = "bsub -J RNA_CSV_TO_TXT___{1} -o RNA_CSV_TO_TXT___{0}.out -w \"ended({1})\" -n 2 -M 8 ".format(run, rna_metrics_job_name) + csv_2_txt
			print(bsub_csv_2_txt)
			call(bsub_csv_2_txt, shell = True)
			rename_bam_files = "/igo/work/nabors/tools/venvpy3/bin/python /igo/work/igo/igo-demux/scripts/rename_rna_bam_files.py {}".format(rna_directory)
			bsub_rename_bam_files = "bsub -J RENAME_RNA_BAM_FILES___{0} -o RENAME_RNA_BAM_FILES___{0}.out -w \"ended(RNA_CSV_TO_TXT___{0})\" -cwd \"{1}\" -n 2 -M 8 ".format(run, rna_directory) + rename_bam_files
			print(bsub_rename_bam_files)
			call(bsub_rename_bam_files, shell = True)
			
		if os.path.isdir(dragen_directory):
			# create the MD, AM and WGS data files, put them back into the directory
			dragen_job_name = run + "___DRAGEN___*"
			csv_2_txt = "/igo/work/nabors/tools/venvpy3/bin/python /igo/work/igo/igo-demux/scripts/dragen_parse_for_am_and_md.py {} {}".format(dragen_directory, work_directory)
			bsub_csv_2_txt = "bsub -J DRAGEN_CSV_TO_TXT___{0} -o DRAGEN_CSV_TO_TXT___{0}.out -w \"ended({1})\" -n 2 -M 8 ".format(run, dragen_job_name) + csv_2_txt
			print(bsub_csv_2_txt)
			call(bsub_csv_2_txt, shell = True)
			
			
		# move all data file to the work directory
		move_all_data_files = "/igo/work/nabors/tools/venvpy3/bin/python3 /igo/work/igo/igo-demux/scripts/move_all_data_files.py {}".format(done_directory)
		bsub_move_all_data_files = "bsub -J move_all_data_files___{0} -o move_all_data_files___{0}.out -w \"ended({0}___*)\" -cwd \"{1}\" -n 2 -M 8 ".format(run, work_directory) + move_all_data_files
		print(bsub_move_all_data_files)
		call(bsub_move_all_data_files, shell = True)

		# move all data file to the work directory
		push_data_to_ngs_and_lims = "/igo/work/nabors/tools/venvpy3/bin/python3 /igo/work/igo/igo-demux/scripts/push_data_to_ngs_and_lims.py {} {}".format(sequencer, run)
		bsub_push_data_to_ngs_and_lims = "bsub -K -J push_data_to_ngs_and_lims___{0} -o push_data_to_ngs_and_lims___{0}.out -w \"ended(move_all_data_files___{0})\" -cwd \"{1}\" -n 2 -M 8 ".format(run, work_directory) + push_data_to_ngs_and_lims
		print(bsub_push_data_to_ngs_and_lims)
		call(bsub_push_data_to_ngs_and_lims, shell = True)


	
def main(sample_sheet):
	# Initaite objects
	get_data = GetSampleData()
	launch_metrics = LaunchMetrics()
	get_run = GetRun()
	
	# let's process the data from the sample sheet
	run = get_run.get_run(sample_sheet)
	all_samples = get_data.get_samples(sample_sheet, run)
	launch_metrics.launch_metrics(all_samples, run)
			
			
############# MAIN ROUTINE
if __name__ == "__main__":
	
	# grab the sample sheet as an argument
	sample_sheet = sys.argv[1]
	
	main(sample_sheet)
	
	
	
# if rna_or_dragen_directory_present:   # this would be set to true
# rna_or_dragen_data_files.append(list(glob.iglob("*WGS.txt")))
# rna_or_dragen_data_files.append(list(glob.iglob("*MD.txt")))
# rna_or_dragen_data_files.append(list(glob.iglob("*AM.txt")))
# rna_or_dragen_data_files = rna_or_dragen_data_files[0] + rna_or_dragen_data_files[1] + rna_or_dragen_data_files[2]
# for data_file in rna_or_dragen_data_files:
# shutil.move(data_file, done_directory)
	
# push_txt_files = "python3 /igo/work/igo/igo-demux/scripts/push_to_ngs_and_lims.py " 
# bsub_mv_all_txt = "bsub -K -J PUSH_DATA___" + run + " -o " + "PUSH_DATA___" + run + ".out -w \"ended(" + run + "___*)\" -n 2 -M 8 " + push_txt_files
# print(bsub_mv_all_txt)
# call(bsub_mv_all_txt, shell = True)
	
	