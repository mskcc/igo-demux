#!/usr/bin/env python3

import os
import re
from subprocess import call
import sys
import csv
import pickle
from dataclasses import dataclass
from collections import OrderedDict
import scripts.generate_run_params
import time
import shutil

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
DO_NOT_PROCESS = ["HumanWholeGenome", "10X_Genomics", "DLP", "MissionBio"]
# this list contains the headers of the columns.  we will access the data using these listings
data_headers = list()

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
		fastq_dir = "/igo/staging/FASTQ/" + run + "/" + row[data_headers.index("Sample_Project")] + "/Sample_" + row[data_headers.index("Sample_ID")] + "/"
		fastqs  = os.listdir(fastq_dir)
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
		# create output directoories
		parent_dir = "/igo/staging/stats/"
		work_dir = parent_dir + run + "/"
		rna_dir = work_dir + "RNA/"
		mwgs_dir = work_dir + "mWGS/"
		
		# check for the work directory
		if not os.path.isdir(work_dir):
			# create the directory if it is not already there
			os.mkdir(work_dir)
			# shutil.rmtree(work_dir)
		os.chdir(work_dir)
		
		#
		# define a list for all RNA samples.  cannot use the metric file name for the RNA bams
		rna_samples_and_gtags = list()
		
		for sample in all_samples:
			# grab the sample parameters (bait set, type, gtag, etc)
			sample_params = self.get_params(sample.genome, sample.recipe)
			# process the RNA data seperately
			if (sample_params["TYPE"] == "RNA"):
				if not os.path.isdir(rna_dir):
					os.mkdir(rna_dir)
				self.rna_alignment_and_metrics(sample, run, sample_params, rna_dir)
				rna_samples_and_gtags.append([sample.sample_id, sample_params["GTAG"]])
				continue
			# process MouseWholeGenome separately
			if (sample.recipe == "MouseWholeGenome"):
				if not os.path.isdir(mwgs_dir):
					os.mkdir(mwgs_dir)
				self.mwgs(sample, run, sample_params)
				continue
			# do the bam alignment
			bams_by_lane, BWAJobNameHeader = self.alignment_to_genome(sample, run, sample_params, work_dir)
			# launch the Picard tools
			self.launch_picard(bams_by_lane, run, sample, sample_params, BWAJobNameHeader)
			# launch rename of RNA metric files
		self.post_data_files(run, rna_samples_and_gtags, work_dir, rna_dir, mwgs_dir)
			
			
	# grab the parameters needed for bwa-mem and picard
	@staticmethod
	def get_params(genome, recipe):
		#
		# call outside scripts and return the parameter data
		sample_params = scripts.generate_run_params.main(["--recipe", recipe, "--species", genome])
		return(sample_params)
	
	
	# let's align the fastqs to the genome!	
	@staticmethod
	def alignment_to_genome(sample, run, sample_params, work_dir):
		#
		# BIG_NODES = "-m \"is01 is02 is03 is04 is05 is06 is07 is08\" -n 60 -M 8 "
		BWAJobNameHeader = run + "___BWA_MEM___"
		os.chdir(work_dir)
		bams_by_lane = list()
		for fq_pair in sample.all_fastqs.lanes:
			fastq_dir = "/igo/staging/FASTQ/" + run + "/" + sample.project + "/Sample_" + sample.sample_id + "/"
			fastq_by_lane = fq_pair.r1[:-16]
			bam_by_lane = fastq_by_lane + ".bam"
			bams_by_lane.append(bam_by_lane)
			fq_for_bwa = [fq_pair.r1, fq_pair.r2]
			if fq_for_bwa[1] is None:
				fq_for_bwa.pop(1)
			bwa_mem = "\"/igoadmin/opt/common/CentOS_7/bwa/bwa-0.7.17/bwa mem -M -t 40 " + sample_params["REFERENCE"] + " " + " ".join(fastq_dir + fq for fq in fq_for_bwa) + " | /igoadmin/opt/common/CentOS_7/samtools/samtools-1.9/bin/samtools view -bS - > " + bam_by_lane + "\""
			bsub_bwa_mem = "bsub -J "+  BWAJobNameHeader + fastq_by_lane + " -o " + BWAJobNameHeader + fastq_by_lane + ".out -n 40 -M 8 " + bwa_mem
			print(bsub_bwa_mem)
			call(bsub_bwa_mem, shell = True)
		return(bams_by_lane, BWAJobNameHeader)
	
	
	# processing the RNA data: Using DRAGEN for the alignment and CollectRNASeqMetrics Picard tool
	@staticmethod
	def rna_alignment_and_metrics(sample, run, sample_params, rna_dir):
		# 
		os.chdir(rna_dir)
		
		PICARD_RNA = "java -Dpicard.useLegacyParser=false -jar /igo/home/igo/resources/picard2.23.2/picard.jar CollectRnaSeqMetrics "
		# prjct = sample.project.split("_")[1]
		prjct = sample.project[8:]
		# gtag = sample_params["GTAG"]
		# if (gtag == "GRCh38"):
		# 	gtag = "hg38"
		# else:
		# 	gtag = "GRCm38"
		RNADragenJobNameHeader = run + "___RNA_DRAGEN___"
		metric_file = run + "___P" + prjct + "___" + sample.sample_id + "___" + sample_params["GTAG"]
		fastq_list = "/igo/staging/FASTQ/" + run + "/Reports/fastq_list.csv "
		launch_dragen_rna = "/opt/edico/bin/dragen -f -r /staging/ref/RNA/" + sample_params["GTAG"]  +  " --fastq-list " + fastq_list + " --fastq-list-sample-id " + sample.sample_id   + " -a " + sample_params["GTF"] + " --enable-map-align true --enable-sort=true --enable-bam-indexing true --enable-map-align-output true --output-format=BAM --enable-rna=true --enable-duplicate-marking true --enable-rna-quantification true " + " --output-file-prefix " + sample.sample_id + " --output-directory ./" 
		bsub_launch_dragen_rna = "bsub -J " + RNADragenJobNameHeader + sample.sample_id + " -o " + RNADragenJobNameHeader + sample.sample_id + ".out -m id01 -q dragen -n 48 -M 4 " + launch_dragen_rna
		print(bsub_launch_dragen_rna)
		call(bsub_launch_dragen_rna, shell = True)
		
		# run Picard RNA metrics tools
		RNAMetricsJobNameHeader = run + "___RNA_METRICS___"
		rnaseq = PICARD_RNA + "--RIBOSOMAL_INTERVALS " + sample_params["RIBOSOMAL_INTERVALS"] + " --STRAND_SPECIFICITY NONE --REF_FLAT " + sample_params["REF_FLAT"] + "  --INPUT " + sample.sample_id + ".bam  " + "--OUTPUT " + metric_file + "___RNA.txt"
		bsub_rnaseq = "bsub -J " + RNAMetricsJobNameHeader + sample.sample_id + " -o " + RNAMetricsJobNameHeader + sample.sample_id + ".out -w \"done(" + RNADragenJobNameHeader + sample.sample_id + ")\" -n 8 -M 8 " + rnaseq
		print(bsub_rnaseq)
		call(bsub_rnaseq, shell = True)
		
		
	@staticmethod
	def mwgs(sample, run, sample_params, mwgs_dir):
		#
		os.chdir(mwgs_dir)
		
		MWGSDragenJobNameHeader = run + "___DRAGEN_mWGS___"
		
		# create metrics file name
		prjct = sample.project[8:]
		metric_file = run + "___P" + prjct + "___" + sample.sample_id + "___" + sample_params["GTAG"]
		fastq_list = "/igo/staging/FASTQ/" + run + "/Reports/fastq_list.csv "
		launch_dragen_mwgs= "/opt/edico/bin/dragen --ref-dir /staging/ref/" + sample_params["GTAG"]  +  " --fastq-list " + fastq_list + " --fastq-list-sample-id " + sample.sample_id + " --intermediate-results-dir /staging/temp --output-directory ./" + " --output-file-prefix " + metric_file + ' --enable-duplicate-marking true'
		bsub_launch_dragen_mwgs = "bsub -J " +  MWGSDragenJobNameHeader  + sample.sample_id + " -o " + MWGSDragenJobNameHeader + sample.sample_id + '.out -m "id01" -q dragen -n 48 -M 4 ' + launch_dragen_mwgs
		print(bsub_launch_dragen_mwgs)
		call(bsub_launch_dragen_mwgs, shell = True)
		
		
	# launch the picrd tools to process the bams
	@staticmethod
	def launch_picard(bams_by_lane, run, sample, sample_params, BWAJobNameHeader):
		#
		# BIG_NODES = " -m \"is01 is02 is03 is04 is05 is06 is07 is08\" -n 60 -M 8 "
		PICARD = "java -Dpicard.useLegacyParser=false -jar /igo/home/igo/resources/picard2.23.2/picard.jar "
		
		# prjct = sample.project.split("_")[1]
		prjct = sample.project[8:]
		metric_file = run + "___P" + prjct + "___" + sample.sample_id + "___" + sample_params["GTAG"]
	
		# merge bams
		MergeBamsJobNameHeader = run + "___MERGE_BAMS___"
		merge_bams = PICARD + "MergeSamFiles --OUTPUT " + sample.sample_id + ".merged.bam " + " ".join("--INPUT " + i for i in bams_by_lane)
		bsub_merge =  "bsub -w \"ended(" + BWAJobNameHeader + sample.sample_id + "*)\" -J " + MergeBamsJobNameHeader + sample.sample_id + " -o " +  MergeBamsJobNameHeader + sample.sample_id + ".out -n 40 -M 8 "
		bsub_merge_bams = bsub_merge + merge_bams
		print(bsub_merge_bams)
		call(bsub_merge_bams, shell = True)
		
		# add or replace read groups
		AORRGJobNameHeader = run + "___AORRG___"
		add_or_replace = PICARD + "AddOrReplaceReadGroups --SORT_ORDER coordinate --CREATE_INDEX true --INPUT " + sample.sample_id + ".merged.bam  " + "--OUTPUT " + sample.sample_id + ".bam  " + "--RGID " + sample.sample_id + "  --RGLB " + sample.sample_id + " --RGPL illumina --RGPU " + sample.project + " --RGSM " + sample.sample_id + " --RGCN GCL@MSKCC"
		bsub_add_or_replace = "bsub -J " + AORRGJobNameHeader + sample.sample_id + " -o " + AORRGJobNameHeader + sample.sample_id + ".out -w \"done(" + MergeBamsJobNameHeader + sample.sample_id + ")\" -n 40 -M 8 " + add_or_replace
		print(bsub_add_or_replace)
		call(bsub_add_or_replace, shell = True)
		
		# mark duplicates
		MarkDupsJobNameHeader = run + "___MARK_DUPLICATES___"
		mark_dup = PICARD + "MarkDuplicates --CREATE_INDEX true --METRICS_FILE " + metric_file + "___MD.txt  " + "--OUTPUT " + sample.sample_id + "___MD.bam  " + "--INPUT " + sample.sample_id + ".bam"
		bsub_mark_dup = "bsub -J " + MarkDupsJobNameHeader + sample.sample_id + " -o " + MarkDupsJobNameHeader + sample.sample_id + ".out -w \"done(" + AORRGJobNameHeader + sample.sample_id + ")\" -n 40 -M 8 " + mark_dup
		print(bsub_mark_dup)
		call(bsub_mark_dup, shell = True)
		
		# alignment summary
		AlignmentJobNameHeader = run + "___ALIGNMENT_SUMMARY___"
		alignment = PICARD + "CollectAlignmentSummaryMetrics --REFERENCE_SEQUENCE " + sample_params["REFERENCE"] + " --INPUT " + sample.sample_id + "___MD.bam  " + "--OUTPUT " + metric_file + "___AM.txt"
		bsub_alignment = "bsub -J " + AlignmentJobNameHeader + sample.sample_id + " -o " + AlignmentJobNameHeader + sample.sample_id + ".out -w \"done(" + MarkDupsJobNameHeader + sample.sample_id + ")\" -n 8 -M 8 " + alignment
		print(bsub_alignment)
		call(bsub_alignment, shell = True)
		
		# determining if we need CollectHsMetrics Picard tool
		if ("BAITS" in sample_params.keys()):
			HsMetricsJobNameHeader = run + "___HS_METRICS___"
			hs_metrics = PICARD + "CollectHsMetrics --INPUT " + sample.sample_id + "___MD.bam  " + " --OUTPUT " + metric_file + "___HS.txt" +  " --REFERENCE_SEQUENCE " + sample_params["REFERENCE"] + " --BAIT_INTERVALS " + sample_params["BAITS"] + " --TARGET_INTERVALS " + sample_params["TARGETS"]
			bsub_hs_metrics = "bsub -J " + HsMetricsJobNameHeader + sample.sample_id + " -o " + HsMetricsJobNameHeader + sample.sample_id + ".out -w \"done(" + MarkDupsJobNameHeader + sample.sample_id + ")\" -n 8 -M 8 " + hs_metrics
			print(bsub_hs_metrics)
			call(bsub_hs_metrics, shell = True)
			
		# let's determine if we need WGS stats
		if (sample_params["TYPE"] == "WGS"):
			WGSMetricsJobNameHeader = run + "___WGS_METRICS___"
			bsub_wait_wgs = "bsub -w \"done(" + MarkDupsJobNameHeader + sample.sample_id + ")\" -J " + WGSMetricsJobNameHeader + sample.sample_id + " -o " + WGSMetricsJobNameHeader + sample.sample_id + ".out -n 8 -M 8 "
			collect_wgs = PICARD + "CollectWgsMetrics --INPUT " + sample.sample_id + "___MD.bam " + "--OUTPUT " + metric_file + "___WGS.txt --REFERENCE_SEQUENCE " + sample_params["REFERENCE"]
			bsub_collect_wgs = bsub_wait_wgs + collect_wgs
			print(bsub_collect_wgs)
			call(bsub_collect_wgs, shell = True)
			
		
	# let's gather the txt data files and move them
	@staticmethod
	def post_data_files(run, rna_samples_and_gtags, work_dir, rna_dir, mwgs_dir):
		#
		sequencer = run.split("_")[0]
		os.chdir(work_dir)
		
		# check to see if this run had any rna samples
		# if rna_samples_and_gtags:
		if os.path.isdir(rna_dir):
			os.chdir(rna_dir)
			rna_samples_and_gtags_file = open('rna_samples_and_gtags.pickle', 'ab')
			pickle.dump(rna_samples_and_gtags, rna_samples_and_gtags_file)
			#
			# create the MD, AM and WGS data files, put them back into the directory 
			RNAMetricsJobName = run + "___RNA_METRICS___*"
			csv_2_txt = "/igo/work/igo/igo-demux/scripts/dragen_csv_2_txt.py ./ ./"
			bsub_csv_2_txt = "bsub -J RNA_CSV_TO_TXT___" + run + " -o RNA_CSV_TO_TXT___" + run + ".out -w \"ended(" + RNAMetricsJobName + ")\" -n 2 -M 8 " + csv_2_txt
			print(bsub_csv_2_txt)
			call(bsub_csv_2_txt, shell = True)
			rename_txt_files = "/igo/work/igo/igo-demux/scripts/rename_txt_files.py " + rna_dir
			bsub_rename_txt_files = "bsub -J RENAME_RNA_TXT_FILES___" + run + " -o RENAME_RNA_TXT_FILES___" + run + ".out -w \"ended(RNA_CSV_TO_TXT___" + run + ")\" -n 2 -M 8 " + rename_txt_files
			print(bsub_rename_txt_files)
			call(bsub_rename_txt_files, shell = True)
			
		if os.path.isdir(mwgs_dir):
			os.chdir(mwgs_dir)
			#
			# create the MD, AM and WGS data files, put them back into the directory
			MWGSDragenJobName = run + "___DRAGEN_mWGS___*"
			csv_2_txt = "/igo/work/igo/igo-demux/scripts/dragen_csv_2_txt.py ./ " + work_dir
			bsub_csv_2_txt = "bsub -J mWGS_CSV_TO_TXT___" + run + " -o mWGS_CSV_TO_TXT___" + run + ".out -w \"ended(" + MWGSDragenJobName + ")\" -n 2 -M 8 " + csv_2_txt
			print(bsub_csv_2_txt)
			call(bsub_csv_2_txt, shell = True)
		
		mv_txt_files = "python3 /igo/work/igo/igo-demux/scripts/mv_txt_files.py " + work_dir
		bsub_mv_all_txt = "bsub -J MOVE_TXT_FILES___" + run + " -o " + "MOVE_TXT_FILES___" + run + ".out -w \"ended(" + run + "___*)\" -n 2 -M 8 " + mv_txt_files
		print(bsub_mv_all_txt)
		call(bsub_mv_all_txt, shell = True)
		
	
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
	
	
	
# extra code
# work_dir_rna = work_dir + "/RNA/"
# make_work_rna_dir = "mkdir -p " + work_dir_rna
# call(make_work_rna_dir, shell = True)
# work_dir_mWGS = work_dir + "/mWGS/"
# make_work_mWGS_dir = "mkdir -p " + work_dir_mWGS
# call(make_work_mWGS_dir, shell = True)
