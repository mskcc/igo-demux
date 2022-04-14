#!/usr/bin/env python3

import os
import re
from subprocess import call
import sys
import csv
from dataclasses import dataclass
from collections import OrderedDict
import generate_run_params

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
	sample_dir: str
	sample_id: str
	genome: str
	recipe: str
	project: str
	all_fastqs: str

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
		csv_sample_data = list()
		testing = list()
		with open(sample_sheet) as csv_file:
			csv_reader = csv.reader(csv_file, delimiter = ",")
			got_data = False
			# once we find the "Lane" header, let's begin to store the data
			for row in csv_reader:
				if (row[0] == "Lane"):
					got_data = True
				elif (row[0] != "Lane") and got_data:
					# do not process hWGS, mWGS, DLP or 10X samples.  they have their own processes
					if ((row[4] == "HumanWholeGenome") or (row[4] == "MouseWholeGenome") or (row[4] == "DLP") or ("10X_Genomics" in row[4])):
						continue
					self.all_sample_ids.append(row[2])
					# check for duplicate samples in the sample sheet
					self.duplicate_sample = self.check_this_sample(row[2], self.all_sample_ids)
					if not self.duplicate_sample:
						# go to a routine to pair the reads.  and return them
						self.all_lanes = self.get_fastqs(self, row, sample_sheet, run)
						sr = Sample(row[1], row[2], row[3], row[4], row[8], self.all_lanes)
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
		fastq_dir = "/igo/staging/FASTQ/" + run + "/" + row[8] + "/" + row[1] + "/"
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
		
	def launch_metrics(self, all_samples, run):
		#
		# create output directoories
		work_dir = "/igo/staging/stats/" + run
		work_dir_rna = "/igo/staging/stats/" + run + "/RNA/"
		make_work_dir = "mkdir -p " + work_dir
		make_work_rna_dir = "mkdir -p " + work_dir_rna
		print(work_dir)
		print(work_dir_rna)
		call(make_work_dir, shell = True)
		call(make_work_rna_dir, shell = True)
		#
		for sample in all_samples:
			# grab the sample parameters (bait set, type, gtag, etc)
			sample_params = self.get_params(sample.genome, sample.recipe)
			# process the RNA data seperately
			if (sample_params["TYPE"] == "RNA"):
				self.rna_alignment_and_metrics(sample, run, sample_params, work_dir_rna)
				continue
			# do the bam alignment
			bams_by_lane = self.alignment_to_genome(sample, run, sample_params, work_dir)
			# launch the Picard tools
			picard_data = self.launch_picard(bams_by_lane, run, sample, sample_params)
			
	# grab the parameters needed for bwa-mem and picard
	@staticmethod
	def get_params(genome, recipe):
		#
		parameter_placement = list
		recipe_and_genome = ["--recipe", recipe, "--species", genome]
		# calll outside scripts and return the parameter data
		sample_params = generate_run_params.main(recipe_and_genome)
		return(sample_params)
	
	# let's align the fastqs to the genome!	
	@staticmethod
	def alignment_to_genome(sample, run, sample_params, work_dir):
		#
		BIG_NODES = "-m \"is01 is02 is03 is04 is05 is06 is07 is08\" -n 60 -M 8 "
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
			bwa_mem = "\"/igoadmin/opt/common/CentOS_7/bwa/bwa-0.7.17/bwa mem -M -t 60 " + sample_params["REFERENCE"] + " " + " ".join(fastq_dir + fq for fq in fq_for_bwa) + " | /igoadmin/opt/common/CentOS_7/samtools/samtools-1.9/bin/samtools view -bS - > " + bam_by_lane + "\""
			bsub_bwa_mem = "bsub -J bwa_mem___" + fastq_by_lane + " -o " + "bwa_mem___" + fastq_by_lane + ".out " + BIG_NODES + bwa_mem
			print(bsub_bwa_mem)
			call(bsub_bwa_mem, shell = True)
		return(bams_by_lane)	
	
	# processing the RNA data: Using DRAGEN for the alignment and CollectRNASeqMetrics Picard tool
	@staticmethod
	def rna_alignment_and_metrics(sample, run, sample_params, work_dir_rna):
		# 
		os.chdir(work_dir_rna)
		PICARD_RNA = "java -Dpicard.useLegacyParser=false -jar /igo/home/igo/resources/picard2.23.2/picard.jar CollectRnaSeqMetrics "
		prjct = sample.project.split("_")[1]
		metric_file = run + "___P" + prjct + "___" + sample.sample_id + "___" + sample_params["GTAG"]
		fastq_list = "/igo/staging/FASTQ/" + run + "/Reports/fastq_list.csv "
		launch_dragen_rna = "/opt/edico/bin/dragen -f -r /staging/ref/RNA/" + sample_params["GTAG"]  +  " --fastq-list " + fastq_list + " --fastq-list-sample-id " + sample.sample_id   + " -a " + sample_params["GTF"] + " --enable-map-align true --enable-sort=true --enable-bam-indexing true --enable-map-align-output true --output-format=BAM --enable-rna=true --enable-duplicate-marking true --enable-rna-quantification true " + " --output-file-prefix " + sample.sample_id + " --output-directory " + work_dir_rna
		bsub_launch_dragen_rna = "bsub -J DRAGEN_RNA___" + sample.sample_id + " -o " + "DRAGEN_RNA___" + sample.sample_id + '.out -m "id01" -q dragen -n 48 -M 4 ' + launch_dragen_rna
		print(bsub_launch_dragen_rna)
		call(bsub_launch_dragen_rna, shell = True)
		if sample_params["TYPE"] == "RNA":
			rnaseq = PICARD_RNA + "--RIBOSOMAL_INTERVALS " + sample_params["RIBOSOMAL_INTERVALS"] + " --STRAND_SPECIFICITY NONE --REF_FLAT " + sample_params["REF_FLAT"] + "  --INPUT " + sample.sample_id + ".bam  " + "--OUTPUT " + metric_file + "___RNA.txt"
			bsub_rnaseq = "bsub -J rnaseq___" + sample.sample_id + " -o " + "rnaseq___" + sample.sample_id + ".out -w \"done(DRAGEN_RNA___" + sample.sample_id + ")\" -n 8 -M 8 " + rnaseq
			print(bsub_rnaseq)
			call(bsub_rnaseq, shell = True)
		
	# launch the picrd tools to process the bams
	@staticmethod
	def launch_picard(bams_by_lane, run, sample, sample_params):
		# 
		BIG_NODES = " -m \"is01 is02 is03 is04 is05 is06 is07 is08\" -n 60 -M 8 "
		PICARD = "java -Dpicard.useLegacyParser=false -jar /igo/home/igo/resources/picard2.23.2/picard.jar "
		#
		prjct = sample.project.split("_")[1]
		metric_file = run + "___P" + prjct + "___" + sample.sample_id + "___" + sample_params["GTAG"]
		#
		# merge bams
		bsub_merge =  "bsub -w \"ended(bwa_mem___" + sample.sample_id + "*)\" -J merge_bam_files___" + sample.sample_id + " -o merge_bam_files___" + sample.sample_id + ".out " + BIG_NODES
		merge_bams = PICARD + "MergeSamFiles --OUTPUT " + sample.sample_id + ".merged.bam " + " ".join("--INPUT " + i for i in bams_by_lane)
		bsub_merge_bams = bsub_merge + merge_bams
		print(bsub_merge_bams)
		call(bsub_merge_bams, shell = True)
		# add or replace read groups
		add_or_replace = PICARD + "AddOrReplaceReadGroups --SORT_ORDER coordinate --CREATE_INDEX true --INPUT " + sample.sample_id + ".merged.bam  " + "--OUTPUT " + sample.sample_id + ".bam  " + "--RGID " + sample.sample_id + "  --RGLB " + sample.sample_id + " --RGPL illumina  --RGPU CUSTOM-BAM  --RGSM " + sample.sample_id + " --RGCN GCL@MSKCC"
		bsub_add_or_replace = "bsub -J add_or_replace_read_groups___" + sample.sample_id + " -o " + "add_or_replace_read_groups___" + sample.sample_id + ".out -w \"done(merge_bam_files___" + sample.sample_id + ")\" " + BIG_NODES + add_or_replace
		print(bsub_add_or_replace)
		call(bsub_add_or_replace, shell = True)
		#
		# mark duplicates
		mark_dup = PICARD + "MarkDuplicates --CREATE_INDEX true --METRICS_FILE " + metric_file + "___MD.txt  " + "--OUTPUT " + sample.sample_id + ".md.bam  " + "--INPUT " + sample.sample_id + ".bam"
		bsub_mark_dup = "bsub -J mark_dup___" + sample.sample_id + " -o " + "mark_dup___" + sample.sample_id + ".out -w \"done(add_or_replace_read_groups___" + sample.sample_id + ")\" " + BIG_NODES + mark_dup
		print(bsub_mark_dup)
		call(bsub_mark_dup, shell = True)
		#
		# alignment summary
		alignment = PICARD + "CollectAlignmentSummaryMetrics --REFERENCE_SEQUENCE " + sample_params["REFERENCE"] + " --INPUT " + sample.sample_id + ".md.bam  " + "--OUTPUT " + metric_file + "___AM.txt"
		bsub_alignment = "bsub -J alignment_summary___" + sample.sample_id + " -o " + "alignment_summary___" + sample.sample_id + ".out -w \"done(mark_dup___" + sample.sample_id + ")\" -n 8 -M 8 " + alignment
		print(bsub_alignment)
		call(bsub_alignment, shell = True)
		# determining if we need CollectHsMetrics Picard tool
		if ("BAITS" in sample_params.keys()):
			hs_metrics = PICARD + "CollectHsMetrics --INPUT " + sample.sample_id + ".md.bam  " + " --OUTPUT " + metric_file + "___HS.txt" +  " --REFERENCE_SEQUENCE " + sample_params["REFERENCE"] + " --BAIT_INTERVALS " + sample_params["BAITS"] + " --TARGET_INTERVALS " + sample_params["TARGETS"]
			bsub_hs_metrics = "bsub -J hs_metrics___" + sample.sample_id + " -o " + "hs_metrics___" + sample.sample_id + ".out -w \"done(mark_dup___" + sample.sample_id + ")\" -n 8 -M 8 " + hs_metrics
			print(bsub_hs_metrics)
			call(bsub_hs_metrics, shell = True)
			
	
def main():
	
	# grab the sample sheet as an argument
	sample_sheet = sys.argv[1]
	
	# Initaite objects
	get_data = GetSampleData()
	launch_metrics = LaunchMetrics()
	get_run = GetRun()
	
	# let's process the data from the sample sheet
	run = get_run.get_run(sample_sheet)
	all_samples = get_data.get_samples(sample_sheet, run)
	all_metrics = launch_metrics.launch_metrics(all_samples, run)
			
			
############# MAIN ROUTINE
if __name__ == "__main__":
	main()
	