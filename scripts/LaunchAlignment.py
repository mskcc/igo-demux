#!/usr/bin/env python3

import os
from subprocess import call
import sys
import csv
import pickle
from dataclasses import dataclass
from collections import OrderedDict
import glob
import shutil
import pathlib
import scripts.generate_run_params


# Global Variable : we do not want to process these experiments in this script
DO_NOT_PROCESS = ["10X_Genomics", "DLP"]
# These recipes will be evaluated using DRAGEN because of their larger size of fastqs
RUN_ON_DRAGEN = ["MissionBio", "SingleCellCNV", "MouseWholeGenome", "ChIPSeq", "AmpliconSeq", "HumanWholeGenome"]
# this list contains the headers of the columns.  we will access the data using these listings
PICARD_VERSION = "2_23_2"
PICARD_JAR = "/igo/home/igo/resources/picard2.23.2/picard.jar "
PICARD_AND_JAR = "java -Dpicard.useLegacyParser=false -jar {} ".format(PICARD_JAR)


# class to launch all Metrics
class LaunchAlignment(object):
	#
	def __init__(self):
		self.bam = ""
		self.rna_samples = list()
		
	def launch_alignment(self, all_samples, run, project_directory):
		#
		# create output directories
		parent_directory = "/igo/staging/stats"
		work_directory = "{}/{}/".format(parent_directory, run)
		rna_directory = "{}RNA/".format(work_directory)
		dragen_directory = "{}DRAGEN/".format(work_directory)
		
		# create work directory	
		pathlib.Path(work_directory).mkdir(parents = True, exist_ok = True)
		
		# get the section of igo storage where the fastqs are located
		igo_storage_location = project_directory.split("/")[2]
		
		# check to see if data is from delivery directory.  if so, we need to change "staging" to "delivery" in fastq_list.csv 
		if igo_storage_location == "delivery":
			fastq_list = "/igo/staging/stats/{}/fastq_list.csv".format(run)
			create_delivery_fastq_list = "sed 's/staging/delivery/g' /igo/delivery/FASTQ/{0}/Reports/fastq_list.csv > /igo/staging/stats/{0}/fastq_list.csv".format(run)
			print(create_delivery_fastq_list)
			call(create_delivery_fastq_list, shell = True)
		else: 
			fastq_list = "/igo/staging/FASTQ/{}/Reports/fastq_list.csv ".format(run)
		
		for sample in all_samples:
			# test to see if there are some samples that this script will not process
			if any(s in sample.recipe for s in DO_NOT_PROCESS):
				continue
			# grab the sample parameters (bait set, type, gtag, etc)
			sample_parameters = self.get_parameters(sample.genome, sample.recipe)
			# process the RNA data seperately
			if (sample_parameters["TYPE"] == "RNA"):
				pathlib.Path(rna_directory).mkdir(parents = True, exist_ok = True)
				self.rna_alignment_and_metrics(sample, run, sample_parameters, rna_directory, work_directory, fastq_list)
				continue
			# check to see if we need to run the samples on dragen
			if any(s in sample.recipe for s in RUN_ON_DRAGEN):
				pathlib.Path(dragen_directory).mkdir(parents = True, exist_ok = True)
				self.dragen(sample, run, sample_parameters, work_directory, dragen_directory, fastq_list)
				continue
			# check for methylated samples
			if ((sample.recipe == "MethylCaptureSeq") or (sample.recipe == "WholeGenomeBisulfiteSequencing")):
				pathlib.Path(dragen_directory).mkdir(parents = True, exist_ok = True)
				self.dragen_methylation(sample, run, sample_parameters, work_directory, dragen_directory, fastq_list)
				continue
			# do the bam alignment
			bams_by_lane = self.alignment_to_genome(self, sample, run, sample_parameters, work_directory, igo_storage_location)
			# launch the Picard tools
			self.launch_picard(bams_by_lane, run, sample, sample_parameters, work_directory)
		
		
	# grab the parameters needed for bwa-mem and picard
	@staticmethod
	def get_parameters(genome, recipe):
		#
		# call outside scripts and return the parameter data
		sample_parameters = scripts.generate_run_params.main(["--recipe", recipe, "--species", genome])
		return(sample_parameters)
	
	
	# let's align the fastqs to the genome!	
	@staticmethod
	def alignment_to_genome(self, sample, run, sample_parameters, work_directory, igo_storage_location):
		#
		# BIG_NODES = "-m \"is01 is02 is03 is04 is05 is06 is07 is08\" -n 60 -M 8 "
		# aorrg_bams_by_lane = list()
		bwa_job_name_header = "{}___BWA_MEM___".format(run)
		os.chdir(work_directory)
		bams_by_lane = list()
		for fq_pair in sample.all_fastqs.lanes:
			fastq_directory = "/igo/{}/FASTQ/{}/{}/Sample_{}/".format(igo_storage_location, run, sample.project, sample.sample_id)
			fastq_by_lane = fq_pair.r1[:-16]
			bam_by_lane = "{}.bam".format(fastq_by_lane)
			bwa_mem_job_name = "{}{}".format(bwa_job_name_header, fastq_by_lane)
			bams_by_lane.append(bam_by_lane)
			fq_for_bwa = [fq_pair.r1, fq_pair.r2]
			if fq_for_bwa[1] is None:
				fq_for_bwa.pop(1)
			bwa_mem = "\"/igoadmin/opt/common/CentOS_7/bwa/bwa-0.7.17/bwa mem -M -t 40 {} {} | /igoadmin/opt/common/CentOS_7/samtools/samtools-1.9/bin/samtools view -bS - > {} \"".format(sample_parameters["REFERENCE"], " ".join(fastq_directory + fq for fq in fq_for_bwa), bam_by_lane)
			bsub_bwa_mem = "bsub -J {0} -o {0}.out -cwd \"{1}\" -n 40 -M 8 {2}".format(bwa_mem_job_name, work_directory, bwa_mem)
			print(bsub_bwa_mem)
			call(bsub_bwa_mem, shell = True)
			# do add or replace read group after alignment by bwa-mem
			# aorrg_bams_by_lane.append(self.add_or_replace_read_groups(bam_by_lane, sample, bwa_mem_job_name, run))
		return(bams_by_lane)
	
	
	# processing the RNA data: Using DRAGEN for the alignment and CollectRNASeqMetrics Picard tool
	@staticmethod
	def rna_alignment_and_metrics(sample, run, sample_parameters, rna_directory, work_directory, fastq_list):
		# 
		os.chdir(rna_directory)
		
		# get the correct path for the reference
		if (sample_parameters["GTAG"] == "GRCh38"):
			rna_path = "/igo/work/igo/dragen_hash_tables/hg38_alt_masked_graph_v2+cnv+graph+rna-8-1644018559"
		else:
			rna_path = "igo/work/igo/dragen_hash_tables/{}".format(sample_parameters["GTAG"])
			
		rna_dragen_job_name_header = "{}___RNA_DRAGEN___".format(run)
		
		
		launch_dragen_rna = "/opt/edico/bin/dragen -f -r {} --fastq-list {} --fastq-list-sample-id {} -a {} --intermediate-results-dir /staging/temp --enable-map-align true --enable-sort true --enable-bam-indexing true --enable-map-align-output true --output-format BAM --enable-rna true --enable-duplicate-marking true --enable-rna-quantification true --output-file-prefix {} --output-directory {} ".format(rna_path, fastq_list, sample.sample_id, sample_parameters["GTF"], sample.sample_id, rna_directory)
		bsub_launch_dragen_rna = "bsub -J {0}{1} -o {0}{1}.out -cwd \"{2}\" -m \"id01 id02 id03\" -q dragen -n 48 -M 4 {3}".format(rna_dragen_job_name_header, sample.sample_id, rna_directory, launch_dragen_rna)
		print(bsub_launch_dragen_rna)
		call(bsub_launch_dragen_rna, shell = True)
		
				
		
	@staticmethod
	def dragen(sample, run, sample_parameters, work_directory, dragen_directory, fastq_list):
		#
		os.chdir(dragen_directory)
		dragen_job_name_header = "{}___DRAGEN___".format(run)
		
		# get the correct path for the reference
		if (sample_parameters["GTAG"] == "GRCh38"):
			dragen_path = "/igo/work/igo/dragen_hash_tables/hg38_alt_masked_graph_v2+cnv+graph+rna-8-1644018559"
		else:
			dragen_path = "igo/work/igo/dragen_hash_tables/{}".format(sample_parameters["GTAG"])
			
		metric_file_prefix = "{}___P{}___{}___{}".format(run, sample.project[8:], sample.sample_id, sample_parameters["GTAG"])
		launch_dragen = "/opt/edico/bin/dragen --ref-dir {} --fastq-list {} --fastq-list-sample-id {} --intermediate-results-dir /staging/temp --output-directory {} --output-file-prefix {} --enable-sort true --enable-duplicate-marking true".format(dragen_path, fastq_list, sample.sample_id, dragen_directory, sample.sample_id)
		bsub_launch_dragen = "bsub -J {0}{1} -o {0}{1}.out -cwd \"{2}\" -m \"id01 id02 id03\" -q dragen -n 48 -M 4 {3}".format(dragen_job_name_header, sample.sample_id, dragen_directory, launch_dragen)
		print(bsub_launch_dragen)
		call(bsub_launch_dragen, shell = True)
		
		
		
	@staticmethod
	def dragen_methylation(sample, run, sample_parameters, work_directory, dragen_directory, fastq_list):
		#
		os.chdir(dragen_directory)
		dragen_methylation_job_name_header = "{}___DRAGEN_METHYLATION___".format(run)
		
		# get the correct path for the reference
		if (sample_parameters["GTAG"] == "GRCh38"):
			dragen_path = "/igo/work/igo/dragen_hash_tables/hg38_methylated"
		else:
			dragen_path = "/igo/work/igo/dragen_hash_tables/grcm39_methylated"
			
		metric_file_prefix = "{}___P{}___{}___{}".format(run, sample.project[8:], sample.sample_id, sample_parameters["GTAG"])
		launch_dragen_methylation = "/opt/edico/bin/dragen --enable-methylation-calling true --methylation-protocol directional --ref-dir {} --fastq-list {} --fastq-list-sample-id {} --intermediate-results-dir /staging/temp --output-directory {} --output-file-prefix {} --enable-sort true --enable-duplicate-marking true".format(dragen_path, fastq_list, sample.sample_id, dragen_directory, sample.sample_id)
		bsub_launch_dragen = "bsub -J {0}{1} -o {0}{1}.out -cwd \"{2}\" -m \"id01 id02 id03\" -q dragen -n 48 -M 4 {3}".format(dragen_methylation_job_name_header, sample.sample_id, dragen_directory, launch_dragen_methylation)
		print(bsub_launch_dragen)
		call(bsub_launch_dragen, shell = True)
		
			
			
	# launch the picrd tools to process the bams
	@staticmethod
	def launch_picard(bams_by_lane, run, sample, sample_parameters, work_directory):
		#
		metric_file_prefix = "{}___P{}___{}___{}".format(run, sample.project[8:], sample.sample_id, sample_parameters["GTAG"])
		
		# merge bams
		# special header for add_or_replace_read_groups
		# aorrg_job_name_header = run + "___AddOrReplaceReadGroups___"
		bwa_mem_job_header = "{}___BWA_MEM___".format(run)
		merge_bams_job_name_header = "{}___MERGE_BAMS___".format(run)
		merge_bams = "{} MergeSamFiles --SORT_ORDER coordinate --CREATE_INDEX true --OUTPUT {}.merged.bam {}".format(PICARD_AND_JAR, sample.sample_id, " ".join("--INPUT " + i for i in bams_by_lane)) 		
		bsub_merge_bams =  "bsub -w \"ended({0}{1}*)\" -J {2}{1} -o {2}{1}.out -cwd \"{3}\" -n 40 -M 8 {4}".format(bwa_mem_job_header, sample.sample_id, merge_bams_job_name_header, work_directory, merge_bams)
		print(bsub_merge_bams)
		call(bsub_merge_bams, shell = True)
		
		# mark duplicates
		mark_duplicates_job_name_header = "{}___MARK_DUPLICATES___".format(run)
		mark_dup = "{0} MarkDuplicates --CREATE_INDEX true --METRICS_FILE {1}___{2}___MD.txt --OUTPUT {3}___MD.bam  --INPUT {3}.merged.bam".format(PICARD_AND_JAR, metric_file_prefix, PICARD_VERSION, sample.sample_id)
		bsub_mark_dup = "bsub -J {0}{2} -o {0}{2}.out -w \"done({1}{2})\" -cwd \"{3}\" -n 40 -M 8 {4}".format(mark_duplicates_job_name_header, merge_bams_job_name_header, sample.sample_id, work_directory, mark_dup)
		print(bsub_mark_dup)
		call(bsub_mark_dup, shell = True)
		
		