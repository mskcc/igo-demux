import os
import re
from subprocess import call
import sys
import csv
from dataclasses import dataclass
from collections import OrderedDict
import scripts.generate_run_params


@dataclass
class Fastqs:
	r1: str
	r2: str

@dataclass
class Lanes:
	lanes: list()
	# lanes: list[Fastqs]


@dataclass
class Sample:
	sample_dir: str
	sample_id: str
	genome: str
	recipe: str
	project: str
	all_fastqs: str

class GetSampleData:
	#
	def __init__(self):
		self.all_sample_ids = list()
		self.all_samples = list()
		self.duplicate_sample = ""
		self.all_fastqs = list()
		
		
	def get_samples(self, sample_sheet, run):
		csv_sample_data = list()
		testing = list()
		with open(sample_sheet) as csv_file:
			csv_reader = csv.reader(csv_file, delimiter = ",")
			got_data = False
			for row in csv_reader:
				if (row[0] == "Lane"):
					got_data = True
				elif (row[0] != "Lane") and got_data:
					# exclude samples with HumanWholeGenome recipe
					if (row[4] != "HumanWholeGenome"):
						self.all_sample_ids.append(row[2])
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
	
	
	@staticmethod
	def check_this_sample(sample_id, all_sample_ids):
		#
		x = 0
		duplicate_sample = False
		x = all_sample_ids.count(sample_id)
		if (x > 1):
			duplicate_sample = True
		return(duplicate_sample)
		
	@staticmethod
	def get_fastqs(self, row, sample_sheet, run):
		#
		# get run from the sample sheet
		fastq_dir = "/igo/staging/FASTQ/" + run + "/" + row[8] + "/Sample_" + row[1] + "/"
		fastqs  = os.listdir(fastq_dir)
		run_type = self.determine_run_type(fastqs)
		#if run_type == "pe":
		fastqs = self.pair_fastqs_by_lane(fastqs, run_type)
		fastqs = self.convert_fastqs_to_class(fastqs)
		return(fastqs)
	
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
			# fastq_pair = tuple(fastq_pair)
				fastqs.remove(fastq_pair[0]), fastqs.remove(fastq_pair[1])
				pairs.append(fastq_pair)
			pairs.sort()
			# print(pairs)
		return(pairs)
	
	@staticmethod
	def determine_run_type(fastqs):
		r2 = [i for i in fastqs if "_R2_001.fastq" in i]
		if len(r2) == 0:
			return("se")
		else:
			return("pe")
		
	@staticmethod
	def convert_fastqs_to_class(fastqs):
		#
		paired_fastqs = list()
		for sample in fastqs:
			paired = Fastqs(sample[0], sample[1])
			paired_fastqs.append(paired)
			# lane = sample[0].split("_")[-3]
		all_fastqs = Lanes(paired_fastqs)
		return(all_fastqs)
		
		
class GetRun(object):
	
	def __init__(self):
		self.run = ""
		
	def get_run(self, sample_sheet):
		self.run = sample_sheet.split("/")[5][19:-4]
		return(self.run)
		

class LaunchMetrics(object):
	#
	def __init__(self):
		self.bam = ""
		
	def launch_metrics(self, all_samples, run):
		#
		# print(all_samples)
		for sample in all_samples:
			if ((sample.recipe == "HumanWholeGenome") or (sample.recipe == "DLP")):
				continue
			else:
				sample_params = self.get_params(sample.genome, sample.recipe)
				# if sample_params["TYPE"] == "RNA":
					# call separate RNA routine
					# continue
				bams_by_lane = self.alignment_to_genome(sample, run, sample_params)
				picard_data = self.launch_picard(bams_by_lane, run, sample, sample_params)
			
		
	@staticmethod
	def get_params(genome, recipe):
		#
		parameter_placement = list
		recipe_and_genome = ["--recipe", recipe, "--species", genome]
		sample_params = scripts.generate_run_params.main(recipe_and_genome)
		print(sample_params)
		# print(type(sample_params))
		return(sample_params)
	
		
		
	@staticmethod
	def alignment_to_genome(sample, run, sample_params):
		#
		# BIG_NODES = "-m \"is01 is02 is03 is04 is05 is06 is07 is08\" -n 60 -M 8 "
		work_dir = "/igo/staging/stats/" + run
		make_work_dir = "mkdir -p " + work_dir
		print(make_work_dir)
		call(make_work_dir, shell = True)
		os.chdir(work_dir)
		# print(all_samples)
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
			bsub_bwa_mem = "bsub -J bwa_mem___" + fastq_by_lane + " -o " + "bwa_mem___" + fastq_by_lane + ".out -n 40 -M 8 " + bwa_mem
			print(bsub_bwa_mem)
			call(bsub_bwa_mem, shell = True)
		return(bams_by_lane)	
	
	
	@staticmethod	
	def launch_picard(bams_by_lane, run, sample, sample_params):
		# 
		# BIG_NODES = " -m \"is01 is02 is03 is04 is05 is06 is07 is08\" -n 60 -M 8 "
		PICARD = "java -Dpicard.useLegacyParser=false -jar /igo/home/igo/resources/picard2.23.2/picard.jar "
		# os.chdir("/Users/naborsd/MSKCC/FASTQ/WORK/")
		# print(bams_by_lane)
		
		prjct = sample.project.split("_")[1]
		metric_file = run + "___P" + prjct + "___" + sample.sample_id + "___" + sample_params["GTAG"]
		
		# merge bams
		bsub_merge =  "bsub -w \"ended(bwa_mem___" + sample.sample_id + "*)\" -J merge_bam_files___" + sample.sample_id + " -o merge_bam_files___" + sample.sample_id + ".out -n 40 -M 8 "
		merge_bams = PICARD + "MergeSamFiles --OUTPUT " + sample.sample_id + ".merged.bam " + " ".join("--INPUT " + i for i in bams_by_lane)
		bsub_merge_bams = bsub_merge + merge_bams
		print(bsub_merge_bams)
		call(bsub_merge_bams, shell = True)
		#
		# add or replace read groups
		add_or_replace = PICARD + "AddOrReplaceReadGroups --SORT_ORDER coordinate --CREATE_INDEX true --INPUT " + sample.sample_id + ".merged.bam  " + "--OUTPUT " + sample.sample_id + ".bam  " + "--RGID " + sample.sample_id + "  --RGLB " + sample.sample_id + " --RGPL illumina  --RGPU CUSTOM-BAM  --RGSM " + sample.sample_id + " --RGCN GCL@MSKCC"
		bsub_add_or_replace = "bsub -J add_or_replace_read_groups___" + sample.sample_id + " -o " + "add_or_replace_read_groups___" + sample.sample_id + ".out -w \"done(merge_bam_files___" + sample.sample_id + ")\" -n 40 -M 8 " + add_or_replace
		print(bsub_add_or_replace)
		call(bsub_add_or_replace, shell = True)
		#
		# mark duplicates
		mark_dup = PICARD + "MarkDuplicates --CREATE_INDEX true --METRICS_FILE " + metric_file + "___MD.txt  " + "--OUTPUT " + sample.sample_id + "___MD.bam  " + "--INPUT " + sample.sample_id + ".bam"
		bsub_mark_dup = "bsub -J mark_dup___" + sample.sample_id + " -o " + "mark_dup___" + sample.sample_id + ".out -w \"done(add_or_replace_read_groups___" + sample.sample_id + ")\" -n 40 -M 8 " + mark_dup
		print(bsub_mark_dup)
		call(bsub_mark_dup, shell = True)
		#
		# alignment summary
		alignment = PICARD + "CollectAlignmentSummaryMetrics --REFERENCE_SEQUENCE " + sample_params["REFERENCE"] + " --INPUT " + sample.sample_id + "___MD.bam  " + "--OUTPUT " + metric_file + "___AM.txt"
		bsub_alignment = "bsub -J alignment_summary___" + sample.sample_id + " -o " + "alignment_summary___" + sample.sample_id + ".out -w \"done(mark_dup___" + sample.sample_id + ")\" -n 8 -M 8 " + alignment
		print(bsub_alignment)
		call(bsub_alignment, shell = True)
	
		#check for HS or RNA
		if sample_params["TYPE"] == "RNA":
			rnaseq = PICARD + "CollectRnaSeqMetrics --RIBOSOMAL_INTERVALS " + sample_params["RIBOSOMAL_INTERVALS"] + " --STRAND_SPECIFICITY NONE --REF_FLAT " + sample_params["REF_FLAT"] + "  --INPUT " + sample.sample_id + "___MD.bam  " + "--OUTPUT " + metric_file + "___RNA.txt"
			bsub_rnaseq = "bsub -J rnaseq___" + sample.sample_id + " -o " + "rnaseq___" + sample.sample_id + ".out -w \"done(mark_dup___" + sample.sample_id + ")\" -n 8 -M 8 " + rnaseq
			print(bsub_rnaseq)
			call(bsub_rnaseq, shell = True)
			
		if ("BAITS" in sample_params.keys()):
			hs_metrics = PICARD + "CollectHsMetrics --INPUT " + sample.sample_id + "___MD.bam  " + " --OUTPUT " + metric_file + "___HS.txt" +  " --REFERENCE_SEQUENCE " + sample_params["REFERENCE"] + " --BAIT_INTERVALS " + sample_params["BAITS"] + " --TARGET_INTERVALS " + sample_params["TARGETS"]
			bsub_hs_metrics = "bsub -J hs_metrics___" + sample.sample_id + " -o " + "hs_metrics___" + sample.sample_id + ".out -w \"done(mark_dup___" + sample.sample_id + ")\" -n 8 -M 8 " + hs_metrics
			print(bsub_hs_metrics)
			call(bsub_hs_metrics, shell = True)
			
			
	
def main(sample_sheet):

	# sample_sheet = ["/Users/naborsd/MSKCC/FASTQ/SampleSheets/SampleSheet_211216_MICHELLE_0466_BH3WT5DSX3.csv", "/Users/naborsd/MSKCC/FASTQ/SampleSheets/SampleSheet_220111_RUTH_0056_AH3T2TDMXY.csv", "/Users/naborsd/MSKCC/FASTQ/SampleSheets/SampleSheet_211001_MICHELLE_0444_AHCWGNDSX2.csv", "/Users/naborsd/MSKCC/FASTQ/SampleSheets/SampleSheet_220221_AYYAN_0112_000000000-K4P42.csv", "/Users/naborsd/MSKCC/FASTQ/SampleSheets/SampleSheet_220221_HERC_0811_HIPHOP1973.csv", "/Users/naborsd/MSKCC/FASTQ/SampleSheets/SampleSheet_220317_RUTH_0078_BHGLK7DSX3.csv", "/Users/naborsd/MSKCC/FASTQ/SampleSheets/SampleSheet_220317_RUTH_0079_AHGLGWDSX3.csv"]
	
	get_data = GetSampleData()
	launch_metrics = LaunchMetrics()
	get_run = GetRun()
	
	run = get_run.get_run(sample_sheet)
	all_samples = get_data.get_samples(sample_sheet, run)
	all_metrics = launch_metrics.launch_metrics(all_samples, run)

	# TODO copy txt files to DONE folder and update ngsstats database and LIMS
	# upload_stats_cmd = "RUNNAME={} /igo/work/igo/igo-demux/scripts/upload_stats.sh".format(sequencer_and_run)
        # subprocess.run(upload_stats_cmd, shell=True)
	
	# TODO email that stats have completed
	
	#print(len(all_samples))
	#for s in all_samples:
		#print(s.sample_id)
		#rint(s.project)
		#for l in s.all_fastqs.lanes:
			#print(l.r1)
			#print(l.r2)
					
############# MAIN ROUTINE
if __name__ == "__main__":
	# grab the sample sheet as an argument
	sample_sheet = sys.argv[1]
	main(sample_sheet)
