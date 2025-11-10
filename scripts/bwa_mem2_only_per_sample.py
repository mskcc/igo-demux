#!/usr/bin/env python3

import sys
import os
import glob
import argparse
from subprocess import call

# ###########  WGS script to run Juan Medina's container for the BAM

def retrieve_downsampled_fastqs(sample_id, project_id, dragen_directory):
	"""
	"""
	# grab downsampled fastqs
	downsampled_fastq_listing = "{}/{}*.fastq".format(dragen_directory, sample_id)
	downsampled_fastqs = glob.glob(downsampled_fastq_listing)
	downsampled_fastqs.sort()

	print(downsampled_fastq_listing)
	print(downsampled_fastqs)
	
	# rename those fastqs we just located
	downsampled_fastq_r1 = "{}/{}.downsample_1.fastq".format(dragen_directory, sample_id)
	downsampled_fastq_r2 = "{}/{}.downsample_2.fastq".format(dragen_directory, sample_id)
	print(downsampled_fastq_r1)
	print(downsampled_fastq_r2)
	
	os.rename(downsampled_fastqs[0], downsampled_fastq_r1)
	os.rename(downsampled_fastqs[1], downsampled_fastq_r2)

	output_directory = "{}/{}".format(dragen_directory[:-7], project_id)
	pathlib.Path(output_directory).mkdir(parents = True, exist_ok = True)

	return(downsampled_fastq_r1, downsampled_fastq_r2, output_directory)


def submit_jobs_to_cluster(downsampled_fastq_r1, downsampled_fastq_r2, output_directory, dragen_job_name_header, sample_id):
	"""
	"""
	run = output_directory.split("/")[5]
	
	generate_ppg_bam = "/igo/work/nabors/tools/bwamem2/bwa_mem2.pl -fragment 10 -reference /igo/work/nabors/tools/bwamem2/ref/gr37.fasta -threads 32 -map_threads 32 -sample {} -outdir {} {} {}".format(sample_id, output_directory, downsampled_fastq_r1, downsampled_fastq_r2)
	bsub_generate_ppg_bam = "bsub -J {0}___BWA_MEM2_for_PPG___{1} -o {0}___BWA_MEM2_for_PPG___{1}.out -n32 -M8 {2}".format(run, sample_id, generate_ppg_bam)
	print(bsub_generate_ppg_bam)
	call(bsub_generate_ppg_bam, shell = True)

	return()


def main(sample_directory, work_directory):
	"""
	"""
	# grab downsampled fastqs
	downsampled_fastq_r1, downsampled_fastq_r2, output_directory = retrieve_downsampled_fastqs(sample_id, project_id, dragen_directory)
	
	# go ahead and start the biild of the PEDPEG bams
	submit_jobs_to_cluster(downsampled_fastq_r1, downsampled_fastq_r2, output_directory, dragen_job_name_header, sample_id)
	
	
############# MAIN ROUTINE
if __name__ == '__main__':
	
	sample_id = sys.argv[1]
	project_id = sys.argv[2]
	dragen_directory = sys.argv[3][:-1]
	dragen_job_name_header = sys.argv[4]
	
	main(sample_id, project_id, dragen_directory, dragen_job_name_header)
		
	
	
