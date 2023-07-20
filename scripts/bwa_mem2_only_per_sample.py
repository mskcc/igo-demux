from __future__ import print_function

import re
import sys
import os
import glob
import argparse

from os.path import join
from os.path import basename
from os.path import abspath
from os.path import isdir
from subprocess import call

# Define some paths and other constants up here

FASTQ_REGEX = r"(([_.]R{0}[_.].+)|([_.]R{0}\.)|(_{0}\.))f(ast)?q(\.gz)?$"

def format_file_name(file_name):
	"""Return destination file name."""
	file_name = basename(file_name)
	for index in [1, 2]:
		if re.search(FASTQ_REGEX.format(index), file_name):
			suffix = "_{}.fastq".format(index)
			letter_index_fastq = r"[_.]R{}([_.])?\.f(ast)?q".format(index)
			number_index_fastq = r"[_.]{}([_.])?\.f(ast)?q".format(index)
			letter_index_any_location = r"[_.]R{}[_.]".format(index)
			file_name = re.sub(letter_index_fastq, ".fastq", file_name)
			file_name = re.sub(number_index_fastq, ".fastq", file_name)
			file_name = re.sub(letter_index_any_location, "_", file_name)
			file_name = re.sub(r"[_.]f(ast)?q", suffix, file_name)

	return(file_name)

def symlink_input_data(outdir, sequencing_data):
	"""Symlink fastq files using pcapcore expected naming."""
	links_dir = join(outdir, "input_data")
	sym_links = []

	if not isdir(links_dir):
		os.makedirs(links_dir)

	for src in map(abspath, sequencing_data):
		dst = join(links_dir, format_file_name(basename(src)))
		sym_links.append(dst)

		try:
			os.unlink(dst)
			os.symlink(src, dst)
		except OSError:
			os.symlink(src, dst)

	return sym_links

# ###########  WGS script to run Juan Medina's container for the BAM

def submit_jobs_to_cluster(sample_directory, work_directory, new_fastqs):
	
	# grab sample_id and the run
	sample_id = sample_directory.split("/")[6][7:]
	run = sample_directory.split("/")[4]
	
	print(sample_id)
	print(run)
	
	generate_ppg_bam = "/igo/work/nabors/tools/bwamem2/bwa_mem2.pl -fragment 10 -reference /igo/work/nabors/tools/bwamem2/ref/gr37.fasta -threads 32 -map_threads 32 -sample {} -outdir {} {}".format(sample_id, work_directory, " ".join(new_fastqs))
	bsub_generate_ppg_bam = "bsub -J {0}___BWA_MEM2___{1} -o {0}___BWA_MEM2___{1}.out -m \"is01 is02 is03 is04 is05 is06 is07 is08\" -n 32 -M 8 {2}".format(run, sample_id, generate_ppg_bam)
	print(bsub_generate_ppg_bam)
	call(bsub_generate_ppg_bam, shell = True)



def main(sample_directory, work_directory):
	
	#change to sample directory
	os.chdir(sample_directory)
	
	# grab the sample fastqs
	sample_fastqs = glob.glob("*fastq.gz")
	print(sample_fastqs)
	new_fastqs = symlink_input_data(work_directory, sample_fastqs)
	new_fastqs.sort()
	# go ahead and start the biild of the PEDPEG bams
	submit_jobs_to_cluster(sample_directory, work_directory, new_fastqs)
	
	

############# MAIN ROUTINE
if __name__ == '__main__':
	
	sample_directory = sys.argv[1]
	work_directory = sys.argv[2]
	
	main(sample_directory, work_directory)
		
	
	