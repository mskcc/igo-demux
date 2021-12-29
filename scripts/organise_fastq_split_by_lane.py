import re
import sys
import glob	
import os
import linecache
from subprocess import call

run = sys.argv[1]

x = 0

os.chdir(run)

for project in glob.iglob("Project_*"):
	# print(project)
	project_dir = run + "/" +  project
	os.chdir(project_dir)
	for fastq in glob.iglob("*fastq.gz"):
		sample_header = "Sample_"
		sample_name = re.search("(.+?)_IGO_(.+?)_(\d+)", fastq)[0]
		fastq_folder = sample_header + sample_name
		fastqs_folders = next(os.walk("."))[1]
		if fastq_folder not in fastqs_folders:
			os.mkdir(fastq_folder, 0o775)
			print(fastq_folder)
			all_fastqs_of_one_sample = sample_name + "*" + ".fastq.gz"
			move_fastq_2_folder = "mv " + all_fastqs_of_one_sample + " " + fastq_folder
			print(move_fastq_2_folder)
			x += 1
			print(x)
			call(move_fastq_2_folder, shell = True)

	

# regex. regular expression
# ^ = start at beginning of string. $ = start at end of string
# findall gives the output in a list, so [0] gets the 0th element out of the list
# look below
# SAMPLE = re.findall(r'^(.+?)___',FILE)[0]

