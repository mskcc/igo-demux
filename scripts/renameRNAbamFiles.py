#!/usr/bin/env python3

import re
import sys
import glob	
import os
import shutil


def main(rna_dir):
	#
	os.chdir(rna_dir)
	
	work_dir = rna_dir[:-5]
	
	txt_files = list()
	
	txt_files.append(list(glob.iglob("*WGS.txt")))
	txt_files.append(list(glob.iglob("*MD.txt")))
	txt_files.append(list(glob.iglob("*AM.txt")))
	txt_files.append(list(glob.iglob("*RNA.txt")))
	txt_files = txt_files[0] + txt_files[1] + txt_files[2] + txt_files[3]
	
	for txt_file in txt_files:
		shutil.move(txt_file, work_dir)
		
	# rename the bam files 
	for bam in glob.iglob("*.bam"):
		renamedBam = bam.split("___")[2] + bam[-4:]
		os.rename(bam, renamedBam)
		
		
		
############# MAIN ROUTINE #########
if __name__ == "__main__":
	
	rna_dir = sys.argv[1]
	
	main(rna_dir)