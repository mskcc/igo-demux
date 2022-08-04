#!/usr/bin/env python3

import os
import sys
import glob
from subprocess import call
import shutil
import time

def main(work_directory):
	
	os.chdir(work_directory)
	
	
	# move MD, AM, HS, WGS and RNA txt files to /igo/stats/DONE
	VALID_SUFFIXES = ["_MD.txt", "_AM.txt", "_HS.txt", "RNA.txt", "WGS.txt"]
	
	run = os.getcwd().split("/")[4]
	sequencer = run.split("_")[0]

	for txt_file in glob.iglob("*txt"):
		if (txt_file[-7:] in VALID_SUFFIXES):
			DONE_directory = "/igo/stats/DONE/" + sequencer + "/" + txt_file
			shutil.move(txt_file, DONE_directory)
		
	delphi_endpoint = "curl http://delphi.mskcc.org:8080/ngs-stats/picardstats/updaterun/" + sequencer + "/" + run
	print(delphi_endpoint)
	call(delphi_endpoint, shell = True)
	
	# give time for NGS to post data
	time.sleep(100)
	
	lims_endpoint = "curl -k https://igolims.mskcc.org:8443/LimsRest/updateLimsSampleLevelSequencingQc?runId=" + run
	print(lims_endpoint)
	call(lims_endpoint, shell = True)
	
	
	
	
	
############# MAIN ROUTINE
if __name__ == "__main__":
	
	work_directory = sys.argv[1]
	
	main(work_directory)
	
	
# bsub_delphi_endpoint = "bsub -J PUSH_TO_DELPHI___" + run + " -o PUSH_TO_DELPHI___" + run + ".out -w \"ended(MOVE_ALL_RNA_TXT_FILES___*)\" -n 2 -M 8 " + delphi_endpoint
# bsub_lims_endpoint = "bsub -J PUSH_TO_LIMS___" + run + " -o PUSH_TO_LIMS___" + run + ".out -w \"ended(PUSH_TO_DELPHI___*)\" -n 2 -M 8 " + lims_endpoint
