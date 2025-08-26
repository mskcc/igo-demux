#!/usr/bin/env python3
from GetTCRProjectData import GetTCRProjectData
from GetRun import GetRun
from LaunchTCRSeq import LaunchTCRSeq
import os
import sys


def main(sample_sheet):
	"""
	"""
	
	# Initiate objects
	get_run = GetRun()
	get_tcr_data = GetTCRProjectData()
	launch_tcr = LaunchTCRSeq()
	
	# grab the run name from the sample sheet
	run = get_run.get_run(sample_sheet)
	
	# get the list of the TCRSeq projects in this run
	tcr_projects_to_launch = get_tcr_data.get_tcr_projects(sample_sheet)
	
	#check to see if we had TCRSeq projects
	if len(tcr_projects_to_launch) == 0:
		print("NO TCRSeq Projects in run {}".format(run))
		sys.exit()
	
	# launch all tcr projects associated with this run
	launch_tcr.launch_tcrseq_analysis(run, tcr_projects_to_launch)
	


############# MAIN ROUTINE
if __name__ == "__main__":
	
	# sample_sheet = "/igo/work/igo/SampleSheetCopies/SampleSheet_250124_FAUCI_0264_A22GCM5LT3.csv"
	# sample_sheet = "/igo/work/igo/SampleSheetCopies/SampleSheet_250808_BONO_0043_A22V5WLLT4.csv"
	# sample_sheet = "/igo/work/igo/SampleSheetCopies/SampleSheet_250725_BONO_0037_A22NWJ7LT4.csv"
	# sample_sheet = "/igo/work/igo/SampleSheetCopies/SampleSheet_250821_AMELIE_0084_AAGMHCCM5.csv"
	sample_sheet = "/igo/work/igo/SampleSheetCopies/SampleSheet_250822_BONO_0048_ABCDEFGHIJ.csv"
	# sample_sheet = sys.argv[1]
	
	main(sample_sheet)