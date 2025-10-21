#!/usr/bin/env python3
from scripts.GetTCRProjectData import GetTCRProjectData
from scripts.GetRun import GetRun
from scripts.LaunchTCRSeq import LaunchTCRSeq
import os
import sys


def main(argv):
	"""
	"""
	
	# Initiate objects
	get_run = GetRun()
	get_tcr_data = GetTCRProjectData()
	launch_tcr = LaunchTCRSeq()
	
	# grab the run name from the sample sheet
	run = get_run.get_run(argv)
	
	# get the list of the TCRSeq projects in this run
	tcr_projects_to_launch = get_tcr_data.get_tcr_projects(argv)
	
	#check to see if we had TCRSeq projects
	if len(tcr_projects_to_launch) == 0:
		print("NO TCRSeq Projects in run {}".format(run))
		sys.exit()
	
	# launch all tcr projects associated with this run
	launch_tcr.launch_tcrseq_analysis(run, tcr_projects_to_launch)
	


############# MAIN ROUTINE
if __name__ == "__main__":
	
	sample_sheet = sys.argv[1]
	
	main(sys.argv)
