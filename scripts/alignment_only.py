from scripts.GetSampleData import GetSampleData
from scripts.GetRun import GetRun
from scripts.LaunchAlignment import LaunchAlignment
import sys
import os

	
def main(argv):
	
	# Initiate objects
	get_run = GetRun()
	get_data = GetSampleData()
	launch_alignment = LaunchAlignment()
	
	# check the type for input, then proceed accordingly
	if type(argv) is list:
		project_directory = argv[0]
		recipe = argv[1]
		genome = argv[2]
		run = get_run.get_run(project_directory)
		all_samples = get_data.get_samples_project_dir(project_directory, run, recipe, genome)
	else:
		project_directory = argv
		run = get_run.get_run(argv)
		all_samples = get_data.get_samples_ss(argv, run)
		
	# lets start the alignment and other metrics
	launch_alignment.launch_alignment(all_samples, run, project_directory)
	
	
############# MAIN ROUTINE
if __name__ == "__main__":
	
	# grab the sample sheet or project directory, recipe and genome
	main(sys.argv)
	
	

	
	