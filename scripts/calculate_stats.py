from scripts.GetSampleData import GetSampleData
from scripts.GetRun import GetRun
from scripts.LaunchMetrics import LaunchMetrics
from scripts.PostStatsToIgodbAndLims import PostStatsToIgodbAndLims
import sys
import os

	
def main(argv):
	
	# Initiate objects
	get_run = GetRun()
	get_data = GetSampleData()
	launch_metrics = LaunchMetrics()
	post_data = PostStatsToIgodbAndLims()
	
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
	launch_metrics.launch_metrics(all_samples, run, project_directory)
	
	# post data text files
	post_data.post_data_files(run)
	
	
############# MAIN ROUTINE
if __name__ == "__main__":
	
	# grab the sample sheet or project directory, recipe and genome
	main(sys.argv)
	
	

	