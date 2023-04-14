from GetSampleData import GetSampleData
from GetRun import GetRun
from LaunchMetrics import LaunchMetrics
from PostStatsToIgodbAndLims import PostStatsToIgodbAndLims
import sys
import os

	
def main(sample_sheet):
	
	# Initiate objects
	get_run = GetRun()
	get_data = GetSampleData()
	launch_metrics = LaunchMetrics()
	post_data = PostStatsToIgodbAndLims()
	
	print(sample_sheet)
	
	# let's get the sequencing run
	run = get_run.get_run(sample_sheet)
	print(run)
	
	# let process all of those samples
	all_samples = get_data.get_samples(sample_sheet, run)
	print(all_samples)
	
	# lets start the alignment and other metrics
	launch_metrics.launch_metrics(all_samples, run)
	
	# post data text files
	post_data.post_data_files(run)
	
	
############# MAIN ROUTINE
if __name__ == "__main__":
	
	# grab the sample sheet as an argument
	sample_sheet = sys.argv[1]
	
	main(sample_sheet)

	
	
	