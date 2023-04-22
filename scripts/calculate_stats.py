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
		run = get_run.get_run(argv)
		all_samples = get_data.get_samples_ss(argv, run)

	# lets start the alignment and other metrics
	launch_metrics.launch_metrics(all_samples, run)
	
	# post data text files
	post_data.post_data_files(run)
	
	
############# MAIN ROUTINE
if __name__ == "__main__":
	
	# grab the sample sheet or project directory, recipe and genome
	main(sys.argv)
	
	
	
	
	
	
	
	
# if rna_or_dragen_directory_present:   # this would be set to true
# rna_or_dragen_data_files.append(list(glob.iglob("*WGS.txt")))
# rna_or_dragen_data_files.append(list(glob.iglob("*MD.txt")))
# rna_or_dragen_data_files.append(list(glob.iglob("*AM.txt")))
# rna_or_dragen_data_files = rna_or_dragen_data_files[0] + rna_or_dragen_data_files[1] + rna_or_dragen_data_files[2]
# for data_file in rna_or_dragen_data_files:
# shutil.move(data_file, done_directory)
	
# push_txt_files = "python3 /igo/work/igo/igo-demux/scripts/push_to_ngs_and_lims.py " 
# bsub_mv_all_txt = "bsub -K -J PUSH_DATA___" + run + " -o " + "PUSH_DATA___" + run + ".out -w \"ended(" + run + "___*)\" -n 2 -M 8 " + push_txt_files
# print(bsub_mv_all_txt)
# call(bsub_mv_all_txt, shell = True)
	
	
	