from scripts.GetSampleData import GetSampleData
from scripts.GetRun import GetRun
from scripts.LaunchMetrics import LaunchMetrics
from scripts.PostStatsToIgodbAndLims import PostStatsToIgodbAndLims
import sys
import os

	
def main(args_entered):
	
	# Initiate objects
	get_run = GetRun()
	get_data = GetSampleData()
	launch_metrics = LaunchMetrics()
	post_data = PostStatsToIgodbAndLims()
	
	ss_or_project_dir = args_entered[0]
	print(ss_or_project_dir)
	
	# let's get the sequencing run
	run = get_run.get_run(ss_or_project_dir)
	print(run)
	
	if len(args_entered) == 3:
		recipe = args_entered[1]
		genome = args_entered[2]
		all_samples = get_data.get_samples_project_dir(ss_or_project_dir, run, recipe, genome)
	else:
		all_samples = get_data.get_samples_ss(ss_or_project_dir, run)

	# lets start the alignment and other metrics
	launch_metrics.launch_metrics(all_samples, run)
	
	# post data text files
	post_data.post_data_files(run)
	
	
############# MAIN ROUTINE
if __name__ == "__main__":
	
	# grab the sample sheet as an argument
	args_entered = sys.argv
	
	# take off first element of this list
	# args_entered.pop(0)
	
	main(args_entered)
	
	
	
	
	
	
	
	
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
	
	
	