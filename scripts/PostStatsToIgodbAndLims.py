#!/usr/bin/env python3

import os
from subprocess import call
import sys
import csv
import pickle
from dataclasses import dataclass
from collections import OrderedDict
import glob
import shutil
import pathlib
import scripts.generate_run_params

class PostStatsToIgodbAndLims:
	# let's gather the txt data files and move them
	
	def post_data_files(self, run):
		#
		sequencer = run.split("_")[0]
		done_directory = "/igo/stats/DONE/{}/".format(sequencer)
		work_directory = "/igo/staging/stats/{}".format(run)
	
		os.chdir(work_directory)
		
		# move all data file to the work directory
		move_all_data_files = "/igo/work/nabors/tools/venvpy3/bin/python3 /igo/work/igo/igo-demux/scripts/move_all_data_files.py {}".format(done_directory)
		bsub_move_all_data_files = "bsub -J move_all_data_files___{0} -o move_all_data_files___{0}.out -w \"ended({0}___*)\" -cwd \"{1}\" -n 2 -M 8 {2}".format(run, work_directory, move_all_data_files)
		print(bsub_move_all_data_files)
		call(bsub_move_all_data_files, shell = True)

		# move all data file to the work directory
		push_data_to_ngs_and_lims = "/igo/work/nabors/tools/venvpy3/bin/python3 /igo/work/igo/igo-demux/scripts/push_data_to_ngs_and_lims.py {} {}".format(sequencer, run)
		bsub_push_data_to_ngs_and_lims = "bsub -K -J push_data_to_ngs_and_lims___{0} -o push_data_to_ngs_and_lims___{0}.out -w \"ended(move_all_data_files___{0})\" -cwd \"{1}\" -n 2 -M 8 {2}".format(run, work_directory, push_data_to_ngs_and_lims)
		print(bsub_push_data_to_ngs_and_lims)
		call(bsub_push_data_to_ngs_and_lims, shell = True)


	
