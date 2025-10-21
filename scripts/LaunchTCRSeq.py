#!/usr/bin/env python3

import os
from subprocess import call
import sys
import glob
from dataclasses import dataclass, replace
import pathlib
import pox

# class to launch all Metrics
class LaunchTCRSeq(object):
	#
	def __init__(self):
		
		self.combined_fastq_directory = ""
		self.tcrseq_manifest_files_present = True
		
	def launch_tcrseq_analysis(self, run, tcr_projects_to_launch):
		"""
		"""	
		
		for tcr_project in tcr_projects_to_launch:
			#
			project_fastq_directory = "/igo/staging/FASTQ/{}/{}".format(run, tcr_project.project_id)
			
			tcrseq_manifest_files_present = True
			
			# grab the Manifests files for projects
			tcrseq_manifest_files_path = "/rtssdc/mohibullahlab/LIMS/TCRseqManifest/*{}*".format(tcr_project.project_id)
			tcrseq_manifest_files = glob.glob(tcrseq_manifest_files_path)
			
			# check to see if project has TCRSeq manifest file
			if len(tcrseq_manifest_files) == 0:
				tcrseq_manifest_files_present = False
				no_tcrseq_manifest_file_message = "!!! NO TCRSeq Manifest files for {}".format(tcr_project.project_id)
				print(no_tcrseq_manifest_file_message)
				continue
			else:
				# copy the manifest files
				self.copy_tcrseq_manifest(tcr_project.project_id, tcrseq_manifest_files_path)
				# this will combine the fastq data to run
				self.combined_fastq_directory = self.combine_tcrseq_fastq_data(tcr_project.project_id)
				os.chdir("/igo/work/igo/TCR_preprocessing")
				launch_trcseq_analysis = "/igo/work/igo/mini3/bin/conda run -p /igo/work/igo/mini3/envs/TCRSeq python /igo/work/igo/TCR_preprocessing/manager.py -p {} -m /igo/work/nabors/TCRSeq/Manifest/{} -species {}".format(self.combined_fastq_directory, tcr_project.project_id, tcr_project.genome)
				print(launch_trcseq_analysis)
				call(launch_trcseq_analysis, shell = True)
			
			
			
	@staticmethod
	def copy_tcrseq_manifest(tcr_project_id, tcrseq_manifest_files_path):
		"""
		"""
		
		# create directory for the manifest files
		tcrseq_analysis_manifest_directory = "/igo/staging/TCRSeq_Analysis/manifest/{}".format(tcr_project_id)
		pathlib.Path(tcrseq_analysis_manifest_directory).mkdir(parents = True, exist_ok = True)
		set_permissions_on_tcrseq_manifest_files = "chmod -R 775 {}".format(tcrseq_analysis_manifest_directory)
		
		# copy manifest files to the analysis directory
		copy_tcrseq_manifest_file = "cp -rv {} {}".format(tcrseq_manifest_files_path, tcrseq_analysis_manifest_directory)
		print(copy_tcrseq_manifest_file)
		call(copy_tcrseq_manifest_file, shell = True)
		
		return()
	
	
	
	@staticmethod
	def combine_tcrseq_fastq_data(tcr_project_id):
		"""
		"""
		# let's look for this project in another run folder
		total_runs_for_project = pox.find(tcr_project_id, recurse = 3, root = "/igo/staging/FASTQ")
	
		valid_runs_for_project = [x for x in total_runs_for_project if ("REFERENCE" not in x)]
		
		combined_fastq_directory = ""
		
		if (len(valid_runs_for_project) > 1):
			combined_fastq_directory = "/igo/staging/TCRSeq_Analysis/combined_fastqs_for_TCRSeq_analysis/{}".format(tcr_project_id)
			pathlib.Path(combined_fastq_directory).mkdir(parents = True, exist_ok = True)
			#
			for valid_run_for_project in valid_runs_for_project:
				copy_fastqs_to_combined_directory = "cp -rv {}/*/*.fastq.gz {}".format(valid_run_for_project, combined_fastq_directory)
				print(copy_fastqs_to_combined_directory)
				call(copy_fastqs_to_combined_directory, shell = True)
			#
		else:
			combined_fastq_directory = valid_runs_for_project[0]	
			
		return(combined_fastq_directory)
	
	
	
	
	
	
	
	
	
	
	
	
	