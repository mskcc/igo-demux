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

# setting up the data classes for the sample sheet structure for launching the metrics
@dataclass
class Fastqs:
	r1: str
	r2: str
	
@dataclass
class Lanes:
	lanes: list()
	
@dataclass
class Sample:
	sample_id: str
	genome: str
	recipe: str
	project: str
	all_fastqs: str
	
# this list contains the headers of the columns.  we will access the data using these listings
data_headers = list()

class GetSampleData:
	#
	def __init__(self):
		self.all_sample_ids = list()
		self.all_samples = list()
		self.duplicate_sample = ""
		self.all_fastqs = list()
		
		
	# grab the project from either the staging or delivery FASTQ location
	def get_samples_project_dir(self, project_directory, run, recipe, genome):
	
		global data_headers
			
		sample_ids = os.listdir(project_directory)
			
		for sample_id in sample_ids:
			#
			# go to a routine to pair the reads.  and return them
			self.all_lanes = self.get_fastqs(self, sample_id, project_directory)
			# get project
			project = project_directory.split("/")[-1]
			sr = Sample(sample_id[7:], genome, recipe, project, self.all_lanes)
			self.all_samples.append(sr)
		return(self.all_samples)
		
		
	# let's grab the sample sheet and run info to store into the data classes
	# the sample sheet is in csv file format
	def get_samples_ss(self, sample_sheet, run):
		
		global data_headers
		
		csv_sample_data = list()
		testing = list()
		with open(sample_sheet) as csv_file:
			csv_reader = csv.reader(csv_file, delimiter = ",")
			got_data = False
			# once we find the "Lane" header, let's store the header row and the data
			for row in csv_reader:
				if (row[0] == "Lane"):
					got_data = True
					data_headers = row
					print(data_headers)
				elif (row[0] != "Lane") and got_data:
					# do not process hWGS, DLP, 10X samples or MissionBio - they have their own processes - Do not PROCESS any POOLEDNORMAL Samples
					if ("POOLEDNORMAL" in row[data_headers.index("Sample_Project")]):
						continue
					self.all_sample_ids.append(row[data_headers.index("Sample_ID")])
					# check for duplicate samples in the sample sheet
					self.duplicate_sample = self.check_this_sample(row[data_headers.index("Sample_ID")], self.all_sample_ids)
					if not self.duplicate_sample:
						# go to a routine to pair the reads.  and return them
						project_directory = "/igo/staging/FASTQ/{}/{}".format(run, row[data_headers.index("Sample_Project")])
						sample_id = "Sample_{}".format(row[data_headers.index("Sample_ID")])
						self.all_lanes = self.get_fastqs(self, sample_id, project_directory)
						sr = Sample(row[data_headers.index("Sample_ID")], row[data_headers.index("Sample_Plate")], row[data_headers.index("Sample_Well")], row[data_headers.index("Sample_Project")], self.all_lanes)
						self.all_samples.append(sr)
				else:
					continue
		return(self.all_samples)
	
	
	# checking for duplicate reads in the sample sheet
	@staticmethod
	def check_this_sample(sample_id, all_sample_ids):
		#
		x = 0
		duplicate_sample = False
		x = all_sample_ids.count(sample_id)
		if (x > 1):
			duplicate_sample = True
		return(duplicate_sample)
	
	# let's start the process of obtaining the fastqs	
	@staticmethod
	def get_fastqs(self, sample_id, project_directory):
		#
		# get run from the sample sheet
		fastq_directory = "{}/{}/".format(project_directory, sample_id)
		fastqs  = os.listdir(fastq_directory)
		# right here, let's eliminate all fastqs from the list that isn't R1 or R2
		fastqs = [x for x in fastqs if ("_R1_001.fastq" in x) or ("_R2_001.fastq" in x)]
		# check the run type: PE or SE
		run_type = self.determine_run_type(fastqs)
		# put the fastqs in order by lane
		fastqs = self.pair_fastqs_by_lane(fastqs, run_type)
		# put the fastqs in the Fastq and Lane data classes
		fastqs = self.convert_fastqs_to_class(fastqs)
		return(fastqs)
	
	# pairing the fastqs according to lane
	@staticmethod
	def pair_fastqs_by_lane(fastqs, run_type):
		""" PAIR FASTQS OF THE SAMPLE BY LANE """
		# grab the number of lanes used
		if run_type == "se":
			pairs = list()
			for fq in fastqs:
				se_pair = [fq, None]
				pairs.append(se_pair)
		else:
			lane_size = int(len(fastqs) / 2)
			pairs = list()
			for i in range(0, lane_size, 1):
				fastq1 = fastqs[0][:-15]
				fastq_pair = [i for i in fastqs if fastq1 in i]
				fastq_pair.sort()
				fastqs.remove(fastq_pair[0]), fastqs.remove(fastq_pair[1])
				pairs.append(fastq_pair)
			pairs.sort()
		return(pairs)
	
	# get that run type: PE or SE
	@staticmethod
	def determine_run_type(fastqs):
		r2 = [i for i in fastqs if "_R2_001.fastq" in i]
		if len(r2) == 0:
			return("se")
		else:
			return("pe")
		
	# storing the fastqs in to the data class
	@staticmethod
	def convert_fastqs_to_class(fastqs):
		#
		paired_fastqs = list()
		for sample in fastqs:
			paired = Fastqs(sample[0], sample[1])
			paired_fastqs.append(paired)
		all_fastqs = Lanes(paired_fastqs)
		return(all_fastqs)
	
	