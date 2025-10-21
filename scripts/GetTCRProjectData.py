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

	
@dataclass
class TCRProject:
	project_id: str
	genome: str
	
# this list contains the headers of the columns.  we will access the data using these listings
data_headers = list()

class GetTCRProjectData:
	#
	def __init__(self):
		
		self.all_tcr_projects = list()

	# let's grab the sample sheet and run info to store into the data classes
	# the sample sheet is in csv file format
	def get_tcr_projects(self, sample_sheet):
		
		TCRSEQ_RECIPES = ["TCRSeq", "TCR Sequencing", "TCRseq-IGO-alpha", "TCRseq-IGO-beta", "TCR_IGO"]
		
		tcr_projects_and_genomes = set()
		
		with open(sample_sheet) as csv_file:
			csv_reader = csv.reader(csv_file, delimiter = ",")
			got_data = False
			# once we find the "Lane" header, let's store the header row and the data
			for row in csv_reader:
				if (row[0] == "Lane"):
					got_data = True
				elif (row[0] != "Lane") and got_data:
					# let's go thru the sample sheet and grab the TCRSeq related samples
					if any(recipe in row[3] for recipe in TCRSEQ_RECIPES):
						tcr_projects_and_genomes.add((row[7], row[2]))
								
		for tcr_entry in tcr_projects_and_genomes:
			tcr_to_dataclass_record = TCRProject(tcr_entry[0], tcr_entry[1])
			self.all_tcr_projects.append(tcr_to_dataclass_record)
						
		return(self.all_tcr_projects)
	
	

	


	
	
	