#!/usr/bin/env python3

# class to grab the run information from the sample sheet
class GetRun(object):
	
	def __init__(self):
		self.run = ""
		
	def get_run(self, ss_or_project_dir):
		if ss_or_project_dir[-3:] == "csv":
			self.run = ss_or_project_dir.split("/")[5][19:-4]
		else:
			self.run = ss_or_project_dir.split("/")[4]
		return(self.run)
	