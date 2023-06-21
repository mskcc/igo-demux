#!/usr/bin/env python3

import glob
from subprocess import call
import sys


def main(done_directory):
	
	for data_file in glob.iglob("*.txt"):
		copy_data_file = "cp -rv {} {}".format(data_file, done_directory)
		print(copy_data_file)
		call(copy_data_file, shell = True)
	
	
if __name__ == "__main__":
	
	done_directory = sys.argv[1]
	
	main(done_directory)








	