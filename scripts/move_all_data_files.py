#!/usr/bin/env python3

import glob
import shutil
import sys


def main(done_directory):
	
	for data_file in glob.iglob("*.txt"):
		shutil.move(data_file, done_directory)
	
	

if __name__ == "__main__":
	
	done_directory = sys.argv[1]
	
	main(done_directory)








	