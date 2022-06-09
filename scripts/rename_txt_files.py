import re
import sys
import glob	
import os
# import linecache
import pickle
import shutil


def main(rna_dir):
	#
	os.chdir(rna_dir)
	
	work_dir = rna_dir[:-4]
	
	run = rna_dir.split("/")[4]
	
	rna_samples_and_gtags_file = open("rna_samples_and_gtags.pickle", "rb")
	rna_samples_and_gtags = pickle.load(rna_samples_and_gtags_file)

	gap = "___"
	
	txt_files = list()
	
	txt_files.append(list(glob.iglob("*WGS.txt")))
	txt_files.append(list(glob.iglob("*MD.txt")))
	txt_files.append(list(glob.iglob("*AM.txt")))
	txt_files = txt_files[0] + txt_files[1] + txt_files[2]

	for txt_file in txt_files:
		#
		sample = re.findall(r"^(.+?)___", txt_file)[0]
		project = re.findall("IGO_(.+?)_\d", sample)[0]
		suffix = re.findall(r"___(.+?)$", txt_file)[0]
		for rna_sample_and_gtag in rna_samples_and_gtags:
			if (sample == rna_sample_and_gtag[0]):
				renamed_txt_file = run + gap + "P" + project + gap + sample + gap + rna_sample_and_gtag[1] + gap + suffix
				print(txt_file + "   >>>   " + renamed_txt_file )
				os.rename(txt_file, renamed_txt_file)
				shutil.move(renamed_txt_file, work_dir)
			else:
				continue
			
	# mv RNA	 txt files
	for rna_txt_file in glob.iglob("*RNA.txt"):
		shutil.move(rna_txt_file, work_dir)
		
	

############# MAIN ROUTINE #########
if __name__ == "__main__":
	
	rna_dir = sys.argv[1]
	
	main(rna_dir)

	
