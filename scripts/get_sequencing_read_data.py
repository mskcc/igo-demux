#!/usr/bin/env python3

import os
import sys
import argparse
from xml.etree import ElementTree


def get_sequencing_read_data(sequencer_path):
	""" method to check the read lengths in the RunInfo.xml file """
	
	# compare tag
	atac_reads = [[51, 8, 24, 50], [51, 8, 16, 51], [51, 8, 24, 51]]
	
	
	# set variables to parge the XML file
	run_info_file = sequencer_path + '/RunInfo.xml'
	tree = ElementTree.parse(run_info_file)
	root = tree.getroot()
	reads_tag = list()
	
	# get the run data!
	for child in root.iter('Read'):
		read_data = list()
		read_data.append(child.attrib.get('Number', None))
		read_data.append(int(child.attrib.get('NumCycles', None)))
		read_data.append(child.attrib.get('IsIndexedRead', None))
		reads_tag.append(read_data)
		
	# create a list to check the run data	
	detected_reads = [y[1] for y in reads_tag]
	
	# check the read lengths.  if they equal the atacReads list, then set atac variable to True
	# if ((detected_reads == atac_reads[0]) or (detected_reads == atac_reads[1]) or (detected_reads == atac_reads[2]))
	# making a simpler if statement to check for reads
	if (detected_reads in atac_reads):
		atac = True
		use_bases_mask = "Y" + str(reads_tag[0][1]) + ",I" + str(reads_tag[1][1]) + ",Y" + str(reads_tag[2][1]) + ",Y" + str(reads_tag[3][1])
	else:
		atac = False
		use_bases_mask = [reads_tag[0][1], reads_tag[-1][1]]
		
	return(atac, use_bases_mask)
		
		
def main(sequencer_path):
	
	atac, use_bases_mask = get_sequencing_read_data(sequencer_path)
	
	return(atac, use_bases_mask)
	
	
		
if __name__ == '__main__':
	
	# grab the sequncer run path  as an argument
	sequencer_path = sys.argv[1]
	
	main(sequencer_path)