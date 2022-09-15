#!/usr/bin/env python3

import os
import sys
import argparse
from xml.etree import ElementTree


def get_sequencing_read_data(sequencer_path):
	""" method to check the read lengths in the RunInfo.xml file """
	
	# compare tag
	atac_reads = [[51, 8, 24, 50], [51, 8, 16, 50]]
	
	
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
	if ((detected_reads == atac_reads[0]) or (detected_reads == atac_reads[1])):
		atac = True
	else:
		atac = False
		
	return(atac)
		
		
def main(sequencer_path):
	
	atac = get_sequencing_read_data(sequencer_path)
	
	return(atac)
	
	
		
if __name__ == '__main__':
	
	# grab the sequncer run path  as an argument
	sequencer_path = sys.argv[1]
	
	main(sequencer_path)