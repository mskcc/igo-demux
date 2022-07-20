#!/usr/bin/env python3

import os
import sys
import argparse
from xml.etree import ElementTree


def getSequencingReadData(sequencerPath):
	""" method to check the read lengths in the RunInfo.xml file """
	
	# compare tag
	atacReads = [51, 8, 24, 50]
	
	# set variables to parge the XML file
	runInfoFile = sequencerPath + '/RunInfo.xml'
	tree = ElementTree.parse(runInfoFile)
	root = tree.getroot()
	readsTag = list()
	
	# get the run data!
	for child in root.iter('Read'):
		readData = list()
		readData.append(child.attrib.get('Number', None))
		readData.append(int(child.attrib.get('NumCycles', None)))
		readData.append(child.attrib.get('IsIndexedRead', None))
		readsTag.append(readData)
		
	# create a list to check the run data	
	detectedReads = [y[1] for y in readsTag]
	
	# check the read lengths.  if they equal the atacReads list, then set atac variable to True
	if (detectedReads == atacReads):
		atac = True
		
	return(atac)
		
		
def main(sequencerPath):
	
	atac = getSequencingReadData(sequencerPath)
	
	return(atac)
	
	
		
if __name__ == '__main__':
	
	# grab the sequncer run path  as an argument
	sequencerPath = sys.argv[1]
	
	main(sequencerPath)