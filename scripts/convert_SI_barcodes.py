import csv
import pandas as pd
import numpy as np
import argparse
import os
import glob
from xml.etree import ElementTree
from subprocess import call
from scripts.cellranger_indexes import *

def convert_SI_barcodes(tenx_ss):
  """ function to convert SI barcodes from sample sheet to the 10X quad barcodes from the cellranger_indexes.py script """
  
  # create new data frame for special sample sheet for the quad barcodes
  quad_ss_data = pd.DataFrame(columns = header)
  
  for ss_record in tenx_ss_data:
    # get the quad from the imported variables
    quad_list = vars()[tenx_ss_data["index"]]
    for quad_index in quad_list:
      quad_ss_data.loc[x] = self.ss_df_data.loc[x]
      quad_ss_data["index"].loc[x] = quad_index
    
  

  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
#  def get_run_data(run_info_location):
#	# add RunInfo.xml constant
#	run_info_file = run_info_location + '/RunInfo.xml'
#	tree = ElementTree.parse(run_info_file)
#	root = tree.getroot()
#	reads_tag = list()
#	# get the run data!
#	for child in root.iter('Read'):
#		read_data = list()
#		read_data.append(child.attrib.get('Number', None))
#		read_data.append(int(child.attrib.get('NumCycles', None)))
#		read_data.append(child.attrib.get('IsIndexedRead', None))
#		reads_tag.append(read_data)
#	use_bases_mask = ' --use-bases-mask y' + str(reads_tag[0][1] - 1) + 'n,i' + str(reads_tag[1][1]) + ',n' + str(reads_tag[2][1]) + ',y' + str(reads_tag[3][1] - 1) + 'n '
#	return use_bases_mask
