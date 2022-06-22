import pandas as pd
from SampleSheet import SampleSheet
from scripts.cellranger_indexes import *

def convert_SI_barcodes(samplesheet):
    """ function to convert SI barcodes from sample sheet to the 10X quad barcodes from the cellranger_indexes.py script """
    """ from here, we will need to convert the variables so it can be used in SampleSheet.py """
    
    print("Converting 10X samplesheet with SI barcodes to their real barcodes")

    # create new data frame for special sample sheet for the quad barcodes
    quad_ss_data = pd.DataFrame(columns=samplesheet.df_ss_data.columns.values)
    
    # row_position will make sure we will skip down to the correct rows when creating the new sample sheet rows
    row_position = 0
    for x in range(0, len(samplesheet.df_ss_data["index"]), 1):
        # get the quad from the imported variables
        si_barcode = samplesheet.df_ss_data["index"].loc[x].replace("-", "_")
        quad_list = globals()[si_barcode]  # lookup "SI-" barcode in the global variable list
        # loop thru the quad set of barcodes and use these to replace the SI barcodes
        for y in range(0, len(quad_list), 1):
            quad_ss_data.loc[row_position] = samplesheet.df_ss_data.loc[x]
            quad_ss_data["index"].loc[row_position] = quad_list[y]
            row_position += 1

    return SampleSheet(samplesheet.df_ss_header, quad_ss_data, samplesheet.path)
  

  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
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
