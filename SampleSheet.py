import pandas
import re
import os
from copy import deepcopy
from scripts.cellranger_indexes import *

"""
Reads an IGO LIMS generated sample sheet .csv and splits the sample sheet if necessary to generate sample sheets ready for 
Illumina DRAGEN demuxes with the correct options set for 10X & DLP samples.
"""
class SampleSheet:

    """
    Overloaded constructor either should have 1 argument with the path to a sample sheet
    or three arguments, the [Header] data frame, the [Data] data frame and the path
    """
    def __init__(self, *args):
        if len(args) == 1: # this is the path to the sample sheet
            self.path = args[0]
            self.read_csv(self.path)
        else:
            self.df_ss_header = args[0]
            self.df_ss_data = args[1]
            self.path = args[2]

        # dictionary of project->recipe - NOTE, not accurate for MICHELLE_0485 where 08822_PC has HumanWholeGenome & RNASeq_RiboDeplete
        self.project_dict = pandas.Series(self.df_ss_data['Sample_Well'].values,index=self.df_ss_data['Sample_Project']).to_dict()
        self.sample_dict = pandas.Series(self.df_ss_data['Sample_Well'].values,index=self.df_ss_data['Sample_ID']).to_dict()
        self.project_set = set(self.df_ss_data['Sample_Project'].tolist())  # Sample_Project column has the projects, convert it to a set

        # set of all recipes in the sample sheet
        self.recipe_set = set(self.df_ss_data['Sample_Well'].tolist())  # Sample_Well column has the recipe, convert it to a set
        # for dual barcode sample sheets concat "index" and "index2" columns
        index_list = self.df_ss_data['index'].tolist()
        if 'index2' in self.df_ss_data.columns:  # check if this is a dual-index run
            index2_list = self.df_ss_data['index2'].tolist()
            self.barcode_list = [a + b for a, b in zip(index_list, index2_list)]
        else:
            self.barcode_list = index_list
        
        barcode_10X = re.compile("SI*")
        # list of all special "SI-*" 10X barcodes
        self.barcode_list_10X = list(filter(barcode_10X.match, self.barcode_list))
        
        # list of index read lengths such as [151,151] for a PE run
        self.read_lengths = []
        index_read1 = int(self.df_ss_header.iat[9,0])
        self.read_lengths.append(index_read1)
        if type(self.df_ss_header.iat[10,0]) != float:
            index_read2 = int(self.df_ss_header.iat[10,0])
            self.read_lengths.append(index_read2)
    
    def read_csv(self, path_to_samplesheet):
        # skip the header and read only data rows in the this dataframe - the [Data] section
        # find row of sample sheet which has the [Data] section
        line_number = 0
        sheet = open(path_to_samplesheet,"r")
        with open(path_to_samplesheet, 'r') as read_obj:
            for line in read_obj:
                line_number += 1
                if "[Data]" in line:
                    break
        sheet.close()
        if line_number <= 1:
            raise Exception("Sample sheet is blank")
        print("[Data] section of sample sheet detected on line: {}".format(line_number))
        self.df_ss_header = pandas.read_csv(path_to_samplesheet,nrows=line_number-1)
        self.df_ss_data = pandas.read_csv(path_to_samplesheet,skiprows=line_number)

    def write_csv(self):
        print("Saving sample sheet to " + self.path)

        # pandas dataframe has blank column headers, make them write correctly to the .csv
        csv = open(self.path, "w")
        csv.write("[Header],,,,,,,,,\n")
        csv.close()

        self.df_ss_header.to_csv(self.path, mode='a',index=False,header=False)
        self.df_ss_data.to_csv(self.path, mode='a',index=False)


    """
    Creates a new [Data] section without 'Lane' information for each sample.
    """
    def remove_lane_information(self):
        #DRAGEN "--no-lane-splitting" requires a sample sheet without lane information
        print("Removing sample sheet lane information.")
        self.df_ss_data.drop_duplicates(inplace=True) 

    """
    Returns a list of sample sheets from splitting the original or the original sample sheet
    if no splitting of samples for demux is necessary
    """
    def split_sample_sheet(self):
        """
        10X 
         if barcodes start with 'SI' like 'SI-NA-C7' take just the barcodes starting with 'SI' and remove index2 called "_10X"
        DLP
         if sample sheet recipes have mixed DLP and other all DLP need to go on a separate sample sheet named "_DLP"
        """
        # if 10x DRAGEN demux add to header CreateFastqForIndexReads,1,,,,,,,
        if any("10X_" in s for s in self.recipe_set):
            self.df_ss_header.loc[len(self.df_ss_header.index)-1] = ["CreateFastqForIndexReads",1,"","","","","","",""]
            self.df_ss_header.loc[len(self.df_ss_header.index)] = ["[Data]","","","","","","","",""]
            print("Added CreateFastqForIndexReads,1 to sample sheet header since 10X samples are present")

        ss_copy = deepcopy(self)

        split_ss_list = [ss_copy, self]  # result list starts with the original sample sheet

        was_split = False
        if "DLP" in self.recipe_set and len(self.recipe_set) > 1:
            print("Copying all DLP samples to a new sample sheet")
            # copy all DLP rows to a new sample sheet
            dlp_data = self.df_ss_data[self.df_ss_data["Sample_Well"].str.match("DLP") == True].copy()
            # and remove DLP samples from the main sample sheet
            self.df_ss_data= self.df_ss_data[self.df_ss_data["Sample_Well"].str.match("DLP") == False].copy()
            # rename DLP sample sheet w/"_DLP.csv"
            dlp_path = os.path.splitext(self.path)[0]+'_DLP.csv'
            header_copy = self.df_ss_header.copy(deep=True)
            dlp_ss = SampleSheet(header_copy, dlp_data, dlp_path)
            split_ss_list.append(dlp_ss)
            was_split = True

        # check if sample sheet has 'SI-*' barcodes and normal barcodes
        if len(self.barcode_list_10X) > 0 and len(self.barcode_list) != len(self.barcode_list_10X):
            print("Copying all 10X SI- barcodes to new sheet and remove index2 column")
            print("Non-DRAGEN demux, must have Sample_ID column with Sample_ prefix")
            tenx_data = self.df_ss_data[ self.df_ss_data["index2"].str.match('^SI-.*') == True ].copy()
            rest_data = self.df_ss_data[ self.df_ss_data["index2"].str.match('^SI-.*') == False ].copy()
            self.df_ss_data = rest_data
            tenx_path = os.path.splitext(self.path)[0]+'_10X.csv'
            # if ATAC because read length is 51,50 () for example DIANA_427 must use cellranger-ATAC mkfastq 
            tenx_ss = SampleSheet(self.df_ss_header, tenx_data, tenx_path)
            # convert SI barcodes to their real barcodes
            tenx_ss_real_barcodes = convert_SI_barcodes(tenx_ss)
            split_ss_list.append(tenx_ss_real_barcodes)
            was_split = True

        if was_split:
            # Rename the original sample sheet
            split_ss_list[0].path = os.path.splitext(self.path)[0]+'_REFERENCE.csv'
            split_ss_list[1].path = os.path.splitext(self.path)[0]+'.csv'
        else:
            split_ss_list = [ss_copy]

        return split_ss_list

def convert_SI_barcodes(samplesheet):
    """ function to convert SI barcodes from sample sheet to the 10X quad barcodes from the cellranger_indexes.py """
    
    print("Converting 10X samplesheet with SI barcodes to their real barcodes")

    # create new data frame for special sample sheet for the quad barcodes
    quad_ss_data = pandas.DataFrame(columns=samplesheet.df_ss_data.columns.values)
    
    # row_position will make sure we will skip down to the correct rows when creating the new sample sheet rows
    row_position = 0
    for x in range(0, len(samplesheet.df_ss_data["index"])):
        # get the quad from the imported variables
        si_barcode = samplesheet.df_ss_data["index"].loc[x].replace("-", "_")
        quad_list = globals()[si_barcode]  # lookup "SI-" barcode in the global variable list
        # loop thru the quad set of barcodes and use these to replace the SI barcodes
        for y in range(0, len(quad_list)):
            quad_ss_data.loc[row_position] = samplesheet.df_ss_data.loc[x]
            quad_ss_data["index"].loc[row_position] = quad_list[y]
            row_position += 1

    quad_ss_data = quad_ss_data.drop(columns=['index2'])
    return SampleSheet(samplesheet.df_ss_header, quad_ss_data, samplesheet.path)