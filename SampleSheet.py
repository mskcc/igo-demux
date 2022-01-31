import pandas
import re
import os
from copy import deepcopy

"""
Reads an IGO LIMS generated sample sheet .csv and splits the sample sheet if necessary to generate sample sheets ready for 
Illumina DRAGEN demuxes with the correct options set for 10X, WGS & DLP samples.
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

        # dictionary of project->recipe
        self.project_dict = pandas.Series(self.df_ss_data['Sample_Well'].values,index=self.df_ss_data['Sample_Project']).to_dict()
        self.project_set = set(self.df_ss_data['Sample_Project'].tolist())  # Sample_Project column has the projects, convert it to a set

        # set of all recipes in the sample sheet
        self.recipe_set = set(self.df_ss_data['Sample_Well'].tolist())  # Sample_Well column has the recipe, convert it to a set
        # for dual barcode sample sheets concat "index" and "index2" columns
        index_list = self.df_ss_data['index'].tolist()
        index2_list = self.df_ss_data['index2'].tolist()
        self.barcode_list = [a + b for a, b in zip(index_list, index2_list)]
        
        barcode_10X = re.compile("SI*")
        # list of all special "SI-*" 10X barcodes
        self.barcode_list_10X = list(filter(barcode_10X.match, self.barcode_list))
        
        # list of index read lengths such as [151,151] for a PE run
        self.read_lengths = []
        index_read1 = int(self.df_ss_header.iat[9,0])
        self.read_lengths.append(index_read1)
        if self.df_ss_header.iat[10,0].isnumeric():
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
        print("[Data] section of sample sheet detected on line: {}".format(line_number))
        assert line_number > 0
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

    def remove_sample_prefix(self):
        self.df_ss_data['Sample_ID'] = self.df_ss_data['Sample_ID'].str.replace('Sample_' , '')

    """
    Creates a new [Data] section without 'Lane' information for each sample.
    """
    def remove_lane_information(self):
        #DRAGEN "--no-lane-splitting" requires a sample sheet without lane information
        print("Removing sample sheet lane information.")
        self.df_ss_data.drop(columns=['Lane', 'Sample_Name'], inplace=True)
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
        PED-PEG & WGS
         if recipe is 'HumanWholeGenome'
        """
        # if 10x DRAGEN demux add to header CreateFastqForIndexReads,1,,,,,,,
        if any("10X_" in s for s in self.recipe_set):
            self.df_ss_header.loc[len(self.df_ss_header.index)-1] = ["CreateFastqForIndexReads",1,"","","","","","",""]
            self.df_ss_header.loc[len(self.df_ss_header.index)] = ["[Data]","","","","","","","",""]
            print("Added CreateFastqForIndexReads,1 to sample sheet header since 10X samples are present")

        # check if the sample sheet is only DLP or WGS samples
        if ("DLP" in self.recipe_set or "HumanWholeGenome" in self.recipe_set) and len(self.recipe_set) == 1:
            print("Adding NoLaneSplitting option to the sample sheet")
            self.df_ss_header.loc[len(self.df_ss_header.index)-1] = ["NoLaneSplitting","true","","","","","","",""]
            self.df_ss_header.loc[len(self.df_ss_header.index)] = ["[Data]","","","","","","","",""]
            self.remove_lane_information()

        ss_copy = deepcopy(self)

        split_ss_list = [ss_copy, self]  # result list starts with the original sample sheet

        was_split = 0
        if "DLP" in self.recipe_set and len(self.recipe_set) > 1:
            print("Copying all DLP samples to a new sample sheet")
            # copy all DLP rows to a new sample sheet
            dlp_data = self.df_ss_data[ self.df_ss_data["Sample_Well"].str.match("DLP") == True ].copy()
            # and remove DLP samples from the main sample sheet
            rest_data = self.df_ss_data[ self.df_ss_data["Sample_Well"].str.match("DLP") == False ].copy()
            self.df_ss_data = rest_data
            # rename DLP sample sheet w/"_DLP.csv"
            dlp_path = os.path.splitext(self.path)[0]+'_DLP.csv'
            header_copy = self.df_ss_header.copy(deep=True)
            header_copy.loc[len(header_copy.index)-1] = ["NoLaneSplitting","true","","","","","","",""]
            header_copy.loc[len(header_copy.index)] = ["[Data]","","","","","","","",""]
            dlp_ss = SampleSheet(header_copy, dlp_data, dlp_path)
            dlp_ss.remove_lane_information()
            split_ss_list.append(dlp_ss)
            was_split = 1

        if "HumanWholeGenome" in self.recipe_set and len(self.recipe_set) > 1:
            print("Copying all HumanWholeGenome samples to a new sample sheet")
            # copy all WGS rows to a new sample sheet
            wgs_data = self.df_ss_data[ self.df_ss_data["Sample_Well"].str.match("HumanWholeGenome") == True ].copy()
            # and remove WGS samples from the main sample sheet
            rest_data = self.df_ss_data[ self.df_ss_data["Sample_Well"].str.match("WGS") == False ].copy()
            self.df_ss_data = rest_data
            # rename WGS sample sheet w/"_wgs.csv"
            wgs_path = os.path.splitext(self.path)[0]+'_wgs.csv'
            header_copy = self.df_ss_header.copy(deep=True)
            # add no lane splitting option to the header
            header_copy.loc[len(header_copy.index)-1] = ["NoLaneSplitting","true","","","","","","",""]
            header_copy.loc[len(header_copy.index)] = ["[Data]","","","","","","","",""]
            wgs_ss = SampleSheet(header_copy, wgs_data, wgs_path)
            wgs_ss.remove_lane_information()
            split_ss_list.append(wgs_ss)
            was_split = 1

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
            split_ss_list.append(tenx_ss)
            was_split = 1

        if was_split:
            # Rename the original sample sheet
            split_ss_list[0].path = os.path.splitext(self.path)[0]+'_REFERENCE.csv'
            split_ss_list[1].path = os.path.splitext(self.path)[0]+'.csv'
        else:
            split_ss_list = [ss_copy]

        return split_ss_list