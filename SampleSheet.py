import pandas
import re
import os
from copy import deepcopy

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
        if len(index2_list) != 0:
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
        self.df_ss_data.drop(columns=['Lane', 'Sample_ID'], inplace=True)
        self.df_ss_data.sort_values("Sample_Name", inplace = True)
        self.df_ss_data.drop_duplicates(subset ="Sample_Name",keep = False, inplace = True) 

    def need_to_split_sample_sheet(self):
        # if DLP is mixed with anything
        if "DLP" in self.recipe_set and len(self.recipe_set) > 1:
            print("Sample sheet must be split due to DLP and other recipes")
            return True

        if len(self.barcode_list_10X) > 0 and len(self.barcode_list) != len(self.barcode_list_10X):
            print("Sample sheet must be split due to 10X barcodes")
            return True

        # TODO WGS & PED-PEG

        return False

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
         if project ID starts with 08822 and recipe is 'HumanWholeGenome' ask Darrell
        """
        # if 10x DRAGEN demux add to header CreateFastqForIndexReads,1,,,,,,,
        if any("10X_" in s for s in self.recipe_set):
            self.df_ss_header.loc[len(self.df_ss_header.index)-1] = ["CreateFastqForIndexReads",1,"","","","","","","",""]
            self.df_ss_header.loc[len(self.df_ss_header.index)] = ["[Data]","","","","","","","","",""]
            print("Added CreateFastqForIndexReads,1 to sample sheet header since 10X samples are present")

        if self.need_to_split_sample_sheet() == False:
            return [self]

        # Reference https://github.com/mskcc/nf-fastq-plus/blob/master/bin/create_multiple_sample_sheets.py

        ss_copy = deepcopy(self)

        split_ss_list = [ss_copy, self]
        if "DLP" in self.recipe_set and len(self.recipe_set) > 1:
            print("Copying all DLP samples to a new sample sheet")
            # copy all DLP rows to a new sample sheet
            # and copy everything else to a different sample sheet
            dlp_data = self.df_ss_data[ self.df_ss_data["Sample_Well"].str.match("DLP") == True ].copy()
            rest_data = self.df_ss_data[ self.df_ss_data["Sample_Well"].str.match("DLP") == False ].copy()
            self.df_ss_data = rest_data
            # rename DLP sample sheet w/"_DLP.csv"
            dlp_path = os.path.splitext(self.path)[0]+'_DLP.csv'
            dlp_ss = SampleSheet(self.df_ss_header, dlp_data, dlp_path)
            split_ss_list.append(dlp_ss)

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

        # Rename the original sample sheet now modified with fewer rows
        split_ss_list[0].path = os.path.splitext(self.path)[0]+'_REFERENCE.csv'
        split_ss_list[1].path = os.path.splitext(self.path)[0]+'_V1.csv'

        return split_ss_list

def test_barcode_read_lengths():
    x = SampleSheet("test/SampleSheet.csv")
    assert (x.read_lengths[0] == 151)
    assert (x.read_lengths[1] == 151)

def test_recipe_set():
    x = SampleSheet("test/SampleSheet.csv")
    assert ("DLP" in x.recipe_set)

def test_barcode_list():
    x = SampleSheet("test/SampleSheet.csv")
    assert ("AAGGACATAACCCCGT" in x.barcode_list)

def test_need_to_split_sample_sheet():
    x = SampleSheet("test/SampleSheet_DLP.csv")
    assert(x.need_to_split_sample_sheet() == True)
    
def test_split():
    x = SampleSheet("test/SampleSheet_DLP.csv")
    ss_list = x.split_sample_sheet()
    path0 = ss_list[0].path
    path1 = ss_list[1].path
    path2 = ss_list[2].path
    path3 = ss_list[3].path
    assert(path0.endswith("_REFERENCE.csv"))
    assert(path1.endswith("_V1.csv"))
    assert(path2.endswith("DLP.csv") or path3.endswith("_DLP.csv"))
    assert(path2.endswith("10X.csv") or path3.endswith("_10X.csv"))
    assert(len(ss_list) == 4)

def test_remove_lane_information():
    x = SampleSheet("test/SampleSheet.csv")
    x.remove_lane_information()
    # TODO test with sample sheet that has multiple rows and duplicates to be removed
    
