import pandas
import re
import os
import pytest

class SampleSheet:
    
    def __init__(self, df_header, df_data, path_to_samplesheet):
        if df_header.empty:
            return  # call read_from_file to add data

        self.df_ss_header = df_header
        self.df_ss_data = df_data
        self.path = path_to_samplesheet

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
    
    def read_from_file(self, path_to_samplesheet):
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
        df_ss_header = pandas.read_csv(path_to_samplesheet,nrows=line_number-1)
        df_ss_data = pandas.read_csv(path_to_samplesheet,skiprows=line_number)

        return SampleSheet(df_ss_header, df_ss_data, path_to_samplesheet)

    def need_to_split_sample_sheet(self):
        # if DLP is mixed with anything
        if "DLP" in self.recipe_set and len(self.recipe_set) > 1:
            print("Sample sheet must be split due to DLP and other recipes")
            return True

        if len(self.barcode_list_10X) > 0 and len(self.barcode_list) != len(self.barcode_list_10X):
            print("Sample sheet must be split due to 10X barcodes")
            return True

        # TODO PED-PEG

        return False

    def write_csv(self, path_to_write):
        print("Saving sample sheet to " + path_to_write)

        # pandas dataframe has blank column headers, make them write correctly to the .csv
        csv = open(path_to_write, "w")
        csv.write("[Header],,,,,,,,,\n")
        csv.close()

        self.df_ss_header.to_csv(path_to_write, mode='a',index=False,header=False)
        self.df_ss_data.to_csv(path_to_write, mode='a',index=False)

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
        if self.need_to_split_sample_sheet() == False:
            return [self]

        # Reference https://github.com/mskcc/nf-fastq-plus/blob/master/bin/create_multiple_sample_sheets.py

        split_ss_list = [self]
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
            print("Copying all 10X SI-barcodes to new sheet and remove index2 column")
            print("Non-DRAGEN demux, must have Sample_ID column with Sample_ prefix")
            tenx_data = self.df_ss_data[ self.df_ss_data["index2"].str.match('^SI-.*') == True ].copy()
            rest_data = self.df_ss_data[ self.df_ss_data["index2"].str.match('^SI-.*') == False ].copy()
            self.df_ss_data = rest_data
            tenx_path = os.path.splitext(self.path)[0]+'_10X.csv'
            # TODO insert line in [Settings]
            # if ATAC because read length is 51,50 () for example DIANA_427 must use cellranger-ATAC mkfastq 
            # and add to the [Header] options for correct DRAGEN demux with index fastqs
            tenx_ss = SampleSheet(self.df_ss_header, tenx_data, tenx_path)
            split_ss_list.append(tenx_ss)
            

        # TODO if 10x DRAGEN demux add to header CreateFastqForIndexReads,1,,,,,,, 

        return split_ss_list

def test_barcode_write():
    x = SampleSheet(pandas.DataFrame(),pandas.DataFrame(),"").read_from_file("test/SampleSheet.csv")
    x.write_csv("test/SampleSheetCopy.csv")

def test_barcode_read_lengths():
    x = SampleSheet(pandas.DataFrame(),pandas.DataFrame(),"").read_from_file("test/SampleSheet.csv")
    assert (x.read_lengths[0] == 151)
    assert (x.read_lengths[1] == 151)

def test_recipe_set():
    x = SampleSheet(pandas.DataFrame(),pandas.DataFrame(),"").read_from_file("test/SampleSheet.csv")
    assert ("DLP" in x.recipe_set)

def test_barcode_list():
    x = SampleSheet(pandas.DataFrame(),pandas.DataFrame(),"").read_from_file("test/SampleSheet.csv")
    assert ("AAGGACATAACCCCGT" in x.barcode_list)

def test_need_to_split_sample_sheet():
    x = SampleSheet(pandas.DataFrame(),pandas.DataFrame(),"").read_from_file("test/SampleSheet_DLP.csv")
    assert(x.need_to_split_sample_sheet() == True)
    
def test_split():
    # TODO add more 
    x = SampleSheet(pandas.DataFrame(),pandas.DataFrame(),"").read_from_file("test/SampleSheet_DLP.csv")
    ss_list = x.split_sample_sheet()
    print(ss_list[2].df_ss_data)
    assert(len(ss_list) == 3)