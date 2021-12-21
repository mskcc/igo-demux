import pandas
import re
import pytest

class SampleSheet:
    
    def __init__(self, path_to_samplesheet):
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

        df = pandas.read_csv(path_to_samplesheet,skiprows=line_number) 

        # set of all recipes in the sample sheet
        self.recipe_set = set(df['Sample_Well'].tolist())  # Sample_Well column has the recipe, convert it to a set
        # for dual barcode sample sheets concat "index" and "index2" columns
        index_list = df['index'].tolist()
        index2_list = df['index2'].tolist()
        self.barcode_list = [a + b for a, b in zip(index_list, index2_list)]
        
        barcode_10X = re.compile("SI*")
        self.barcode_list_10X = list(filter(barcode_10X.match, self.barcode_list))
        # TODO 
        self.read_lengths = [] # such as 151,151 for a PE run
    

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
        #TODO save file
        print("Saving sample sheet to " + path_to_write)

    def split_sample_sheet(self):
        """
        10X 
         if barcodes start with 'SI' like 'SI-NA-C7' take just the barcodes starting with 'SI' and remove index2 called "_10X"
        DLP
         if sample sheet recipes have mixed DLP and other all DLP need to go on a separate sample sheet named "_DLP"
        PED-PEG 
        & WGS
         if project ID starts with 08822 and recipe is 'HumanWholeGenome' ask Darrell
        
        """
        if self.needToSplitSampleSheet() == False:
            return []

        # Reference https://github.com/mskcc/nf-fastq-plus/blob/master/bin/create_multiple_sample_sheets.py

        if "DLP" in self.recipe_set and len(self.recipe_set) > 1:
            print("Copying all DLP samples to a new sample sheet")
            # TODO copy all DLP rows to a new sample sheet
            # TODO and copy everything else to a different sample sheet

        # sample sheet has 'SI-*' barcodes and others
        if len(self.barcode_list_10X) > 0 and len(self.barcode_list) != len(self.barcode_list_10X):
            print("Copying all 10X SI-barcodes to new sheet and remove index2 column")
            # and add to the [Header] options for correct DRAGEN demux with index fastqs

        return [] # TODO split and return the list of sample sheets created


def test_recipe_set():
    x = SampleSheet("test/SampleSheet.csv")
    assert ("DLP" in x.recipe_set)

def test_barcode_list():
    x = SampleSheet("test/SampleSheet.csv")
    assert ("AAGGACATAACCCCGT" in x.barcode_list)

def test_needToSplitSampleSheet():
    x = SampleSheet("test/SampleSheet_DLP.csv")
    assert(x.need_to_split_sample_sheet() == True)
    