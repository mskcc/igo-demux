import pandas
import pytest

class SampleSheet:
    def __init__(self, path_to_samplesheet):
        # skip the header and read only data rows in the this dataframe - the [Data] section
        df = pandas.read_csv(path_to_samplesheet,skiprows=17)  
        # set of all recipes in the sample sheet
        self.recipe_set = set(df['Sample_Well'].tolist())  # Sample_Well column has the recipe, convert it to a set
        # for dual barcode sample sheets concat "index" and "index2" columns
        index_list = df['index'].tolist()
        index2_list = df['index2'].tolist()
        self.barcode_list = [a + b for a, b in zip(index_list, index2_list)]
        
        # TODO 
        self.read_lengths = [] # such as 151,151 for a PE run
    
    def needToSplitSampleSheet(self):
        # TODO DLP, PED-PEG, 10X, etc.
        return False

    def splitSampleSheet(self, path_to_write):
        return [] # TODO split and return the list of sample sheets created


def test_recipe_set():
    x = SampleSheet("test/SampleSheet.csv")
    assert ("DLP" in x.recipe_set)

def test_barcode_list():
    x = SampleSheet("test/SampleSheet.csv")
    assert ("AAGGACATAACCCCGT" in x.barcode_list)
    