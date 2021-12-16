import pandas
import pytest

class SampleSheet:
    def __init__(self, path_to_samplesheet):
        df = pandas.read_csv(path_to_samplesheet,skiprows=17)  # skip the header and read only data rows in the this dataframe
        # set of all recipes in the sample sheet
        self.recipe_set = set(df['Sample_Well'].tolist())  # Sample_Well column has the recipe, convert it to a set
        # TODO for dual barcode could concat "index" and "index2" columns
        self.barcode_list = [] # list of all barcodes in the sample sheet with dual barcodes separated by a "-"
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