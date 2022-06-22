from SampleSheet import SampleSheet
from scripts.convert_SI_barcodes import *
import pytest

def test_read_10X_sample_sheet():
    samplesheet = SampleSheet("test/SampleSheet_10X_SI.csv")
    corrected = convert_SI_barcodes(samplesheet)
    print(corrected.df_ss_data)