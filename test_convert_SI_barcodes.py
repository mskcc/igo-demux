from SampleSheet import SampleSheet,convert_SI_barcodes
import pytest

def test_read_10X_sample_sheet():
    samplesheet = SampleSheet("test/SampleSheet_10X_SI.csv")
    corrected = convert_SI_barcodes(samplesheet)
    # Sample sheet Data section should have 8 rows and include barcodes:
    # SI_GA_G9 = SI_P2_G9 = ["TAGGACGT", "ATCCCACA", "GGAATGTC", "CCTTGTAG"]
    assert(len(corrected.df_ss_data) == 8)
    #print(corrected.df_ss_data.to_string())