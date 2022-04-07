from SampleSheet import SampleSheet
import pytest

def test_read_empty_sample_sheet():
    x = SampleSheet("test/empty_sample_sheet.csv")
    print("Success")

def test_read_blank_sample_sheet():
    with pytest.raises(Exception):
        x = SampleSheet("test/blank_sample_sheet.csv")

def test_read_10X_sample_sheet():
    x = SampleSheet("test/SampleSheet_10X.csv")
    print("Success")

def test_read_SE_sample_sheet():
    x = SampleSheet("test/SampleSheet_PEPE.csv")
    print("Success")

def test_WGS_only_not_split():
    x = SampleSheet("test/DIANA_0434.csv")
    ss_list = x.split_sample_sheet()
    assert(len(ss_list) == 1)

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
    
def test_split():
    x = SampleSheet("test/SampleSheet_DLP.csv")
    ss_list = x.split_sample_sheet()
    path0 = ss_list[0].path
    path1 = ss_list[1].path
    path2 = ss_list[2].path
    path3 = ss_list[3].path
    assert(path0.endswith("_REFERENCE.csv"))
    assert(path1.endswith(".csv"))
    assert(path2.endswith("DLP.csv") or path3.endswith("_DLP.csv"))
    assert("Lane" in ss_list[2].df_ss_data.columns) 
    assert(path2.endswith("10X.csv") or path3.endswith("_10X.csv"))
    assert(len(ss_list) == 4)

# Test when a sample sheet is only DLP lane information is removed and it is demuxed with "NoLaneSplitting" option in the sample sheet
def test_only_DLP_split():
    x = SampleSheet("test/MICHELLE_420_ONLY_DLP.csv")
    ss_list = x.split_sample_sheet()
    assert(len(ss_list) == 1)
    assert("Lane" in ss_list[0].df_ss_data.columns)

def test_remove_sample_prefix():
    x = SampleSheet("test/DIANA_0434.csv")
    print(x.project_dict)
    x.remove_sample_prefix()
    #TODO complete