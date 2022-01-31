from SampleSheet import SampleSheet

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
    assert(path1.endswith("_V1.csv"))
    assert(path2.endswith("DLP.csv") or path3.endswith("_DLP.csv"))
    assert("Lane" not in ss_list[2].df_ss_data.columns)  # Confirm "Lane" column was removed
    assert(path2.endswith("10X.csv") or path3.endswith("_10X.csv"))
    assert(len(ss_list) == 4)

# Test when a sample sheet is only DLP lane information is removed and it is demuxed with "NoLaneSplitting" option in the sample sheet
def test_remove_lane_information_only_DLP():
    x = SampleSheet("test/MICHELLE_420_ONLY_DLP.csv")
    ss_list = x.split_sample_sheet()
    assert(len(ss_list) == 1)
    assert("Lane" not in ss_list[0].df_ss_data.columns)  # Confirm "Lane" column was removed

def test_remove_lane_information():
    x = SampleSheet("test/DIANA_0434.csv")
    x.remove_lane_information()
    assert(len(x.df_ss_data) == 4)

def test_remove_sample_prefix():
    x = SampleSheet("test/DIANA_0434.csv")
    print(x.project_dict)
    x.remove_sample_prefix()
    #TODO complete