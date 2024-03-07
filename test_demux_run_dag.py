from SampleSheet import SampleSheet
import demux_run_dag

def test_WGS_only_not_split():
    x = SampleSheet("/home/runner/work/igo-demux/igo-demux/test/SampleSheet_220304_MICHELLE_0485_BHFN7NDSX3.csv")
    cmd_set = demux_run_dag.build_dragen_cmds(x, "MICHELLE_0485_BHFN7NDSX3")
    assert(len(cmd_set) == 7)


def test_get_dlp_chip():
    # Test that the DLP chip returned is correct even when the run has multiple DLP projects with different chip IDs
    sample_sheet = SampleSheet("/home/runner/work/igo-demux/igo-demux/test/SampleSheet_220412_MICHELLE_0501_BHFNH5DSX3_DLP.csv")
    for project in sample_sheet.project_set:
        chip_id = demux_run_dag.get_dlp_chip(sample_sheet, project)
        if project == "13098":
            assert(chip_id == "128676A")
        if project == "13098_C":
            assert(chip_id == "128680A")


def test_get_dlp_chip_from_sample_name():
    assert(demux_run_dag.get_dlp_chip_from_sample_name("1_cellcover_4C_128749A_37_42") == "128749A")