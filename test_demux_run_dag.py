from SampleSheet import SampleSheet
import demux_run_dag

def test_WGS_only_not_split():
    x = SampleSheet("test/SampleSheet_220304_MICHELLE_0485_BHFN7NDSX3.csv")
    cmd_set = demux_run_dag.build_dragen_cmds(x, "MICHELLE_0485_BHFN7NDSX3")
    assert(len(cmd_set) == 7)


def test_get_dlp_chip():
    sample_sheet = SampleSheet("test/SampleSheet_220412_MICHELLE_0501_BHFNH5DSX3_DLP.csv")
    result = []
    for project in sample_sheet.project_set:
        result.append(demux_run_dag.get_dlp_chip(sample_sheet, project))
        
    assert(["128676A", "128676A", "128680A", "128680A"] == result)