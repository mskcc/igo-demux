from SampleSheet import SampleSheet
import demux_run_dag

def test_WGS_only_not_split():
    x = SampleSheet("test/SampleSheet_220304_MICHELLE_0485_BHFN7NDSX3.csv")
    cmd_set = demux_run_dag.build_dragen_cmds(x, "MICHELLE_0485_BHFN7NDSX3")
    assert(len(cmd_set) == 7)


def test_get_dlp_chip():
    sample_sheet = SampleSheet("test/MICHELLE_420_ONLY_DLP.csv")
    result = demux_run_dag.get_dlp_chip(sample_sheet)
    assert("110IO_DLP_UNSORTED_110720" == result)