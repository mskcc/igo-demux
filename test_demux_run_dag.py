from SampleSheet import SampleSheet
import demux_run_dag

def test_WGS_only_not_split():
    x = SampleSheet("test/DIANA_0434.csv")
    cmd_set = demux_run_dag.build_dragen_cmds(x, "DIANA_0434_AH2V3TDSX3_WGS")
    assert(len(cmd_set) == 3)

def test_build_dragen_cmds():
    sample_sheet = SampleSheet("test/DIANA_0441_WGS.csv")
    cmd_set = demux_run_dag.build_dragen_cmds(sample_sheet, "DIANA_0441_AH2V3TDSX3_WGS")
    assert(cmd_set.pop() == "bsub -J DIANA_0441_AH2V3TDSX3_WGS_PS4268T_IGO_04540_Q_10 -eo /igo/staging/stats/DIANA_0441_AH2V3TDSX3_WGS/PS4268T_IGO_04540_Q_10.out -q dragen -n 48 -M 4 /opt/edico/bin/dragen --ref-dir /staging/ref/GRCh38_graph --enable-duplicate-marking true --enable-map-align-output true --fastq-list /igo/staging/FASTQ/DIANA_0441_AH2V3TDSX3_WGS/Reports/fastq_list.csv --output-directory /igo/staging/stats/DIANA_0441_AH2V3TDSX3_WGS --fastq-list-sample-id PS4268T_IGO_04540_Q_10 --output-file-prefix DIANA_0441_AH2V3TDSX3___P04540_Q___PS4268T_IGO_04540_Q_10")
    print(*cmd_set)

def test_get_dlp_chip():
    sample_sheet = SampleSheet("test/MICHELLE_420_ONLY_DLP.csv")
    result = demux_run_dag.get_dlp_chip(sample_sheet)
    assert("110IO_DLP_UNSORTED_110720" == result)