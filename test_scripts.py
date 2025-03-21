import scripts.cellranger as cellranger
import scripts.get_total_reads_from_demux
import scripts.move_failed_fastqs
from SampleSheet import SampleSheet
import pytest

def test_get_project_id():
    project_id = scripts.move_failed_fastqs.get_project_id("13820_B_1_1_1_1_1")
    assert(project_id == "13820_B")

def test_get_base_igo_id():
    igo_id = scripts.move_failed_fastqs.get_base_igo_id("13820_B_1_1_1_1_1")
    assert(igo_id == "13820_B_1")

def testCellranger_generate_cellranger_cmd():
    sample_ID_list = ["06265_8869_1_IGO_06265_AG_3","Third-Transcriptome_IGO_11969_E_3", "Second_IGO_11969_E_2"]
    fastq_file_list_dict = {'06265_8869_1_IGO_06265_AG_3': ['/igo/staging/FASTQ/DIANA_0453_AHFKJ5DRXY/Project_06265_AG/Sample_06265_8869_1_IGO_06265_AG_3'], 'Third-Transcriptome_IGO_11969_E_3': ['/igo/staging/FASTQ/DIANA_0450_AH3JL3DSX3/Project_11969_E/Sample_Third-Transcriptome_IGO_11969_E_3', '/igo/staging/FASTQ/DIANA_0454_BH555MDMXY/Project_11969_E/Sample_Third-Transcriptome_IGO_11969_E_3'], 'Second_IGO_11969_E_2': ['/igo/staging/FASTQ/DIANA_0453_AHFKJ5DRXY/Project_11969_E/Sample_Second_IGO_11969_E_2', '/igo/staging/FASTQ/DIANA_0450_AH3JL3DSX3/Project_11969_E/Sample_Second_IGO_11969_E_2']}
    genome_dict = {"06265_8869_1_IGO_06265_AG_3":"Human","Third-Transcriptome_IGO_11969_E_3":"Human_GeneticallyModified", "Second_IGO_11969_E_2":"Mouse"}
    cmd = []
    for sample in sample_ID_list:
        if genome_dict[sample] != "Human" and genome_dict[sample] != "Mouse":
            genome_dict[sample] = "Mouse"
        cmd.append(cellranger.generate_cellranger_cmd(sample, "count", genome_dict[sample], fastq_file_list_dict[sample], "DIANA_0453_AHFKJ5DRXY"))
    test_result = ["bsub -J DIANA_0453_AHFKJ5DRXY_Project_06265_AG_06265_8869_1_IGO_06265_AG_3_count_cellranger -o DIANA_0453_AHFKJ5DRXY_Project_06265_AG_06265_8869_1_IGO_06265_AG_3_count_cellranger.out /igo/work/nabors/tools/cellranger-9.0.1/cellranger count --id=Sample_06265_8869_1_IGO_06265_AG_3__count --transcriptome=/igo/work/nabors/genomes/10X_Genomics/GEX/refdata-gex-GRCh38-2020-A --fastqs=/igo/staging/FASTQ/DIANA_0453_AHFKJ5DRXY/Project_06265_AG/Sample_06265_8869_1_IGO_06265_AG_3 --create-bam=true --nopreflight --jobmode=lsf --mempercore=64 --disable-ui --maxjobs=200",   
    "bsub -J DIANA_0453_AHFKJ5DRXY_Project_11969_E_Third-Transcriptome_IGO_11969_E_3_count_cellranger -o DIANA_0453_AHFKJ5DRXY_Project_11969_E_Third-Transcriptome_IGO_11969_E_3_count_cellranger.out /igo/work/nabors/tools/cellranger-9.0.1/cellranger count --id=Sample_Third-Transcriptome_IGO_11969_E_3__count --transcriptome=/igo/work/nabors/genomes/10X_Genomics/GEX/refdata-gex-mm10-2020-A --fastqs=/igo/staging/FASTQ/DIANA_0450_AH3JL3DSX3/Project_11969_E/Sample_Third-Transcriptome_IGO_11969_E_3,/igo/staging/FASTQ/DIANA_0454_BH555MDMXY/Project_11969_E/Sample_Third-Transcriptome_IGO_11969_E_3 --create-bam=true --nopreflight --jobmode=lsf --mempercore=64 --disable-ui --maxjobs=200",
    "bsub -J DIANA_0453_AHFKJ5DRXY_Project_11969_E_Second_IGO_11969_E_2_count_cellranger -o DIANA_0453_AHFKJ5DRXY_Project_11969_E_Second_IGO_11969_E_2_count_cellranger.out /igo/work/nabors/tools/cellranger-9.0.1/cellranger count --id=Sample_Second_IGO_11969_E_2__count --transcriptome=/igo/work/nabors/genomes/10X_Genomics/GEX/refdata-gex-mm10-2020-A --fastqs=/igo/staging/FASTQ/DIANA_0453_AHFKJ5DRXY/Project_11969_E/Sample_Second_IGO_11969_E_2,/igo/staging/FASTQ/DIANA_0450_AH3JL3DSX3/Project_11969_E/Sample_Second_IGO_11969_E_2 --create-bam=true --nopreflight --jobmode=lsf --mempercore=64 --disable-ui --maxjobs=200"]
    
    for i in range (3): 
        assert(cmd[i] == test_result[i])

def testCellranger_get_tag():
    assert(cellranger.get_tag("10X_genomic") == "Skip")
    assert(cellranger.get_tag("SC_Chromium-GEX-3") == "count")

def testCellranger_get_sequencer_runID():
    fastq_path = "/igo/staging/FASTQ/DIANA_0453_AHFKJ5DRXY/Project_06265_AG/Sample_06265_8869_1_IGO_06265_AG_3"
    assert(cellranger.get_sequencer_runID(fastq_path) == ("diana", "DIANA_0453_AHFKJ5DRXY"))

def testGettotalreads():
    sample_list = ["PDX_WD0010_P1_1845_IGO_12754_E_1", "PDX_WD0010_P1_1850_IGO_12754_E_2"]
    total_reads_dict = scripts.get_total_reads_from_demux.get_total_reads(sample_list, "/home/runner/work/igo-demux/igo-demux/test/Demultiplex_Stats.csv")
    print(total_reads_dict)
    assert(total_reads_dict["PDX_WD0010_P1_1845_IGO_12754_E_1"] == 770373032)
    assert(total_reads_dict["PDX_WD0010_P1_1850_IGO_12754_E_2"] == 602357556)

def testGettotalreadsDLP():
    sample_sheet = SampleSheet("/home/runner/work/igo-demux/igo-demux/test/SampleSheet_DLP_multiprojects.csv")
    total_reads_dict = scripts.get_total_reads_from_demux.get_total_reads_DLP(sample_sheet, "/home/runner/work/igo-demux/igo-demux/test/Demultiplex_Stats_DLP.csv" )
    print(total_reads_dict)
    assert(total_reads_dict["Project_11113_L"]["samples"][1] == 3802998466)
    assert(total_reads_dict["Project_11113_L"]["pos_control"][1] == 654555718)
    assert(total_reads_dict["Project_11113_L"]["neg_control"][1] == 247002)
    assert(total_reads_dict["Project_11113_M"]["pos_control"][1] == 599399798)
