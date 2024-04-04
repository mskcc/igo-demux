# work folder
STATS_AREA = "/igo/stats/CELLRANGER/"

# config info 
ACCESS = 0o775
config_dict = {
    "count": {
        "tool": " /igo/work/nabors/tools/cellranger-8.0.0/cellranger count ",
        "genome": {
            "Human": " --transcriptome=/igo/work/nabors/genomes/10X_Genomics/GEX/refdata-gex-GRCh38-2020-A ",
            "Mouse": " --transcriptome=/igo/work/nabors/genomes/10X_Genomics/GEX/refdata-gex-mm10-2020-A "
        }
    },
    "vdj": {
        "tool": " /igo/work/nabors/tools/cellranger-8.0.0/cellranger vdj ",
        "genome": {
            "Human": " --reference=/igo/work/genomes/10X_Genomics/VDJ/refdata-cellranger-vdj-GRCh38-alts-ensembl-7.0.0 ",
            "Mouse": " --reference=/igo/work/genomes/10X_Genomics/VDJ/refdata-cellranger-vdj-GRCm38-alts-ensembl-7.0.0 "
        }
    },
    "atac_count": {
        "tool": " /igo/work/nabors/tools/cellranger-atac-2.1.0/cellranger-atac count ",
        "genome": {
            "Human": " --reference=/igo/work/nabors/genomes/10X_Genomics/ATAC/refdata-cellranger-atac-GRCh38-1.0.1 ",
            "Mouse": " --reference=/igo/work/nabors/genomes/10X_Genomics/ATAC/refdata-cellranger-atac-mm10-1.1.0 "
        }
    },
    "cnv": {
        "tool": " /igo/work/nabors/tools/cellranger-dna-1.1.0/cellranger-dna cnv ",
        "genome": {
            "Human": " --reference=/igo/work/nabors/10X_Genomics/CNV/refdata-GRCh38-1.0.0 ",
            "Mouse": " --reference=/igo/work/nabors/10X_Genomics/CNV/refdata-GRCm38-1.0.0 "
        }
    },
    "multi": {
        "tool": " /igo/work/nabors/tools/cellranger-8.0.0/cellranger multi "
    },
    "arc": {
        "tool": " /igo/work/nabors/tools/cellranger-arc-2.0.2/cellranger-arc count ",
        "genome": {
            "Human": " --reference=/igo/work/nabors/genomes/10X_Genomics/ARC/refdata-cellranger-arc-GRCh38-2020-A-2.0.0 ",
            "Mouse": " --reference=/igo/work/nabors/genomes/10X_Genomics/ARC/refdata-cellranger-arc-mm10-2020-A-2.0.0 "
        }
    },
    "spaceranger": {
        "tool": " /igo/work/nabors/tools/spaceranger-3.0.0/spaceranger count ",
        "genome": {
            "Human": " --transcriptome=/igo/work/nabors/genomes/10X_Genomics/GEX/refdata-gex-GRCh38-2020-A ",
            "Mouse": " --transcriptome=/igo/work/nabors/genomes/10X_Genomics/spatial_gex/refdata-gex-mm10-2020-A "
        },
        "probe": {
            "Human": "/igo/work/nabors/genomes/10X_Genomics/spatial_gex/Visium_Human_Transcriptome_Probe_Set_v1.0_GRCh38-2020-A.csv",
            "Human_CytAssist": "/igo/work/genomes/10X_Genomics/spaceranger/Visium_Human_Transcriptome_Probe_Set_v2.0_GRCh38-2020-A.csv",
            "Mouse": "/igo/work/nabors/tools/spaceranger-3.0.0/external/tenx_feature_references/targeted_panels/Visium_Mouse_Transcriptome_Probe_Set_v1.0_mm10-2020-A.csv",
            "Mouse_HD": "/igo/work/nabors/tools/spaceranger-3.0.0/external/tenx_feature_references/targeted_panels/Visium_Mouse_Transcriptome_Probe_Set_v2.0_mm10-2020-A.csv"
        }
    }
}

# cellranger command line options
OPTIONS = " --create-bam=true --nopreflight --jobmode=lsf --mempercore=64 --disable-ui --maxjobs=200"

# 10X recipe list for different pipelines
COUNT_FLAVORS = ["10X_Genomics_GeneExpression-3", "10X_Genomics_GeneExpression-5"]
VDJ_FLAVORS = ["10X_Genomics_VDJ"]
ATAC_FLAVORS = ["10X_Genomics_ATAC"]
CNV_FLAVORS = ["10X_Genomics_CNV"]
ARC_FLAVORS = ["10X_Genomics_Multiome", "10X_Genomics_Multiome_ATAC", "10X_Genomics_Multiome_GeneExpression"]
SPATIAL_FLAVORS = ["10X_Genomics_Visium"]

# we do not want to PROCESS SAIL (15500) or SCRI (12437) projects
SCRI = "12437"
SAIL = "15500"
DO_NOT_PROCESS = [SCRI, SAIL]

VISIUM_ENDPOINT = "https://igolims.mskcc.org:8443/LimsRest/getConfig?igoId="
original_tiff_images_directory = "/rtssdc/mohibullahlab/IGO_Pipeline_Results/Single_Cell/10X_Genomics/TIFF_Images/"
tiff_images_directory = "/igo/work/igo/TIFF_Images/"
