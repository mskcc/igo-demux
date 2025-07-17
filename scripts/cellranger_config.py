# work folder
STATS_AREA = "/igo/staging/CELLRANGER/"

# config info 
ACCESS = 0o775
config_dict = {
    "count": {
        "tool": " /igo/work/nabors/tools/cellranger-9.0.1/cellranger count ",
        "genome": {
            "Human": " --transcriptome=/igo/work/nabors/genomes/10X_Genomics/GEX/refdata-gex-GRCh38-2020-A ",
            "Mouse": " --transcriptome=/igo/work/nabors/genomes/10X_Genomics/GEX/refdata-gex-mm10-2020-A "
        }
    },
    "vdj": {
        "tool": " /igo/work/nabors/tools/cellranger-9.0.1/cellranger vdj ",
        "genome": {
            "Human": " --reference=/igo/work/genomes/10X_Genomics/VDJ/refdata-cellranger-vdj-GRCh38-alts-ensembl-7.0.0 ",
            "Mouse": " --reference=/igo/work/genomes/10X_Genomics/VDJ/refdata-cellranger-vdj-GRCm38-alts-ensembl-7.0.0 "
        }
    },
    "multi": {
        "tool": " /igo/work/nabors/tools/cellranger-9.0.1/cellranger multi "
    },
    "ocm": {
        "tool": " /igo/work/nabors/tools/cellranger-9.0.1/cellranger multi "
    },
    "arc": {
        "tool": " /igo/work/nabors/tools/cellranger-arc-2.0.2/cellranger-arc count ",
        "genome": {
            "Human": " --reference=/igo/work/nabors/genomes/10X_Genomics/ARC/refdata-cellranger-arc-GRCh38-2020-A-2.0.0 ",
            "Mouse": " --reference=/igo/work/nabors/genomes/10X_Genomics/ARC/refdata-cellranger-arc-mm10-2020-A-2.0.0 "
        }
    },
    "spaceranger": {
        "tool": " /igo/work/nabors/tools/spaceranger-4.0.1/spaceranger count ",
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
ARC_OPTIONS = " --nopreflight --jobmode=lsf --mempercore=64 --disable-ui --maxjobs=200"

# 10X recipe list for different pipelines
TAG_DICT = {
    "SC_Chromium-GEX-3": "count",
    "SC_Chromium-GEX-5": "count",
    "SC_Chromium-TCR": "vdj_t",
    "SC_Chromium-BCR": "vdj_b",
    "SC_Chromium-Multiome": "arc", 
    "SC_Chromium-Multiome-ATAC": "arc", 
    "SC_Chromium-Multiome-GEX": "arc",
    "ST_Visium": "spaceranger",
    "ST_Visium-HD": "spaceranger",
    "SC_Chromium-OCM": "ocm",
    "SC_Chromium-Flex": "flex"
}

# we do not want to PROCESS SAIL (15500) or SCRI (12437) projects, adding new chemistry for SAIL, 16364
SCRI = "12437"
SAIL = "15500"
DO_NOT_PROCESS = [SCRI, SAIL, "16364"]

VISIUM_ENDPOINT = "https://igolims.mskcc.org:8443/LimsRest/getConfig?igoId="
original_tiff_images_directory = "/rtssdc/mohibullahlab/IGO_Pipeline_Results/Single_Cell/10X_Genomics/TIFF_Images/"
tiff_images_directory = "/igo/work/igo/TIFF_Images/"
CONFIG_AREA = "/igo/stats/Multi_config/"
DRIVE_LOCATION = "/igo/work/igo/Cellranger_Multi_Config/"
ORIGIN_DRIVE_LOCATION = "/rtssdc/mohibullahlab/LIMS/LIMS_cellranger_multi/"
