from collections import OrderedDict

# 1) Determined by recipe
TYPE = "TYPE"
# 2) Determined by genome & type (see: genome_reference_mapping)
DEFAULT = "DEFAULT_TYPE"
GENOME = "GENOME"
REFERENCE = "REFERENCE"
REF_FLAT = "REF_FLAT"
GTF = "GTF"
RIBOSOMAL_INTERVALS = "RIBOSOMAL_INTERVALS"
GTAG = "GTAG"
CELLRANGER_ATAC = "CELLRANGER_ATAC"
CELLRANGER_VDJ = "CELLRANGER_VDJ"
CELLRANGER_CNV = "CELLRANGER_CNV"
CELLRANGER_COUNT = "CELLRANGER_COUNT"
CELLRANGER_ARC = "CELLRANGER_ARC"
HAPLOTYPE_MAP = "HAPLOTYPE_MAP"

# 3) Determined by recipe (see: recipe_options_mapping)
BAITS="BAITS"
TARGETS="TARGETS"
MSKQ="MSKQ"
MD="MD"
DGN_REFERENCE="DGN_REFERENCE"
"""
				D E P E N D E N C Y    G R A P H
									+-----------+
			+----------->   opts    +-------------+
			|           +-----------+             |
			|                                     |
+-----+------+    +-----------+      +-----v------+
|   recipe   +---->   type    +------> run_params |
+-----+------+    +-----+-----+      +------^-----+
			|                 |                   |
			|           +-----v-----+             |
			+----------->   refr    +-------------+
									+-----+-----+             |
												^                   |
+------------+          |                   |
|   species  +----------+-------------------+
+------------+

species/refr -> recipe -> opts -> type          (Later options override previous options)

		species/refr, species_genome_mapping: Default reference for organism
		recipe, recipe_overrides: Specific reference required for recipe
		opts, recipe_options_mapping: More specific options for a recipe (Not reference)
		type, sample_type_mapping: WGS, RNA, or DNA
"""

def get_ordered_dic(unordered_dic):
		"""Returns a dictionary ordered by the length of its keys

		Args:
			unordered_dic: Dictionary of stuff to order

		Returns:
			type, OrderedDict: Ordered dictionary by key-length
		"""
		return OrderedDict(sorted(unordered_dic.items(), key=lambda t: -len(t[0])))

""" Mapping of recipes to their type, default should be DNA """
recipe_type_mapping_UNORDERED = {
		"MouseWholeGenome": { TYPE: "WGS" },
		"ShallowWGS": { TYPE: "WGS" },
		"10X_Genomics_WGS": { TYPE: "WGS" },
		"WholeGenomeSequencing": { TYPE: "WGS" },
		"HumanWholeGenome": { TYPE: "WGS" },
		".*RNA.*": { TYPE: "RNA" },
		".*96Well_SmartSeq2": { TYPE: "RNA" },
		".*SMARTer.*": { TYPE: "RNA" },
		"FusionDiscoverySeq": { TYPE: "RNA" },
		".*Ribo.*": { TYPE: "RNA" },
		".*CDH1_RNA.*": { TYPE: "CAPTURE" },
		# FOR NEW ENTRIES
		# "{regex}": { TYPE: type }
		".*": { TYPE: "DNA" }      # DEFAULT
}
recipe_type_mapping = get_ordered_dic(recipe_type_mapping_UNORDERED)

""" Recipes that should have determine recipe (instead of species -> genome logic) """
recipe_overrides = {
		"Human": {
				"ADCC1_v3": "GRCh37",
				"PCFDDR_.*": "hg19",
				"IWG.*": "hg19",
				"CH_v1": "hg19",
				"Twist_Exome": "hg19",
				"IDT_Exome_v1": "hg19",
				"ADCC1_v2": "hg19",
				"RDM": "hg19",
				"myTYPE_V1": "hg19",
				"PanCancerV2": "hg19",
				"MissionBio-Heme": "hg19",
				"WholeExome_v4": "hg19",
				"AmpliSeq": "hg19",
				"HemeBrainPACT_v1": "hg19"
		},
		"Mouse": {
				"M-IMPACT_v1": "mm10",
				"M-IMPACT_v2": "mm10",
				"Twist_mWES": "mm10",
				"10X_Genomics_Multiome": "mm10"
		}
}
""" Mapping of species to their genome-type """
species_genome_mapping_UNORDERED = {
		"Human": "GRCh38",
		"Mouse": "grcm39",
		"Mouse_GeneticallyModified": "grcm39",
		"Drosophilia": "dm3",
		"Zebrafish": "grcz11",
		"Chicken": "galGal4",
		".*uberculosis": "mtubf11",
		"S.Cerevisae": "sccer",
		"other": "sccer",
		"Fungus": "sccer",
		"E.Coli": "ecolik12",
		"Bacteria": "ecolik12",
		"C.Elegans": "ce10",
		"S.Pombe": "pombe",
		"R.norvegicus": "rn6",
		"E.Lambda": "elambda",
		"Plasmid": "ecolik12",
		"Pig": "Sscrofa11.1",
		# FOR NEW ENTRIES
		# "{regex}": "{GENOME}"
	
		".*": "INVALID"     # DEFAULT
}
species_genome_mapping = get_ordered_dic(species_genome_mapping_UNORDERED)

""" Mapping of genome-type to their reference files """
genome_reference_mapping_UNORDERED = {
		"hg19": {
				DEFAULT: {
						GENOME: "/igo/work/genomes/H.sapiens/hg19/BWA_0.7.5a/human_hg19.fa",
						REFERENCE: "/igo/work/genomes/H.sapiens/hg19/human_hg19.fa",
						HAPLOTYPE_MAP: "/home/igo/fingerprint_maps/map_files/hg19_ACCESS.map"
				},
				"RNA": {
						REF_FLAT: "/home/igo/resources/BED-Targets/hg19-Ref_Flat.txt",
						RIBOSOMAL_INTERVALS: "/igo/work//bed-targets/ucsc_hg19_rRNA.iList"
				},
				"miRNA": {
						RIBOSOMAL_INTERVALS: "/home/igo/resources/BED-Targets/ucsc_hg19_rRNA.iList",
						REF_FLAT: "/home/igo/resources/BED-Targets/hg19-Ref_Flat.txt",
						GENOME: "/home/igo/resources/BED-Targets/human_Counting_miRNA_Genome.fasta",
						GTAG: "human_miRNA"
				}
		},
		"grch37": {
				DEFAULT: {
						GENOME: "/igo/work/genomes/H.sapiens/GRCh37/GRCh37.fasta",
						REFERENCE: "/igo/work/genomes/H.sapiens/GRCh37/GRCh37.fasta",
						HAPLOTYPE_MAP: "/home/igo/fingerprint_maps/map_files/GRCh37_ACCESS.map"
				},
				"RNA": {
						GENOME: '/igo/work/nabors/bed_files/GRCh37_RNA_Ensembl/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa',
						REFERENCE: '/igo/work/nabors/bed_files/GRCh37_RNA_Ensembl/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa',
						REF_FLAT: '/igo/work/nabors/bed_files/GRCh37_RNA_Ensembl/refFlat_ensembl.v75.txt',
						RIBOSOMAL_INTERVALS: '/igo/work/nabors/bed_files/GRCh37_RNA_Ensembl/Homo_sapiens.GRCh37.75.rRNA.interval_list'
				},
		},
		"grch38": {
				DEFAULT: {
						GENOME: "/igo/work/nabors/genomes/hg38/hg38.fa",
						REFERENCE: "/igo/work/nabors/genomes/hg38/hg38.fa",
						CELLRANGER_ATAC: "/igo/work/nabors/genomes/10X_Genomics/ATAC/refdata-cellranger-atac-GRCh38-1.0.1",
						CELLRANGER_VDJ: "/igo/work/nabors/genomes/10X_Genomics/VDJ/refdata-cellranger-vdj-GRCh38-alts-ensembl-2.0.0",
						CELLRANGER_CNV: "/igo/work/nabors/10X_Genomics_references/CNV/refdata-GRCh38-1.0.0",
						CELLRANGER_COUNT: "/igo/work/nabors/genomes/10X_Genomics/GEX/refdata-gex-GRCh38-2020-A",
						CELLRANGER_ARC: "/igo/work/nabors/genomes/10X_Genomics/ARC/refdata-cellranger-arc-GRCh38-2020-A-2.0.0",
						HAPLOTYPE_MAP: "/home/igo/fingerprint_maps/map_files/hg38_igo.map"	# TODO - Verify this
				},
				"RNA": {
						GENOME: "/igo/work/nabors/genomes/hg38/hg38.fa",
						REFERENCE: "/igo/work/nabors/genomes/hg38/hg38.fa",
						REF_FLAT: "/igo/work/genomes/H.sapiens/hg38/ref_flat.txt",
						RIBOSOMAL_INTERVALS: "/igo/work/genomes/H.sapiens/hg38/ribosomal.interval_list",
						GTF: "/igo/work/genomes/H.sapiens/hg38/gencode.v40.chr_patch_hapl_scaff.annotation.gtf"
				},
		},
		"grcm38": {
				DEFAULT: {
						GENOME: '/igo/work/nabors/genomes/GRCm38/Mus_musculus.GRCm38.dna.primary_assembly.fa',
						REFERENCE: '/igo/work/nabors/genomes/GRCm38/Mus_musculus.GRCm38.dna.primary_assembly.fa',
						CELLRANGER_ATAC: "/igo/work/nabors/genomes/10X_Genomics/ATAC/refdata-cellranger-atac-mm10-1.1.0",
						CELLRANGER_VDJ: "/igo/work/nabors/genomes/10X_Genomics/VDJ/refdata-cellranger-vdj-GRCm38-alts-ensembl-2.2.0",
						CELLRANGER_CNV: "/igo/work/nabors/10X_Genomics_references/CNV/refdata-GRCm38-1.0.0",
						CELLRANGER_COUNT: "/igo/work/nabors/genomes/10X_Genomics/GEX/refdata-gex-mm10-2020-A"
				},
				"RNA": {
						REF_FLAT: "/igo/work/nabors/genomes/GRCm38/Mus_musculus.GRCm38.99.ref_flat",
						RIBOSOMAL_INTERVALS: "/igo/work/nabors/genomes/GRCm38/Mus_musculus.GRCm38.interval_list",
						GTF: "/igo/work/nabors/genomes/GRCm38/Mus_musculus.GRCm38.99.gtf"
				},
		},
		"grcm39": {
				DEFAULT: {
						GENOME: "/igo/work/nabors/genomes/grcm39/Mus_musculus.GRCm39.dna.primary_assembly.fa",
						REFERENCE: "/igo/work/nabors/genomes/grcm39/Mus_musculus.GRCm39.dna.primary_assembly.fa",
				},
				"RNA": {
						REF_FLAT: "/igo/work/nabors/genomes/grcm39/Mus_musculus.GRCm39.107.ref.flat",
						RIBOSOMAL_INTERVALS: "/igo/work/nabors/genomes/grcm39/Mus_musculus.GRCm39.107.rRNA.intervals",
						GTF: "/igo/work/nabors/genomes/grcm39/Mus_musculus.GRCm39.107.gtf"
				},
		},
		"mm9": {
				DEFAULT: {
						GENOME: "/igo/work/genomes/M.musculus/mm9/BWA_0.7.5a/mouse_mm9.fa",
						REFERENCE: "/igo/work/genomes/M.musculus/mm9/mouse_mm9.fa"
				},
				"RNA": {
						REF_FLAT: "/home/igo/resources/BED-Targets/mm9-Ref_Flat.txt",
						RIBOSOMAL_INTERVALS: "/home/igo/resources/BED-Targets/mm9.ribosomal.interval_file"
				},
				"miRNA": {
						REF_FLAT: "/home/igo/resources/BED-Targets/mm9-Ref_Flat.txt",
						GENOME: "/home/igo/resources/BED-Targets/mouse_Counting_miRNA_Genome.fasta",
						GTAG: "mouse_miRNA"
				}
		},
		"mm10": {
				DEFAULT: {
						GENOME: "/igo/work/genomes/M.musculus/mm10/BWA_0.7.5a/mouse_mm10__All.fa",
						REFERENCE: "/igo/work/genomes/M.musculus/mm10/BWA_0.7.5a/mouse_mm10__All.fa",
						CELLRANGER_ATAC: "/igo/work/nabors/genomes/10X_Genomics/ATAC/refdata-cellranger-atac-mm10-1.1.0",
						CELLRANGER_VDJ: "/igo/work/nabors/genomes/10X_Genomics/VDJ/refdata-cellranger-vdj-GRCm38-alts-ensembl-2.2.0",
						CELLRANGER_CNV: "/igo/work/nabors/10X_Genomics_references/CNV/refdata-GRCm38-1.0.0",
						CELLRANGER_COUNT: "/igo/work/nabors/genomes/10X_Genomics/GEX/refdata-gex-mm10-2020-A",
						CELLRANGER_ARC: "/igo/work/nabors/genomes/10X_Genomics/ARC/refdata-cellranger-arc-mm10-2020-A-2.0.0"
				},
				"RNA": {
						GENOME: '/igo/work/nabors/genomes/GRCm38/Mus_musculus.GRCm38.dna.primary_assembly.fa',
						REFERENCE: '/igo/work/nabors/genomes/GRCm38/Mus_musculus.GRCm38.dna.primary_assembly.fa',
						REF_FLAT: '/igo/work/nabors/genomes/GRCm38/Mus_musculus.GRCm38.99.ref_flat',
						RIBOSOMAL_INTERVALS: '/igo/work/nabors/genomes/GRCm38/Mus_musculus.GRCm38.interval_list',
						GTAG: 'GRCm38'
						# REF_FLAT: "/home/igo/resources/BED-Targets/mm10-Ref_Flat.txt",
						# RIBOSOMAL_INTERVALS: "/home/igo/resources/BED-Targets/mm10.ribosomal.interval_file"
				},
				"miRNA": {
						REF_FLAT: "/home/igo/resources/BED-Targets/mm10-Ref_Flat.txt",
						GENOME: "/home/igo/resources/BED-Targets/mouse_Counting_miRNA_Genome.fasta",
						GTAG: "mouse_miRNA"
				}
		},
		"Sscrofa11.1": {
				DEFAULT: {
						GENOME: "/igo/work/nabors/genomes/pig/ensembl/Sus_scrofa.Sscrofa11.1.dna_sm.toplevel.fa",
						REFERENCE: "/igo/work/nabors/genomes/pig/ensembl/Sus_scrofa.Sscrofa11.1.dna_sm.toplevel.fa",
				},
				"RNA": {
						GENOME: "/igo/work/nabors/genomes/pig/ensembl/Sus_scrofa.Sscrofa11.1.dna_sm.toplevel.fa",
						REFERENCE: "/igo/work/nabors/genomes/pig/ensembl/Sus_scrofa.Sscrofa11.1.dna_sm.toplevel.fa",
						REF_FLAT: "/igo/work/nabors/genomes/pig/ensembl/Sus_scrofa.Sscrofa11.1.109.ref.flat",
						RIBOSOMAL_INTERVALS: "/igo/work/nabors/genomes/pig/ensembl/Sus_scrofa.Sscrofa11.1.109.sorted.rRNA.intervals",
						GTF: "/igo/work/nabors/genomes/pig/ensembl/Sus_scrofa.Sscrofa11.1.109.gtf",
						GTAG: "Sscrofa11.1"
				}
		},
		"rn6": {
				DEFAULT: {
						GENOME: "/igo/work/genomes/R.norvegicus/Rnor_6.0/BWA_0.7.5a/Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa",
						REFERENCE: "/igo/work/genomes/R.norvegicus/Rnor_6.0/Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa"
				},
				"RNA": {
						REF_FLAT: "/home/igo/resources/BED-Targets/Rattus_norvegicus.Rnor_6.0.94.ref_flat",
						RIBOSOMAL_INTERVALS: "/home/igo/resources/BED-Targets/Rattus_norvegicus.Rnor_6.0.94.interval_list"
				}
		},
		"dm3": {
				DEFAULT: {
						GENOME: "/igo/work/genomes/D.melanogaster/dm3/BWA_0.7.5a/dm3.fa",
						REFERENCE: "/igo/work/genomes/D.melanogaster/dm3/BWA_0.7.5a/dm3.fa"
				},
				"RNA": {
						REF_FLAT: "/home/igo/resources/BED-Targets/dm3-Ref_Flat.txt"
				}
		},
		"galGal4": {
				DEFAULT: {
						GENOME: "/igo/work/genomes/G.gallus/BWA_0.7.5a/galGal4.fa",
						REFERENCE: "/igo/work/genomes/G.gallus/BWA_0.7.5a/galGal4.fa"
				},
				"RNA": {
						REF_FLAT: "/home/igo/resources/BED-Targets/galGal4-Ref_Flat.txt"
				}
		},
		"sccer": {
				DEFAULT: {
						GENOME: "/igo/work/genomes/S.cerevisiae/sacCer2/UCSC/BWA_0.7.5a/UCSC_sacCer2.fasta",
						REFERENCE: "/igo/work/genomes/S.cerevisiae/sacCer2/UCSC/BWA_0.7.5a/UCSC_sacCer2.fasta"
				}
		},
		"grcz11": {
				DEFAULT: {
						GENOME: "/igo/work/nabors/genomes/zebrafish/Danio_rerio.GRCz11.dna.primary_assembly.fa",
						REFERENCE: "/igo/work/nabors/genomes/zebrafish/Danio_rerio.GRCz11.dna.primary_assembly.fa"
				},
				"RNA": {
						GENOME: "/igo/work/nabors/genomes/zebrafish/Danio_rerio.GRCz11.dna.primary_assembly.fa",
						REFERENCE: "/igo/work/nabors/genomes/zebrafish/Danio_rerio.GRCz11.dna.primary_assembly.fa",
						REF_FLAT: "/igo/work/nabors/genomes/zebrafish/Danio_rerio.GRCz11.109.refFlat",
						RIBOSOMAL_INTERVALS: "/igo/work/nabors/genomes/zebrafish/Danio_rerio.GRCz11.109.rRNA.intervals",
						GTF: "/igo/work/nabors/genomes/zebrafish/Danio_rerio.GRCz11.109.gtf",
						GTAG: "grcz11"
				}
		},
		"ce10": {
				DEFAULT: {
						GENOME: "/igo/work/genomes/C.elegans/ce10/BWA_0.7.5a/ce10.fa",
						REFERENCE: "/igo/work/genomes/C.elegans/ce10/ce10.fa",
						REF_FLAT: "/home/igo/resources/BED-Targets/cel10-Ref_Flat.txt"
				}
		},
		"mtubf11": {
				DEFAULT: {
						GENOME: "/igo/work/genomes/M.tuberculosis/f11/BWA_0.7.5a/mycobacterium_tuberculosis_f11__finished__4_contigs.fasta",
						REFERENCE: "/igo/work/genomes/M.tuberculosis/f11/BWA_0.7.5a/mycobacterium_tuberculosis_f11__finished__4_contigs.fasta"
				},
				"RNA": {
						REF_FLAT: "/home/igo/resources/BED-Targets/mycobacterium_tuberculosis_f11__RefFlat.txt"
				}
		},
		"mtubh37r": {
				DEFAULT: {
						GENOME: "/igo/work/genomes/M.tuberculosis/h37rv_2/BWA_0.7.5a/mycobacterium_tuberculosis_h37rv2.fasta",
						REFERENCE: "/igo/work/genomes/M.tuberculosis/h37rv_2/BWA_0.7.5a/mycobacterium_tuberculosis_h37rv2.fasta"
				},
				"RNA": {
						REF_FLAT: "/home/igo/resources/BED-Targets/mycobacterium_tuberculosis_h37rv2_RefFlat.txt"
				}
		},
		"ecolik12": {
				DEFAULT: {
						GENOME: "/igo/work/genomes/E.coli/K12/MG1655/BWA_0.7.x/eColi__MG1655.fa",
						REFERENCE: "/igo/work/genomes/E.coli/K12/MG1655/BWA_0.7.x/eColi__MG1655.fa"
				}
		},
		"pseudomonas": {
				DEFAULT: {
						GENOME: "/igo/work/genomes/P.aeruginosa/BWA_0.7.5a/NC_002516.fa",
						REFERENCE: "/igo/work/genomes/P.aeruginosa/BWA_0.7.5a/NC_002516.fa"
				}
		},
		"pombe": {
				DEFAULT: {
						GENOME: "/igo/work/genomes/S.pombe/Ensembl/BWA_0.7.5a/Schizosaccharomyces_pombe.ASM294v2.20.fa",
						REFERENCE: "/igo/work/genomes/S.pombe/Ensembl/Schizosaccharomyces_pombe.ASM294v2.20.fa"
				}
		},
		"ercc": {
				DEFAULT: {
						GENOME: "/home/igo/resources/BED-Targets/ERCC/BWA_0.7.5a/ERCC.fasta"
				}
		},
		"elambda": {
				DEFAULT: {
						GENOME: "/igo/work/genomes/viruses/LAMBDA/Enterobacteriophage_lambda.fa",
						REFERENCE: "/igo/work/genomes/viruses/LAMBDA/Enterobacteriophage_lambda.fa"
				}
		},
		"ct24": {
				DEFAULT: {
						GENOME: '/igo/work/nabors/C_thermophilum/C_thermophilum.annotation.v2.4/C_thermophilum.annotation.v2.4.scaff.fa',
						REFERENCE: '/igo/work/nabors/C_thermophilum/C_thermophilum.annotation.v2.4/C_thermophilum.annotation.v2.4.scaff.fa'
				},
				"RNA": {
						REF_FLAT: '/igo/work/nabors/C_thermophilum/C_thermophilum.annotation.v2.4/C_thermophilum.annotation.v2.4.ref_flat',
						RIBOSOMAL_INTERVALS: '/igo/work/nabors/C_thermophilum/C_thermophilum.annotation.v2.4/C_thermophilum.annotation.v2.4.interval_list'
				}
		},
		".*": {
				DEFAULT: {
						GENOME: "INVALID",
						REFERENCE: "INVALID"
				}
		},
}
genome_reference_mapping = get_ordered_dic(genome_reference_mapping_UNORDERED)

""" Mapping of recipe to additional run options """
recipe_options_mapping_UNORDERED = {
		"IMPACT341": {
				BAITS: "/home/igo/resources/ilist/IMPACT341/b37/picard_baits.interval_list",
				TARGETS: "/home/igo/resources/ilist/IMPACT341/b37/picard_targets.interval_list",
				MSKQ: "yes",
				MD: "yes"
		},
		"IMPACT410.*": {
				BAITS: "/home/igo/resources/ilist/IMPACT410/b37/picard_baits.interval_list",
				TARGETS: "/home/igo/resources/ilist/IMPACT410/b37/picard_targets.interval_list",
				MSKQ: "yes",
				MD: "yes"
		},
		"IMPACT468": {
				# NOTE: interval list file name "IMPACT468_BAITS" is stored in LIMS and passed to pipelines, change file name with caution
				BAITS: "/home/igo/resources/ilist/IMPACT468/b37/IMPACT468_BAITS.interval_list",
				TARGETS: "/home/igo/resources/ilist/IMPACT468/b37/IMPACT468_TARGETS.interval_list",
				MSKQ: "yes",
				MD: "yes"
		},
		"IMPACT505": {
				# NOTE: interval list file name "IMPACT468_BAITS" is stored in LIMS and passed to pipelines, change file name with caution
				BAITS: "/igo/home/igo/resources/ilist/hg38/IMPACT505/IMPACT505_BAITS.baits",
				TARGETS: "/igo/home/igo/resources/ilist/hg38/IMPACT505/IMPACT505_TARGETS.targets",
				MSKQ: "yes",
				MD: "yes"
		},
		"HemePACT_v4": {
				BAITS: "/igo/home/igo/resources/ilist/hg38/HemePACT_v4/HemePACT_v4_BAITS.baits",
				TARGETS: "/igo/home/igo/resources/ilist/hg38/HemePACT_v4/HemePACT_v4_TARGETS.targets",
				MSKQ: "yes",
				MD: "yes"
		},
		"IMPACT_Heme_v2": {
				BAITS: "/igo/home/igo/resources/ilist/hg38/IMPACT-Heme_v2/IMPACT-Heme_v2_BAITS.iList",
				TARGETS: "/igo/home/igo/resources/ilist/hg38/IMPACT-Heme_v2/IMPACT-Heme_v2_TARGETS.iList",
				MSKQ: "yes",
				MD: "yes"
		},
		"M-IMPACT_v1": {
				BAITS: "/home/igo/resources/BED-Targets/IMPACT/MM_IMPACT/mm_IMPACT_v1_mm10_BAITS.iList",
				TARGETS: "/home/igo/resources/BED-Targets/IMPACT/MM_IMPACT/mm_IMPACT_v1_mm10_TARGETS.iList",
				MSKQ: "yes",
				MD: "yes"
		},
		"M-IMPACT_v2": {
				BAITS: "/home/igo/resources/BED-Targets/IMPACT/MM_IMPACT/M-IMPACT_v2.baits",
				TARGETS: "/home/igo/resources/BED-Targets/IMPACT/MM_IMPACT/M-IMPACT_v2.targets",
				MSKQ: "yes",
				MD: "yes"
		},
		"WholeExomeSequencing": {
				BAITS: "/home/igo/resources/ilist/IDT_Exome_v1_FP/b37/IDT_Exome_v1_FP_b37_baits.interval_list",
				TARGETS: "/home/igo/resources/ilist/IDT_Exome_v1_FP/b37/IDT_Exome_v1_FP_b37_targets.interval_list",
				MSKQ: "no",
				MD: "yes"
		},
		"Twist_Exome": {
				# TODO - Delete "Twist_Exome" or change interval lists to be GRCh37
				BAITS: "/home/igo/resources/BED-Targets/Twist/Twist_Exome_Hg19_TARGETS.iList",
				TARGETS: "/home/igo/resources/BED-Targets/Twist/Twist_Exome_Hg19_TARGETS.iList",
				MSKQ: "no",
				MD: "yes"
		},
		"Agilent_v4_51MB_Human": {
				BAITS: "/home/igo/resources/ilist/AgilentExon_51MB_b37_v3/b37/AgilentExon_51MB_b37_v3_baits.interval_list",
				TARGETS: "/home/igo/resources/ilist/AgilentExon_51MB_b37_v3/b37/AgilentExon_51MB_b37_v3_targets.interval_list",
				MSKQ: "no",
				MD: "yes"
		},
		"Agilent_MouseAllExonV1": {
				BAITS: "/home/igo/resources/BED-Targets/Agilent_MouseAllExonV1_mm10_v1_baits.ilist",
				TARGETS: "/home/igo/resources/BED-Targets/Agilent_MouseAllExonV1_mm10_v1_targets.ilist",
				MSKQ: "no",
				MD: "yes"
		},
		"IDT_Exome_v1_FP_Viral_Probes": {
				BAITS: "/home/igo/resources/ilist/IDT_Exome_v1_FP/b37/IDT_Exome_v1_FP_b37_baits.interval_list",
				TARGETS: "/home/igo/resources/ilist/IDT_Exome_v1_FP/b37/IDT_Exome_v1_FP_b37_targets.interval_list",
				MSKQ: "no",
				MD: "yes"
		},
		"IDT_Exome_v2_FP_Viral_Probes": {
			BAITS: "/igo/home/igo/resources/ilist/hg38/IDT_Exome_v2_FP/IDT_Exome_v2_FP_BAITS.baits",
			TARGETS: "/igo/home/igo/resources/ilist/hg38/IDT_Exome_v2_FP/IDT_Exome_v2_FP_TARGETS.targets",
			MSKQ: "no",
			MD: "yes"
		},
		"IDT_Exome_v1": {
				BAITS: "/home/igo/resources/BED-Targets/xgen-exome-research-panel-BAITS.iList",
				TARGETS: "/home/igo/resources/BED-Targets/xgen-exome-research-panel-TARGETS.iList",
				MSKQ: "no",
				MD: "yes"
		},
		"IDT_Exome_v1_FP": {
				BAITS: "/home/igo/resources/ilist/IDT_Exome_v1_FP/b37/IDT_Exome_v1_FP_b37_baits.interval_list",
				TARGETS: "/home/igo/resources/ilist/IDT_Exome_v1_FP/b37/IDT_Exome_v1_FP_b37_targets.interval_list",
				MSKQ: "no",
				MD: "yes"
		},
		"ADCC1_v2": {
			BAITS: "/home/igo/resources/BED-Targets/ADCC1_V2/ADCC1_V2_BAITS.iList",
				TARGETS: "/home/igo/resources/BED-Targets/ADCC1_V2/ADCC1_V2_TARGETS.iList",
				MSKQ: "no",
				MD: "yes"
		},
		"ADCC1_v3": {
				BAITS: "/igo/work/nabors/bed_files/ADCC1_v3/ADCC1_v3_capture_targets_BAITS.interval_list",
				TARGETS: "/igo/work/nabors/bed_files/ADCC1_v3/ADCC1_v3_primary_targets_TARGETS.interval_list",
				MSKQ: "no",
				MD: "yes"
		},
		"RDM": {
				BAITS: "/home/igo/resources/BED-Targets/Rdm_Final_BAIT.iList",
				TARGETS: "/home/igo/resources/BED-Targets/Rdm_Final_TARGET.iList",
				MSKQ: "yes",
				MD: "yes"
		},
		"myTYPE_V1": {
				BAITS: "/home/igo/resources/BED-Targets/MM_MSK_permiss_BAIT.iList",
				TARGETS: "/home/igo/resources/BED-Targets/MM_MSK_permiss_TARGET.iList",
				MSKQ: "no",
				MD: "yes"
		},
		"NR3C1": {
				BAITS: "/home/igo/resources/BED-Targets/mm10_NR3C1_BAITS.iList",
				TARGETS: "/home/igo/resources/BED-Targets/mm10_NR3C1_TARGETS.iList",
				MSKQ: "yes",
				MD: "yes"
		},
		"MSK-ACCESS_v1": {
				BAITS: "/igo/home/igo/resources/ilist/hg38/MSK-ACCESS-v1/MSK-ACCESS-v1_0-probesAllwFP.baits",
				TARGETS: "/igo/home/igo/resources/ilist/hg38/MSK-ACCESS-v1/MSK-ACCESS-v1_0-probesAllwFP.targets",
				MSKQ: "no",
				MD: "yes",
				HAPLOTYPE_MAP: "/home/igo/fingerprint_maps/map_files/hg38_no_chr_ACCESS_unordered.map"
		},
		"PanCancerV2": {
				BAITS: "/home/igo/resources/BED-Targets/PanCancerV2/PanCancerV2_BAITS.iList",
				TARGETS: "/home/igo/resources/BED-Targets/PanCancerV2/PanCancerV2_TARGETS.iList",
				MSKQ: "no",
				MD: "yes"
		},
		"CH_v1": {
				BAITS: "/home/igo/resources/BED-Targets/CH_v1/CH_v1_BAITS.interval_list",
				TARGETS: "/home/igo/resources/BED-Targets/CH_v1/CH_v1_TARGETS.interval_list",
				MSKQ: "no",
				MD: "yes"
		},
		"MissionBio-Heme": {
				BAITS: "/home/igo/resources/BED-Targets/Mission_Bio/AML_BAITS.interval_list",
				TARGETS: "/home/igo/resources/BED-Targets/Mission_Bio/AML_TARGETS.interval_list",
				MSKQ: "no",
				MD: "yes"
		},
		"WholeExome_v4": {
				BAITS: "/home/igo/resources/BED-Targets/IDT_Exome_v1_FP_BAITS.iList",
				TARGETS: "/home/igo/resources/BED-Targets/IDT_Exome_v1_FP_TARGETS.iList",
				MSKQ: "no",
				MD: "yes"
		},
		"AmpliSeq": {
				BAITS: "/home/igo/resources/BED-Targets/AmpliSeq.ComprehensiveCancerPanel/ComprehensiveCancer.dna_manifest.20180509.BAITS.interval_list",
				TARGETS: "/home/igo/resources/BED-Targets/AmpliSeq.ComprehensiveCancerPanel/ComprehensiveCancer.dna_manifest.20180509.TARGETS.interval_list",
				MSKQ: "no",
				MD: "yes"
		},
		"PCFDDR_.*": {
				BAITS: "/home/igo/resources/BED-Targets/PCFDDR_v1/PCFDDR_v1__BAITS.interval_list",
				TARGETS: "/home/igo/resources/BED-Targets/PCFDDR_v1/PCFDDR_v1__TARGETS.interval_list",
				MSKQ: "no",
				MD: "yes"
		},
		"HemeBrainPACT_v1": {
				BAITS: "/home/igo/resources/BED-Targets/HemeBrainPACT_v1/HemeBrainPACT_v1_BAITS.interval_list",
				TARGETS: "/home/igo/resources/BED-Targets/HemeBrainPACT_v1/HemeBrainPACT_v1_TARGETS.interval_list",
				MSKQ: "no",
				MD: "yes"
		},
		"Twist_mWES": {
			BAITS: "/home/igo/resources/ilist/Twist_mWES/Twist_mWES_BAITS.IntervalList",
			TARGETS: "/home/igo/resources/ilist/Twist_mWES/Twist_mWES_TARGETS.IntervalList"
		},
		"Twist_Kentsis_spikeinWES_RK_V3": {
			BAITS: "/home/igo/resources/ilist/Twist_Kentsis_spikeinWES_RK_V3/Twist_Kentsis_spikeinWES_RK_V3_BAITS.intervalList",
			TARGETS: "/home/igo/resources/ilist/Twist_Kentsis_spikeinWES_RK_V3/Twist_Kentsis_spikeinWES_RK_V3_TARGETS.intervalList"
		},
		"UBA1_plus_v1": {
				BAITS: "/home/igo/resources/BED-Targets/UBA1_plus_v1/UBA1_plus_v1.baits",
				TARGETS: "/home/igo/resources/BED-Targets/UBA1_plus_v1/UBA1_plus_v1.targets",
				MSKQ: "no",
				MD: "yes"
		},
		"diaz_CH_v3": {
				BAITS: "/igo/work/nabors/bed_files/diaz_CH_v3/diaz_CH_v3.baits",
				TARGETS: "/igo/work/nabors/bed_files/diaz_CH_v3/diaz_CH_v3.targets",
				MSKQ: "no",
				MD: "yes"
		},
		"CDH1_RNA": {
				BAITS: "/igo/work/nabors/bed_files/CDH1_RNA/CDH1/CDH1.hg38.baits",
				TARGETS: "/igo/work/nabors/bed_files/CDH1_RNA/CDH1/CDH1.hg38.targets",
				MSKQ: "no",
				MD: "yes"
		},
		"08390_Hg19": {
				BAITS: "/igo/work/nabors/bed_files/IMPACT505_PLUS_for_08390/impact505_and_MC1R/hg38/08390-IMPACT505-MC1R-BAITS.hg38.iList",
				TARGETS: "/igo/work/nabors/bed_files/IMPACT505_PLUS_for_08390/impact505_and_MC1R/hg38/08390-IMPACT505-MC1R-TARGETS.hg38.iList",
				MSKQ: "no",
				MD: "yes"
		},
		"CD274_WJ_v1": {
				BAITS: "/igo/work/nabors/bed_files/CD274_WJ_v1/hg38/CD274_WJ_v1_hg38.baits",
				TARGETS: "/igo/work/nabors/bed_files/CD274_WJ_v1/hg38/CD274_WJ_v1_hg38.targets",
				MSKQ: "no",
				MD: "yes"
		},
		"MethylCaptureSeq": {
				BAITS: "/igo/home/igo/resources/ilist/hg38/MethylCaptureSeq/truseq-methyl-capture-epic-manifest-file-hg38.baits.ilist",
				TARGETS: "/igo/home/igo/resources/ilist/hg38/MethylCaptureSeq/truseq-methyl-capture-epic-manifest-file-hg38.targets.ilist",
				MSKQ: "no",
				MD: "yes"
		},
		"HumanWholeGenome": {
				MSKQ: "no",
				MD: "yes",
				HAPLOTYPE_MAP: "", # TODO - Add this
				GENOME: "/igo/work/genomes/H.sapiens/GRCh38.p13/ncbi-genomes-2021-09-23/GCF_000001405.39_GRCh38.p13_genomic.fna",  # References that created DRAGEN reference
				REFERENCE: "/igo/work/genomes/H.sapiens/GRCh38.p13/ncbi-genomes-2021-09-23/GCF_000001405.39_GRCh38.p13_genomic.fna",
				DGN_REFERENCE: "/staging/ref/hg38_alt_masked_graph_v2+cnv+graph+rna-8-1644018559"
		},
		"MouseWholeGenome": {
				MSKQ: "no",
				MD: "yes"
				# TODO
				# sh $DIR/../PicardScripts/LaunchPipelines.sh $RUNTYPE --input /igo/work/FASTQ/$RUNNAME/$PROJECT/ --genome $GENOME --type WGS --md $MARKDUPLICATES --mskq $MSKQ
		},
		"ShallowWGS": {
				MSKQ: "no",
				MD: "yes"
				# TODO
				# sh $DIR/../PicardScripts/LaunchPipelines.sh $RUNTYPE --input /igo/work/FASTQ/$RUNNAME/$PROJECT/ --genome $GENOME --type WGS --md $MARKDUPLICATES --mskq $MSKQ
		},
		"10X_Genomics_WGS": {
				MSKQ: "no",
				MD: "yes"
				# TODO
				# sh $DIR/../PicardScripts/LaunchPipelines.sh $RUNTYPE --input /igo/work/FASTQ/$RUNNAME/$PROJECT/ --genome $GENOME --type WGS --md $MARKDUPLICATES --mskq $MSKQ
		},
		"CustomAmplificationPCR": {
				MSKQ: "no",
				MD: "no"
				# TODO
				# sh $DIR/../PicardScripts/LaunchPipelines.sh $RUNTYPE --input /igo/work/FASTQ/$RUNNAME/$PROJECT/ --genome $GENOME --md $MARKDUPLICATES --mskq $MSKQ
		},
		"AmpliconSeq": {
				MSKQ: "no",
				MD: "yes"
				# TODO
				# sh $DIR/../PicardScripts/LaunchPipelines.sh $RUNTYPE --input /igo/work/FASTQ/$RUNNAME/$PROJECT/ --genome $GENOME --md $MARKDUPLICATES --mskq $MSKQ
		},
		"CRISPRSeq": {
				MSKQ: "no",
				MD: "yes"
				# TODO
				# sh ~/Scripts/PicardScripts/LaunchPipelines.sh $RUNTYPE --input /igo/work/FASTQ/$RUNNAME/$PROJECT/ --genome $GENOME --md $MARKDUPLICATES --mskq $MSKQ
		},
		"FusionDiscoverySeq": {
				MSKQ: "no",
				MD: "yes"
				# TODO
				# sh $DIR/../PicardScripts/LaunchPipelines.sh $RUNTYPE --input /igo/work/FASTQ/$RUNNAME/$PROJECT/ --genome $GENOME --type RNA --md $MARKDUPLICATES --mskq $MSKQ
		},
		".*RNA.*": {
				MSKQ: "no",
				MD: "yes"
				# TODO
				# sh $DIR/../PicardScripts/LaunchPipelines.sh $RUNTYPE --input /igo/work/FASTQ/$RUNNAME/$PROJECT/ --genome $GENOME --type RNA --md $MARKDUPLICATES --mskq $MSKQ
		},
		".*Ribo.*": {
				MSKQ: "no",
				MD: "yes"
				# TODO
				# sh $DIR/../PicardScripts/LaunchPipelines.sh $RUNTYPE --input /igo/work/FASTQ/$RUNNAME/$PROJECT/ --genome $GENOME --type RNA --md $MARKDUPLICATES --mskq $MSKQ
		},
		".*SMARTer.*": {
				MSKQ: "no",
				MD: "yes"
				# TODO
				# sh $DIR/../PicardScripts/LaunchPipelines.sh $RUNTYPE --input /igo/work/FASTQ/$RUNNAME/$PROJECT/ --genome $GENOME --type RNA --md $MARKDUPLICATES --mskq $MSKQ
		},
		".*96Well_SmartSeq2.*": {
				MSKQ: "no",
				MD: "yes"
				# TODO
				# sh $DIR/../PicardScripts/LaunchPipelines.sh $RUNTYPE --input /igo/work/FASTQ/$RUNNAME/$PROJECT/ --genome $GENOME --type RNA --md $MARKDUPLICATES --mskq $MSKQ
		},
		".*HemePACT.*v3.*": {
				BAITS: "/home/igo/resources/ilist/HemePACT_v3/b37/HemePACT_v3_b37_baits.ilist",
				TARGETS: "/home/igo/resources/ilist/HemePACT_v3/b37/HemePACT_v3_b37_targets.ilist",
				MSKQ: "yes",
				MD: "yes"
		},
			"IWG.*": {
				BAITS: "/home/igo/resources/BED-Targets/papaemme_IWG_OID43089_hg19_MHC_RNA_max10_20oct2015_BAITS.iList",
				TARGETS: "/home/igo/resources/BED-Targets/papaemme_IWG_OID43089_hg19_MHC_RNA_max10_20oct2015_TARGETS.iList",
				MSKQ: "no",
				MD: "yes"
		},
		".*08428_mm10.*": {
				BAITS: "/home/igo/resources/BED-Targets/8428_mm10/8428_mm10_BAITS.iList",
				TARGETS: "/home/igo/resources/BED-Targets/8428_mm10/8428_mm10_TARGETS.iList",
				MSKQ: "no",
				MD: "yes"
		},
		# FOR NEW ENTRIES
		# "{regex}": { [KEY]: "{REFERENCE_FILE}", ... }
	
		".*": {
				MSKQ: "no",
				MD: "yes"
		}     # DEFAULT
}
recipe_options_mapping = get_ordered_dic(recipe_options_mapping_UNORDERED)
