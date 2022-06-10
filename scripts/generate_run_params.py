"""Determines run parameters to generate stats for input recipe & species
Args:
		recipe: Project recipe
		species: Project species
Usage: python generate_run_params.py -r HemeBrainPACT_v1 -s Human
"""

import re
import sys
import getopt
from collections import OrderedDict
from scripts.run_param_config import DEFAULT, GTAG, get_ordered_dic, recipe_type_mapping, species_genome_mapping, genome_reference_mapping, recipe_options_mapping, recipe_overrides

def find_mapping(mapping, target):
		"""Retrieves sample type from recipe

		Args:
			mapping (dic): regex-to-value mapping
			target (string): target match for regexes

		Returns:
			value (any): desired mapping of target
		"""
		for key_re, val in mapping.items():
				expr = re.compile(key_re, re.IGNORECASE)
				if(expr.match(target)):
						return val
		return None

def get_sample_type_from_recipe(recipe):
		"""Retrieves sample type from recipe

		Args:
			recipe: Recipe of the project

		Returns:
			sample_type_mapping, dic: Sample type of the project
			For Example:
				{ TYPE: "RNA" } , { TYPE: "DNA" }, { TYPE: "WGS" }
		"""
		return find_mapping(recipe_type_mapping, recipe)

def get_reference_configs(recipe, sample_type, species):
		"""Retrieves sample type from recipe

		Args:
			recipe: Project recipe
			sample_type: Project sample type (E.g. "RNA", "DNA")
			species: Project species
		Returns:
			dic: Contains project params w/ following potential keys,
				DEFAULT
				GENOME
				REFERENCE
				REF_FLAT
				RIBO_INTER
				GTAG

			For Example:
				{'GENOME': '/path/to/hg19.fa', 'REFERENCE': '/path/to/hg19/hg19.fa'}
		"""
		species_overrides = {} if species not in recipe_overrides else recipe_overrides[species]
		genome = find_mapping(species_overrides, recipe)
		if genome == None:
				genome = find_mapping(species_genome_mapping, species)
			
		mapping = find_mapping(genome_reference_mapping, genome)
	
	
		"""
		SHOULD TRANSFORM,
			type: RNA
			DEFAULT: {
					GENOME: "/igo/work/genomes/M.musculus/mm9/BWA_0.7.5a/mouse_mm9.fa",
					REFERENCE: "/igo/work/genomes/M.musculus/mm9/mouse_mm9.fa"
			},
			"RNA": {
					REF_FLAT: "/home/igo/resources/BED-Targets/mm9-Ref_Flat.txt",
					RIBO_INTER: "/home/igo/resources/BED-Targets/mm9.ribosomal.interval_file"
			}
		TO,
			{ 
					GENOME: "/igo/work/genomes/M.musculus/mm9/BWA_0.7.5a/mouse_mm9.fa",
					REFERENCE: "/igo/work/genomes/M.musculus/mm9/mouse_mm9.fa"
					REF_FLAT: "/home/igo/resources/BED-Targets/mm9-Ref_Flat.txt",
					RIBO_INTER: "/home/igo/resources/BED-Targets/mm9.ribosomal.interval_file"
			}
		"""
		genome_configs = mapping[DEFAULT].copy() # Base configuration for all recipes
		overrides = {} if sample_type not in mapping else mapping[sample_type]
		genome_configs.update(overrides)
		if 'GTAG' not in genome_configs:
				genome_configs.update( { GTAG: genome } )
		# GTAG = genome if 'GTAG' not in mapping else mapping[ 'GTAG' ]
		# genome_configs.update( { GTAG: GTAG } ) # Need to add the GTAG mapping
			
		return genome_configs

def get_recipe_options(recipe):
		"""Retrieves additional options for the given recipe

		Args:
			recipe: Project recipe
		Returns:
			dic: Contains recipe params w/ following potential keys,
				BAIT
				TARGET
				MSKQ
				MARKDUPLICATES

			For Example:
				{'BAIT': '/path/to/HemeBrainPACT_v1_BAITS.interval_list',
					'TARGET': '/path/to/HemeBrainPACT_v1_TARGETS.interval_list'}
		"""
		return find_mapping(recipe_options_mapping, recipe)


def main(argv):
		recipe = ''
		species = ''
		try:
				opts, args = getopt.getopt(argv,"hr:s:",["recipe=","species="])
		except getopt.GetoptError:
				print('usage: setup_stats.py -r <recipe> -s <species>')
				sys.exit(2)
		if(len(argv) == 0):
				print('usage: setup_stats.py -r <recipe> -s <species>')
				sys.exit(2)
		for opt, arg in opts:
				if opt == '-h':
						print('usage: setup_stats.py -r <recipe> -s <species>')
						sys.exit()
				elif opt in ("-r", "--recipe"):
						recipe = arg
				elif opt in ("-s", "--species"):
						species = arg
				else:
						print('usage: setup_stats.py -r <recipe> -s <species>')
						sys.exit(2)
					
		# Order is important because some assignments are dependent on previous assignments (see run_param_config.py)
		sample_type_dic = get_sample_type_from_recipe(recipe)
		sample_type = list(sample_type_dic.values())[0]     # { TYPE: "RNA" } -> "RNA"
		refr = get_reference_configs(recipe, sample_type, species).copy()
		opts = get_recipe_options(recipe)
		# TODO - Special Cases
		# "MethylCaptureSeq": sh $DIR/../PicardScripts/Methylseq.sh /igo/work/FASTQ/$RUNNAME/$PROJECT/
		# 10X_Genomics_*: sh $DIR/../PicardScripts/LaunchPipelines.sh $RUNTYPE --input /igo/work/FASTQ/$RUNNAME/$PROJECT/ --genome $GENOME --md $MARKDUPLICATES
		# DLP) echo "<<< DLP SAMPLES finished demuxing >>>"
	
		# Consolidate options
		refr.update(opts)
		refr.update(sample_type_dic)
	
		run_params = get_ordered_dic(refr) # Want to print in same order
		output=""
		for k,v in run_params.items():
				output="{} {}".format(output, "{}={}".format(k, v))
		# print(output)   # Required to output to stdout for nextflow script
		# return output
		return(run_params)

if __name__ == "__main__":
		main(sys.argv[1:])
	
	
