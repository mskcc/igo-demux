import glob
import sys
import re
import json
from platformdirs import os
import requests
import subprocess
from sys import path_importer_cache, stdout
from time import process_time

"""""
from getSampleManifests import get_sample_manifests
from getSampleManifests import get_igo_id_mappings
from getSampleManifests import get_igo_id
"""

sample_id_manifests = {
    # IGO_ID -> { "igoId": IGO_ID, "cmoPatientId": {CMO_PATIENT_ID}, "tumorOrNormal": {"Tumor"/"Normal"} }
    "05841_N_1007": {
        "igoId": "05841_N_1007",
        "sampleName": "05841_N_1007",
        "cmoPatientId": "",     # We want this to be blank because we don't want the same name as the "05841_N_7" sample
        "tumorOrNormal": "Tumor"
    },
    "06457_E_10013": {
        "igoId": "06457_E_10013",
        "sampleName": "06457_E_10013",
        "cmoPatientId": "",     # We want this to be blank because we don't want the same name as the real sample
        "tumorOrNormal": "Tumor"
    },
    "09641_Z_10040": {
        "igoId": "09641_Z_10040",
        "sampleName": "09641_Z_10040",
        "cmoPatientId": "",     # We want this to be blank because we don't want the same name as the real sample
        "tumorOrNormal": "Tumor"
    },
    "09641_Z_10056": {
        "igoId": "09641_Z_10056",
        "sampleName": "09641_Z_10056",
        "cmoPatientId": "",     # We want this to be blank because we don't want the same name as the real sample
        "tumorOrNormal": "Tumor"
    },
    "09641_Z_10083": {
        "igoId": "09641_Z_10083",
        "sampleName": "09641_Z_10083",
        "cmoPatientId": "",     # We want this to be blank because we don't want the same name as the real sample
        "tumorOrNormal": "Tumor"
    }
}

igoId_to_correctedCmoPatientId_map = {
    # { IGO_ID -> CORRECTED }
}

MAPPED_FIELDS = ["cmoPatientId", "tumorOrNormal"]

def fingerprint(project_id):
    project_id_p = "P"+project_id
    lims_host = 'igolims:8443'
    STATS_DIR = '/igo/staging/stats/'
    REFERENCE_SEQUENCE_DIR_38 = '/igo/work/genomes/H.sapiens/GRCh38.p13/GRCh38.p13.dna.primary.assembly.fa'
    REFERENCE_SEQUENCE_DIR_37 = '/igo/work/genomes/H.sapiens/hg19/human_hg19.fa'
    REFERENCE_DIR_ILLUMINA = '/igo/work/genomes/H.sapiens/hg38/hg38.fa' #Illumina Alt Aware graph reference file for WGS
    HAPLOTYPE_MAP_38 = '/home/igo/fingerprint_maps/map_files/hg38_igo.map'
    HAPLOTYPE_MAP_37 = '/home/igo/fingerprint_maps/map_files/hg19_nochr.map'
    t_start = process_time()

    vcfs = []
    
    #find all bams
    input_bams = set()
    print("Finding ___MD.bams for project {}".format(project_id_p))
    subprocess.call("cd /igo/staging/stats/", shell=True)
    for fileName in glob.glob('/igo/staging/stats/**/*___' + project_id_p + '___*___MD.bam', recursive=False):
        input_bams.add(fileName.split('/')[len(fileName.split('/')) - 1])
        print(fileName + " Added")

    if (len(input_bams) == 0):
        print('No Picard .bams found, searching again')
        for fileName in glob.glob('/igo/staging/stats/**/*___' + project_id_p + '___*.bam', recursive=False):
            print(fileName)
            baseName = os.path.basename(fileName)
            input_bams.add(baseName)
    
    print("Number of BAM files: ", len(input_bams))
    igo_ids=list(set(list(map(lambda bam: get_igo_id(bam), input_bams))))
    print("number of igo ids: " , len(igo_ids))
    sample_manifest = get_sample_manifests(igo_ids, lims_host)
    igo_id_mappings = get_igo_id_mappings(sample_manifest, MAPPED_FIELDS)

    processedIgoIds = []
    extractFingerprint_start = process_time()
    command0 = 'mkdir /igo/staging/stats/VCF/vcf_{}'.format(project_id)
    subprocess.call(command0, shell=True)
    HAPLOTYPE_MAP = ''
    REFERENCE_SEQUENCE = ''
    if 'h37___' in next(iter(input_bams)).lower():
        HAPLOTYPE_MAP = HAPLOTYPE_MAP_37
        REFERENCE_SEQUENCE = REFERENCE_SEQUENCE_DIR_37
    else:
        HAPLOTYPE_MAP = HAPLOTYPE_MAP_38
        REFERENCE_SEQUENCE = REFERENCE_SEQUENCE_DIR_38

    for bam in input_bams:
        regex = "IGO_([a-zA-Z0-9_]*?)___"
        igoId = re.findall(regex, bam)[0]
        patient_id = igo_id_mappings[igoId]['cmoPatientId']
        print("patient_id: " + patient_id)
        if (igoId not in processedIgoIds):
            processedIgoIds.append(igoId)
        else:
            print('Processed igo id, continuing to the next bam.')
            continue      

        EXECUTION_DIR = STATS_DIR + 'VCF/'
        output_vcf = EXECUTION_DIR + 'vcf_' + project_id + '/' + patient_id + '__' + project_id + '__' + igoId + '.vcf'
        runFolder = bam.split('___')[0]
        bam = STATS_DIR + runFolder + '/' + bam
        command1 = 'bsub -J "extract_fingerprint_{}" /home/igo/resources/gatk-4.1.9.0/gatk ExtractFingerprint --HAPLOTYPE_MAP \'{}\'  --INPUT \'{}\' --OUTPUT \'{}\' --REFERENCE_SEQUENCE \'{}\' --SAMPLE_ALIAS \'{}\''.format(igoId, HAPLOTYPE_MAP, bam, output_vcf, REFERENCE_SEQUENCE, patient_id)
        subprocess.call(command1, shell=True)
        print("Running extract fingerprint: " + command1)
        vcfs.append(output_vcf)

    extractFingerprint_finish = process_time()
    print('Elapsed time to extract fingerprint for all bams: ', extractFingerprint_finish - extractFingerprint_start)
    # TODO redirect bsub output to file and check return value of the command?
    command_crosscheck = 'bsub -n 9 -M 6 -w "ended(extract_fingerprint_*)" -J "CrosscheckFingerprint_{}" /home/igo/resources/gatk-4.1.9.0/gatk CrosscheckFingerprints LOD_THRESHOLD=-5.0 CROSSCHECK_BY=FILE NUM_THREADS=30 OUTPUT=/igo/staging/stats/VCF/vcf_{}/crosscheck_fingerprint_{}.tsv HAPLOTYPE_MAP=\'{}\' INPUT='.format(project_id, project_id, project_id, HAPLOTYPE_MAP)    
    vcfInputs = " INPUT=".join(vcfs)
    
    command_crosscheck += vcfInputs
    subprocess.call(command_crosscheck, shell=True)
    print("Running cross-check fingerprint: " + command_crosscheck)

    crosscheckFingerprint_stop = process_time()
    print('Elapsed time to CrosscheckFingerprints: ', crosscheckFingerprint_stop - extractFingerprint_finish)

    # Save to ngs-stats database
    done_path = '/igo/stats/DONE/crosscheck_metrics/{}'.format(project_id)
    if not os.path.exists(done_path):
        os.makedirs(done_path)
    # Example copy to file name:  /igo/stats/DONE/crosscheck_metrics/08236_J/08236_J.crosscheck_metrics
    copy_command = 'bsub -w "done(CrosscheckFingerprint_{})" cp /igo/staging/stats/VCF/vcdf_{}/crosscheck_fingerprint_{}.tsv /igo/stats/DONE/crosscheck_metrics/{}/{}.crosscheck_metrics'.format(project_id, project_id, project_id, project_id, project_id)
    subprocess.call(copy_command, shell=True)

    t_stop = process_time()
    print("Elapsed time to fingerprint in totall: ", t_stop - t_start) 
        

#################################
def get_project_id(file_name):
    """ Extracts project ID from intput BAM filename.

    :param file_name: string    e.g. "/PITT_0452_AHG2THBBXY_A1___P10344_C___13_cf_IGO_10344_C_20___hg19___MD.bam"
    :return: string             e.g. "P10344_C"
    """
    regex = "(?<=___)P[0-9]{5}[_A-Z,a-z]*(?=___)"  # Valid project ID is "P" + 5 numbers + (optional) [ "_" + 2 letters]
    matches = re.findall(regex, file_name)
    if len(matches) == 0:
        print("ERROR: Could not find Project ID in filename: %s with regex: \"%s\"" % (file_name, regex))
        sys.exit(1)
    if len(matches) > 1:
        print("WARNING: More than one match: %s" % str(matches))

    return matches[0]

def get_igo_id(file_name):
    """ Extracts IGO ID from intput BAM filename.

    :param file_name: string    e.g. "/PITT_0452_AHG2THBBXY_A1___P10344_C___13_cf_IGO_10344_C_20___hg19___MD.bam"
    :return: string             e.g. "10344_C_20"
    """
    regex = "IGO_([a-zA-Z0-9_]*?)___"
    matches = re.findall(regex, file_name)
    if len(matches) == 0:
        print("ERROR: Could not find IGO ID in filename: %s with regex: \"%s\"" % (file_name, regex))
        sys.exit(1)
    if len(matches) > 1:
        print("WARNING: More than one match: %s" % str(matches))

    return matches[0]


def get_sample_manifests(igo_ids, lims_host):
    """ Retrieves list of metadata (manifest) for the given IGO ID

    :param igo_ids: string[] List of target IGO IDs
    :return:
    """

    max_query = 10  # Maximum number of IDs that can be queried at one time
    idx = 0

    print("Retrieving sample manifests from %d igo_ids (max_query: %d)" % (len(igo_ids), max_query))

    sample_manifests = []
    while len(igo_ids) > idx:
        query_min = idx
        query_max = idx + max_query
        query_list = igo_ids[query_min:query_max]
        params = "igoSampleId=%s" % "&igoSampleId=".join(query_list)

        url = "https://%s/LimsRest/api/getSampleManifest?%s" % (lims_host, params)
        print("Retrieving manifests from samples [%d,%d): %s" % (query_min, query_max, ", ".join(query_list)))

        try:
            resp = requests.get(url, auth=('pms', 'tiagostarbuckslightbike'), verify=False)
        except:
            print("Request Failed: %s" % url)
            # We need to return a list of mappings, even if empty
            return [ { "igoId": id, "cmoPatientId": "NO_CMO_PID", "tumorOrNormal": "NA", "sampleName": id } for id in igo_ids  ]

        content = json.loads(resp.content)
        if (resp.status_code != 200 or len(content) == 0):
            # Warning
            print("getSampleManifest failed to return data for IGO IDs: %s. Service Response: %s" %
                  (str(query_list), resp.status_code))

        # returned order of manifests should be the same as @query_list
        for resp_idx in range(len(query_list)):
            igo_id = query_list[resp_idx]
            manifest = content[resp_idx]
            if is_invalid_manifest(manifest) and igo_id in sample_id_manifests:
                overridden_manifest = sample_id_manifests[igo_id]
                content[resp_idx] = overridden_manifest

        sample_manifests.extend(content)
        idx += max_query

    return sample_manifests

def is_invalid_manifest(sample_manifest):
    """ Checks that the sample_manifest has valid values for the required fields

    :param sample_manifests: Object     entry from the getSampleManifest API response
    """
    if not sample_manifest[ "igoId" ]:
        return True
    for field in MAPPED_FIELDS:
        if not sample_manifest[ field ]:
            return True

    return False

def get_igo_id_mappings(sample_manifests, mapped_fields):
    """ Returns mapping of igoIds to relevant information contained in sample_manifests

    :param sample_manifests: Object[]   List of sample manifests for IGO IDs
    :param mapped_fields: string[]      List of relevant fields for IGO ID
    :return: mapping Object
        e.g. {
                'IGO_ID_1':{
                    'cmoPatientId':'ID',
                    'tumorOrNormal':'Tumor'
                },
                'IGO_ID_2':{
                    'cmoPatientId':'ID_2',
                    'tumorOrNormal':'Normal'
                },
                ...
        }
    """
    dic = {}

    for manifest in sample_manifests:
        key_value = manifest["igoId"]

        if key_value in dic:
            print("Warning: overwriting manifest entry for %s" % key_value)

        entry = {}
        for field in mapped_fields:
            entry[field] = manifest[field]

        if key_value in igoId_to_correctedCmoPatientId_map:
            entry["cmoPatientId"] = igoId_to_correctedCmoPatientId_map[key_value]

        if entry["cmoPatientId"] == "":
            # TODO - A bit confusing. "mapped_fields" isn't a good idea - can't handle this situation well, hard to read
            entry["cmoPatientId"] = manifest["sampleName"]

        dic[key_value] = entry

    return dic       

# Runs fingerpritning for the project given on the command line
if __name__ == '__main__':
    if (sys.argv[1] == ''):
        print("Usage: fingerprinting_dag.py PROJECT_ID")
        exit(1)

    project_id = sys.argv[1]
    fingerprint(project_id)
