import glob
import sys
import re
import json
import requests
import subprocess
from sys import path_importer_cache, stdout
from time import process_time
from rich import print_json

from sqlalchemy import true
import config
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
    lims_host = 'igo-lims04.mskcc.org:8443'
    STATS_DIR = '/igo/staging/stats/'
    REFERENCE_SEQUENCE_DIR = '/igo/work/genomes/H.sapiens/GRCh38.p13/GRCh38.p13.dna.primary.assembly.fa'
    HAPLOTYPE_MAP = '/home/igo/fingerprint_maps/map_files/hg38_igo.map'
    t_start = process_time()

    vcfs = []
    
    #find all bams
    input_bams = set()
    print("Finding bams of the run argument...")
    subprocess.call("cd /igo/staging/stats/", shell=True)
    for fileName in glob.glob('./**/*___' + project_id + '___*___MD.bam', recursive=True):
        input_bams.add(fileName.split('/')[len(fileName.split('/') - 1)])
        print(fileName + " Added")

    
    if(len(input_bams) == 0):
        REFERENCE_SEQUENCE_DIR = '/igo/work//genomes/H.sapiens/hg38/hg38.fa' # TODO: illumina ref dir
        for fileName in glob.glob('./**/*___' + project_id + '___*___.bam', recursive=True):
            input_bams.add(fileName.split('/')[len(fileName.split('/') - 1)])
            print(fileName + " Added")
            
        
    # else:
    #     for bam in input_bams:
    #         bam = bam.split('IGO')[1] # review this logic
    print("Number of BAM files: ", len(input_bams))
    igo_ids=list(set(list(map(lambda bam: get_igo_id(bam), input_bams))))
    sample_manifest = get_sample_manifests(igo_ids, lims_host)
    igo_id_mappings = get_igo_id_mappings(sample_manifest, MAPPED_FIELDS)
    
    for bam in input_bams:
        regex = "IGO_([a-zA-Z0-9_.-]*?)___"
        igoId = re.findall(regex, bam)

    i = 0
    processedIgoIds = []
    for bam in sorted(input_bams):
        patient_id = igo_id_mappings[igo_ids[i]]['cmoPatientId']
        regex = "IGO_([a-zA-Z0-9_.-]*?)___"
        igoId = re.findall(regex, bam)
        if(igoId[0] not in processedIgoIds):
            processedIgoIds.append(igoId[0])
        else:
            print('Processed igo id, breaking out of loop.')
            continue
        print("patient_id: " + patient_id)
        i += 1 

        EXECUTION_DIR = STATS_DIR
        output_vcf = EXECUTION_DIR + 'VCF/' + patient_id + '_' + project_id + '_' + igoId[0] + '_.vcf'
        runFolder = bam.split('___')[0]
        bam = EXECUTION_DIR + runFolder + '/' + bam

        command1 = '/home/igo/resources/gatk-4.1.9.0/gatk ExtractFingerprint --HAPLOTYPE_MAP \'{}\'  --INPUT \'{}\' --OUTPUT \'{}\' --REFERENCE_SEQUENCE \'{}\' --SAMPLE_ALIAS \'{}\''.format(HAPLOTYPE_MAP, bam, output_vcf, REFERENCE_SEQUENCE_DIR, patient_id)
        print("Running extract fingerprint: " + command1)
        subprocess.call(command1, shell=True)
        vcfs.append(output_vcf)


    command2 = '/home/igo/resources/gatk-4.1.9.0/gatk CrosscheckFingerprints LOD_THRESHOLD=-5.0 CROSSCHECK_BY=FILE NUM_THREADS=30 OUTPUT=crosscheck_fingerprint.tsv HAPLOTYPE_MAP=\'{}\' INPUT='.format(HAPLOTYPE_MAP)
    listOfInputs = []
    for vcf in vcfs:
        listOfInputs.append(vcf)
        
    vcfInputs = " INPUT=".join(listOfInputs)
    command2 += vcfInputs
    print("Running cross-check fingerprint: " + command2)
    subprocess.call(command2, shell=True)

    t_stop = process_time()
    print("Elapsed time to fingerprint: ", t_stop - t_start)
     
        

#################################
def get_project_id(file_name):
    """ Extracts project ID from intput BAM filename.

    :param file_name: string    e.g. "/PITT_0452_AHG2THBBXY_A1___P10344_C___13_cf_IGO_10344_C_20___hg19___MD.bam"
    :return: string             e.g. "P10344_C"
    """
    regex = "(?<=___)P[0-9]{5}[_A-Z,a-z]*(?=___)"  # Valid project ID is "P" + 5 numbers + (optional) [ "_" + 2 letters]
    matches = re.findall(regex, file_name)
    if len(matches) == 0:
        print("ERROR: Could not find IGO ID in filename: %s with regex: \"%s\"" % (file_name, regex))
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
            resp = requests.get(url, auth=(config.LIMS_USER, config.LIMS_PASSWORD), verify=False)
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

#fingerprint('MICHELLE_0474_AH5L2MDSX3')