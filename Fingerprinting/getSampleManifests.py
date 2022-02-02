import json
import os
import re
import sys
import requests

import config

MAPPED_FIELDS = ["cmoPatientId", "tumorOrNormal"]  # These will be joined in the bam. Order matters

# Delimiter for BAM name. Each field in the name can be parsed out on this unique string
BAM_DELIMITER = "___"
igoId_to_correctedCmoPatientId_map = {
    # { IGO_ID -> CORRECTED }
}
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

def fail(err_msg = None):
    """ Exits Application and logs error message """
    print("Usage is 'python ./create_merge_commands ${FILE_TO_WRITE_COMMANDS} ${DIR_TO_WRITE_BAMS} [--t] [--s=/path/to/samtools] ${BAM_FILES}'")
    if(err_msg):
        print("ERROR: " + err_msg)
    sys.exit(1)

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
    :return: string             e.g. "10344_C"
    """
    regex = "IGO_([a-zA-Z0-9_.-]*?)___"
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
     