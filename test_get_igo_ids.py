import sys
import re

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


def test_get_igo_ids(): 
    res = get_igo_id("/PITT_0452_AHG2THBBXY_A1___P10344_C___13_cf_IGO_10344_20___hg19___MD.bam")
    assert(res == '10344_20')