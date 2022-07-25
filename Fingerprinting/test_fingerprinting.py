import pytest
import fingerprinting_dag

def test_get_igo_id():
    igo_id = fingerprinting_dag.get_igo_id("PITT_0452_AHG2THBBXY_A1___P10344_C___13_cf_IGO_10344_C_20___hg19___MD.bam")
    assert(igo_id == "10344_C_20")

    # Current WGS naming convention for .bams
    igo_id = fingerprinting_dag.get_igo_id("MICHELLE_0533_BHWLFVDSX3___P05816_DR___FL001-164SC_IGO_05816_DR_5.bam")
    assert(igo_id == "05816_DR_5")