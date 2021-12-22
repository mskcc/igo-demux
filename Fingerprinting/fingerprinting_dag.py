import subprocess
from sys import path_importer_cache, stdout
from time import process_time
from airflow import DAG
from airflow.models.variable import Variable

from airflow.operators.bash import BashOperator
from airflow.operators.python import PythonOperator

# with DAG(
#     dag_id='fingerprinting',
#     tags=['fingerprinting'],
# ) as dag:

def fingerprint(patient_id, run):
    STATS_DIR = '/igo/staging/stats/'
    REFERENCE_SEQUENCE_DIR = '/igo/work/genomes/H.sapiens/GRCh38.p13/GRCh38.p13.dna.primary.assembly.fa'
    HAPLOTYPE_MAP = "/home/igo/fingerprint_maps/map_files/hg38_igo.map"
    EXECUTION_DIR = STATS_DIR + run
    t_start = process_time()
    vcfs = []
    subprocess.chdir(EXECUTION_DIR) 
    #find all bams
    input_bams = []
    print("Finding bams of the run argument...")
    input_bams = subprocess.run('find . -maxdepth 2 -name *' + run + '*.bam"')
    input_bams = input_bams.stdout.split('\n')
    subprocess.run('mkdir ./VCF/')
    for bam in input_bams:
        output_vcf = STATS_DIR + run + '/VCF/' + patient_id + '_' + run + '.vcf'
        subprocess.chdir(EXECUTION_DIR)    
        command1 = "bsub gatk ExtractFingerprint --HAPLOTYPE_MAP \'{}\'  --INPUT \'{}\' --OUTPUT \'{}\' --REFERENCE_SEQUENCE \'{}\' --SAMPLE_ALIAS \'{}\'".format(HAPLOTYPE_MAP, bam, output_vcf, REFERENCE_SEQUENCE_DIR, patient_id)
        print("Running extract fingerprint: " + command1)
        subprocess.run(command1, shell=True, check=True)
        vcfs.append(output_vcf)

    for vcf in vcfs:
        command2 = "bsub gatk CrosscheckFingerprints LOD_THRESHOLD=-5.0 CROSSCHECK_BY=FILE NUM_THREADS=30 OUTPUT=crosscheck_fingerprint.tsv HAPLOTYPE_MAP=\'{}\' INPUT=\'{}\'".format(HAPLOTYPE_MAP, vcf)
        print("Running cross-check fingerprint: " + command2)
        subprocess.run(command1, shell=True, check=True)

    t_stop = process_time()
    print("Elapsed time to fingerprint: ", t_stop - t_start)
        
"""""
fingerprint=PythonOperator(
    task_id='fingerprinting',
    bash_command=fingerprint,
)
fingerprint
"""""

if __name__ == "__main__":
    #dag.cli()
    fingerprint()
