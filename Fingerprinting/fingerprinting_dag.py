import subprocess
from airflow import DAG
from airflow.models.variable import Variable

from airflow.operators.bash import BashOperator

with DAG(
    dag_id='fingerprinting',
    tags=['fingerprinting'],
) as dag:

    def fingerprint(cmo_id):
        bams = Variable.get("Ready_to_fingerprint")
        bam_path=subprocess._FILE()
        vcf_output_directory=""
        Variable.set("Ready_to_fingerprint", )

        #find all bams
        for bam in bams:
            print("Running extract fingerprint: " + command1)
            haplotype_map=
            input_bam=
            output_vcf=
            reference_seq=
            patient_id=
            command1 = "gatk ExtractFingerprint --HAPLOTYPE_MAP \'{}\'  --INPUT \'{}\' --OUTPUT \'{}\' --REFERENCE_SEQUENCE \'{}\' --SAMPLE_ALIAS \'{}\'".format(haplotype_map, input_bam, output_vcf, reference_seq, cmo_id)
            subprocess.run(command1, shell=True, check=True)

            print("Running cross-check fingerprint: " + command2)
            command2 = "gatk CrosscheckFingerprints LOD_THRESHOLD=-5.0 CROSSCHECK_BY=FILE NUM_THREADS=30 OUTPUT=crosscheck_fingerprint.tsv HAPLOTYPE_MAP=\'{}\' INPUT=\'{}\'".format(haplotype_map, output_vcf)
            subprocess.run(command1, shell=True, check=True)
        

    fingerprint=BashOperator(
        task_id='fingerprinting',
        bash_command=fingerprint,
    )
    fingerprint

if __name__ == "__main__":
    dag.cli()
