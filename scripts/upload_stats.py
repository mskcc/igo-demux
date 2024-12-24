import requests
import sys

def upload_stats(sequencer_and_run):
    sequencer_and_run_prefix = "_".join(sequencer_and_run.split("_")[0:3])
    sequencer_name = sequencer_and_run.split("_")[0]

    DB_ENDPOINT="http://igodb.mskcc.org:8080/ngs-stats/picardstats/updaterun/{}/{}".format(sequencer_name, sequencer_and_run_prefix)
    PICARD_LIMS_ENDPOINT="https://igo-lims02.mskcc.org:8443/LimsRest/updateLimsSampleLevelSequencingQc?runId={}".format(sequencer_and_run_prefix)
    TENX_LIMS_ENDPOINT="https://igo-lims02.mskcc.org:8443/LimsRest/updateTenXSampleLevelStats?runId={}".format(sequencer_and_run_prefix)

    print(DB_ENDPOINT)
    print(PICARD_LIMS_ENDPOINT)
    # push data into database
    r = requests.get(DB_ENDPOINT)
    print(r.text)
    # check if successful
    if r.status_code == requests.codes.ok:
        print("Picard stats posted to stats DB")
    else:
        r.raise_for_status()

    r2 = requests.get(PICARD_LIMS_ENDPOINT, verify=False)
    if r2.status_code == requests.codes.ok:
        print("Update LIMS QC matrics done")
    else:
        r2.raise_for_status()

    r3 = requests.get(TENX_LIMS_ENDPOINT, verify=False)
    if r3.status_code == requests.codes.ok:
        print("Update 10X LIMS QC matrics done")
    else:
        r3.raise_for_status()
    

if __name__ == "__main__":
    sequencer_and_run = sys.argv[1]

    upload_stats(sequencer_and_run)