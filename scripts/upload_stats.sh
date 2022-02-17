#!/bin/bash
# Uploads stat (.txt) files to ngs-stats & LIMS

# Run: 
#   RUNNAME=MICHELLE_0480_AH5KTWDSX3_WGS ./upload_stats.sh

MACHINE=$(echo ${RUNNAME} | cut -d'_' -f1)  # MICHELLE_0347_BHWN55DMXX -> MICHELLE
RUN_NUM=$(echo ${RUNNAME} | cut -d'_' -f2)  # MICHELLE_0347_BHWN55DMXX -> 0347
FLOWCELL=$(echo ${RUNNAME} | cut -d'_' -f3) # MICHELLE_0347_BHWN55DMXX -> BHWN55DMXX
STAT_PREFIX="${MACHINE}_${RUN_NUM}_${FLOWCELL}"

echo "START: $(date +"%D %T")"
DELPHI_ENDPOINT="http://delphi.mskcc.org:8080/ngs-stats/picardstats/updaterun/${MACHINE}/${STAT_PREFIX}"
echo "Updating Picard stats DB: ${DELPHI_ENDPOINT}"
curl "${DELPHI_ENDPOINT}"
sleep 10
LIMS_ENDPOINT="https://igo-lims02.mskcc.org:8443/LimsRest/updateLimsSampleLevelSequencingQc?runId=${STAT_PREFIX}"
printf "\nUpdating LIMS QC Metrics: ${LIMS_ENDPOINT}\n"
curl -k ${LIMS_ENDPOINT}
printf "\nEND: $(date +"%D %T")\n"