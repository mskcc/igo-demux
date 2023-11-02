#!/bin/bash
# Sends summary email for every samplesheet written to the LIMS samplesheet directory

cd /rtssdc/mohibullahlab/LIMS/LIMS_SampleSheets

find /rtssdc/mohibullahlab/LIMS/LIMS_SampleSheets/*.csv -mmin -60 > ~/RunTracker/Run_For_Tracker.txt

# try to make sure the csv file is not still being written by the LIMS
sleep 120 

chmod -R 775 ~/RunTracker/Run_For_Tracker.txt

for x in $(cat ~/RunTracker/Run_For_Tracker.txt); do
  RUNNAME=$(echo $x | awk '{sub(/SampleSheet_/,""); print $0}' | awk '{sub(/.csv/,""); print $0}')
  Content=$(sh ~/Scripts/Automate-Casava/Prepare-HiSeq-Tracker.sh $RUNNAME); 
  echo $RUNNAME | mail -s "New run ready for tracker $Content " mcmanamd@mskcc.org luc@mskcc.org naborsd@mskcc.org cobbsc@mskcc.org
done
