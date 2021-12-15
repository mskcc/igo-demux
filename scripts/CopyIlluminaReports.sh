#!/bin/bash

DIR=/igo/staging/FASTQ
cd $DIR

echo "Searching for recently completed demuxes in $DIR"
FASTQ_DIRS=$(find /igo/staging/FASTQ -mindepth 1 -maxdepth 1 -cmin -4000)
for fastq_dir in ${FASTQ_DIRS}; do
  echo "Processing $fastq_dir" # for example: /igo/staging/FASTQ/MICHELLE_0465_AH3GKKDSX3
  runFullName="$(echo $fastq_dir| cut -d'/' -f 5)"

  echo "For Run $runFullName"
  # for bcl2fastq reports are located in /igo/staging/FASTQ/AYYAN_0105_000000000-K4R7C/Reports/html/000000000-K4R7C/all/all/all/
  htmlName="/Reports/html/*/all/all/all/laneBarcode.html"
  filename=$fastq_dir$htmlName
  
  homedir="/home/igo/html/"
  html="_laneBarcode.html"
  htmlnewfile="_laneBarcode_.html"
  copiedname=$homedir$runFullName$html
  copiednamePlus=$homedir$runFullName$htmlnewfile

  reports="/Reports/"
  dragen_reports_dir=$fastq_dir$reports
  dragen_replay=$fastq_dir"/dragen-replay.json"
  echo $dragen_replay
  
  # Check for bcl2fastq laneBarcode.html existence
  echo "Checking for existence of $filename"
  if [ -f $filename ]; then
    echo "Copying $filename to $copiedname"
    cp -p $filename $copiedname

    /opt/common/CentOS_7/python/python-3.7.1/bin/python3 /igo/work/igo/igo-demux/scripts/barcodelookup.py $copiedname /igo/work/igo/igo-demux/scripts/Barcodes.json
    toDir=/srv/www/sequencing-qc/static/html/FASTQ/
    toName=$toDir$runFullName$html

    touch $copiednamePlus -r $filename #set correct timestamp on new html file
    echo "scp $copiednamePlus to $toName"
    scp -p $copiednamePlus igo@igo:$toName
  elif [ -f $dragen_replay ]; then
    # This was a DRAGEN demux
    toDir=/srv/www/sequencing-qc/static/html/FASTQ/
    toNameLocal=$homedir$runFullName$html
    toNameRemote=$toDir$runFullName$html
    echo "Converting DRAGEN reports in folder $dragen_reports_dir to name $toNameLocal"
    python /igo/work/igo/igo-demux/scripts/dragen_csv_to_html.py $dragen_reports_dir $toNameLocal

    touch $toNameLocal -r $dragen_replay #set correct timestamp on new html file from a DRAGEN demux file
    echo "scp $toNameLocal to $toNameRemote"
    scp -p $toNameLocal igo@igo:$toNameRemote
  else
    echo "Failed to copy demux reports for run $fastq_dir"
  fi
done