
#!/bin/bash

DIR=/igo/staging/FASTQ
cd $DIR

FASTQ_DIRS=$(find /igo/staging/FASTQ -mindepth 1 -maxdepth 1 -cmin -123)
for fastq_dir in ${FASTQ_DIRS}; do
  dirName="$(echo $fastq_dir| cut -d'/' -f 2)"
  runLong="$(echo $fastq_dir| cut -d'_' -f 3)"
  if [[ "$runLong" == *"-"* ]]; then  # dir names are either like PITT_0276_BH2JK3BBXY or VIC_2421_000000000-C3W3G
    runName=$runLong
  elif [[ -z $(echo ${dirname} | grep PEPE) ]]; then
    # Pepe directory structure is slightly different, e.g. /igo/sequencers/pepe/output/210623_PEPE_8_AAAKYHCM5/
    runName=${runLong:0:9}
  else
    runName=${runLong:1:9}
  fi

  echo "For Run $runName whole name $runLong"
  part1="/Reports/html/"
  part2="/all/all/all/laneBarcode.html"
  filename=$fastq_dir$part1$runName$part2
  
  homedir="/home/igo/html/"
  html="_laneBarcode.html"
  htmlnewfile="_laneBarcode_.html"
  copiedname=$homedir$runName$html

  csv="/Reports/Demultiplex_Stats.csv"
  filename_dragen=$fastq_dir$csv
  dragen_replay=$fastq_dir"/dragen-replay.json"
  echo $dragen_replay
  
  # Check for bcl2fastq laneBarcode.html existence
  if [ -f $filename ]; then
    echo "Copying $filename to $copiedname"
    cp -p $filename $copiedname

    /opt/common/CentOS_7/python/python-3.7.1/bin/python3 /igo/work/igo/igo-demux/scripts/barcodelookup.py $copiedname /igo/work/igo/igo-demux/scripts/Barcodes.json
    toDir=/srv/www/sequencing-qc/static/html/FASTQ/
    toName=$toDir$dirName$html
    enrichedName=$homedir$runName$htmlnewfile

    touch $enrichedName -r $filename #set correct timestamp on new html file
    echo "scp $enrichedName to $toName"
    scp -p $enrichedName igo@igo:$toName
  elif [ -f $dragen_replay ]; then
    # This was likely a DRAGEN demux
    cp -p $filename_dragen $copiedname
    python3 /igo/work/igo/igo-demux/scripts/dragen_csv_to_html.py $copiedname

    toDir=/srv/www/sequencing-qc/static/html/FASTQ/
    toName=$toDir$dirName$html
    enrichedName=$homedir$runName$htmlnewfile

    touch $enrichedName -r $filename_dragen #set correct timestamp on new html file
    echo "scp $enrichedName to $toName"
    scp -p $enrichedName igo@igo:$toName
  else
    echo "Failed to copy demux reports"
  fi
done