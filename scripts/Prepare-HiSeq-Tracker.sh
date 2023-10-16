#!/bin/bash
#Example: sh ./Prepare-HiSeq-Tracker.sh MOMO_0125

RUN=$1
SAMPLESHEET="/rtssdc/mohibullahlab/LIMS/LIMS_SampleSheets/SampleShee*"$RUN"*.csv"

#parse sample sheet to parameters

ROWS=$(sed -n '/Reads/,/Settings/p' $SAMPLESHEET | wc -l)

#echo $ROWS

if [[ "$ROWS" < 5  ]]; then
  RUNTYPE="SE"
else
  RUNTYPE="PE"
fi
echo $RUNTYPE

ROWS1=$(sed -n '/Reads/,/Settings/p' $SAMPLESHEET )
echo $ROWS1
RUNNAME=$(echo $RUN | awk '{pos=match($0,"_"); print (substr($0,pos+1,length($0)))}')

#make summary files that contains only Project, species and recipe
awk '{if(found) print} /Lane/{found=1}' $SAMPLESHEET | awk 'BEGIN { FS = "," } ;{printf"%s,%s,,%s,%s\n",$8,$9,$5,$4}' | sort | uniq | awk '{sub(/Project_/,"");print $0}'