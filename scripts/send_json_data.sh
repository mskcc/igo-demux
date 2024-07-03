#!/bin/bash


##Usage sh send_json_data.sh patth_to_json json_file  
args=$@

path_to_json=$1
json_file=$2

cd $path_to_json

echo $path_to_json
json_data=$(cat $json_file)
echo $json_data


curl -d "$json_data" -H "Content-Type: application/json" -X POST "http://igodb.mskcc.org:8080/ngs-stats/saveCellRangerSample"

