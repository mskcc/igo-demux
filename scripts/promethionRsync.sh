#!/bin/bash
# A script to pull incremental backups from igo-prom using rsync
# rsync -avz prom@PC24B217:/data/20211122/* ./promethion/

set -o errexit
set -o nounset
set -o pipefail

readonly SOURCE_DIR="prom@igo-prom:/data/"
readonly DESTINATION_DIR="/igo/delivery/share/promethion/backup"
readonly DATETIME="$(date '+%Y-%m-%d_%H:%M:%S')"
readonly DESTINATION_PATH="${DESTINATION_DIR}/${DATETIME}"

# mkdir -p "${DESTINATION_DIR}"

# -r and -t options–recursive copy and time stamp copy
# -a to preserve the most important attributes of the source files 
# -v option to make the command more verbose

rsync -ravzu \
  "${SOURCE_DIR}/" \
  "${DESTINATION_PATH}"