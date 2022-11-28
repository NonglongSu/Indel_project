#!/bin/bash
#Usage: bash ../Script/90/check_file_num.sh Raw_data/cds/


Target=$1

for subDir in ${Target}*
do
	echo "$(basename ${subDir})"
	ls ${subDir} | wc -l
done
