#!/bin/bash
dir1=$1
dir2=$2

for del in ${dir1}/*.max.txt
do
	seed=$(basename "${del}" .max.txt)
	echo ${seed}
	while read -r file; do
		rm  ${dir2}/${seed}/${file}.fa*
	done < "${del}"
done
