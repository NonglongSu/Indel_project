#!/bin/bash
dir1=$1
dir2=$2

for del in ${dir1}/*
do
	seed=$(basename "${del}" .txt)
	echo ${seed}
	while read -r file; do
		rm  "${dir2}/${seed}/${file}*"
	done < "${del}"
done
