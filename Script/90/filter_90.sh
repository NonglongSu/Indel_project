#!/bin/bash

#file1=QC/multi_no_3.txt

dir1=$1
dir2=$2

for file in ${dir1}/*
do
	if test -f "${file}"
	then
		seed=$(basename "${file}" .txt)
		echo  "${seed} contains CDS not multiple of 3."
		cat ${file} | xargs -I % sh -c 'find ${dir2}/${seed}/ -name %.fasta -delete'
	else
		echo "Nothing changed"
	fi
done





