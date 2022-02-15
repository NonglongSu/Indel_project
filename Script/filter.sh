#!/bin/bash
file1=QC/multi_no_3.txt
file2=nameTmp.txt
file3=nameList.txt
dir=$1

a=$( test -f $file1 ;echo $?)
b=0

if [[ $a == $b ]]
then 
	echo "There are protein coding sequences are not multiple of three."
	ls $dir/ | sed -e 's/.fa//g' > $file2
	cat $file1 | xargs -I % sh -c 'find cds_seq/ -name %.fa -delete '	
	awk 'NR==FNR{a[$0]=1;next}!a[$0]' $file1 $file2 > $file3
	rm $file2	
else
	echo "Nothing needs to be changed."
	ls $dir/ | cut -f 1 -d '.' > $file3

fi

