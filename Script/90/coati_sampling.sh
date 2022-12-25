#!/bin/bash
dir1=$1
dir2=$2
file=$3

while read -r line
do
  pars=$line
done < $file

seed1=$(basename "${dir2}")
echo ${seed1}

for f in ${dir1}/${seed1}/*
do
  seed2=$(basename "${f}")
  ~/coati/builddir/src/coati-sample ${f} -o ${dir2}/${seed2}.json -k 3 -n 100 ${pars}
done
