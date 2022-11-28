#!/bin/bash
# dir1=$1
# dir2=$2
#
# for subD in ${dir1}/*
# do
#   seed1=$(basename "${subD}")
#   #avoid repeated running
#   num1=$(ls ${subD} | wc -l)
#   num2=$(ls ${dir2}/${seed1} | wc -l)
#   if [ ${num1} -eq ${num2} ]; then
#     continue
#   else
#     echo ${seed1}
#   fi
#   for file in ${subD}/*
#   do
#     seed2=$(basename "${file}")
#     ~/coati/builddir/src/coati-alignpair ${file} -o ${dir2}/${seed1}/${seed2} -k 3
#   done
# done

dir1=$1
dir2=$2
file=$3

while read -r line
do
  pars=$line
done < $file

seed1=$(basename "${dir2}")
echo ${seed1}

for file in ${dir1}/${seed1}/*
do
  seed2=$(basename "${file}")
  ~/coati/builddir/src/coati-alignpair ${file} -o ${dir2}/${seed2} -k 3 ${pars}
done
