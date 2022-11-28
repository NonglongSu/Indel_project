#!/bin/bash
#usage: bash ../../../Script/90/create_dir.sh ../cds  (inside align/)


find $1/* -type d > dirs.txt


# while read -r file
# do
#   seed=$(basename ${file})
#   echo $seed >> Species.txt
# done < dirs.txt
# rm dirs.txt

while read -r file
do
  seed=$(basename ${file})
  echo $seed
  mkdir -p ${seed}
done < dirs.txt


rm dirs.txt
