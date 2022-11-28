#!/bin/bash

Dir1=$1
Dir2=$2

for file in ${Dir1}/*
do
  seed=$(basename ${file})
  echo ${seed}
  awk '/^>/ {n++} n>1 {print}' ${file} > ${Dir2}/${seed}
done
