#!/bin/bash

dir=$1

for subD in ${dir}/*
do
  for file in ${subD}/*.mfa
  do
    seed=$(basename ${file} .mfa)
    mv ${file} ${subD}/${seed}.fa
  done
done
