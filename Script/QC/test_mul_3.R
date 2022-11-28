#Record any files that are "not_multiple of 3". 

library(Biostrings)
library(scales)
library(stringr)

#setwd("~/Dropbox (ASU)/Indel_project/Script")

#dir    = "../test_human_mouse_rat/Raw_data/cds_seq/"
#output = "../test_human_mouse_rat/Raw_data/multi_no_3.txt"

find_mul_no_3 = function(dir, ouFile){
  
  allData = list.files(dir, full = TRUE)
  for(i in allData){
    cds   = readDNAStringSet(i, format = "fasta")
    width = width(cds)
    perc  = sum(width %% 3 == 0)/length(cds)
    
    if(perc < 1){
        name = str_extract(basename(i), "[^.]+")
        cat(name, file = ouFile, sep = "\n", append = TRUE)
    }
  }
}

args = commandArgs(trailingOnly = TRUE)
find_mul_no_3(args[1],args[2])
