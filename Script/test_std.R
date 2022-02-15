#Test if real stopcodon exists.

#setwd("~/Dropbox (ASU)/Indel_project/Script")
library(Biostrings)
library(stringr)
library(seqinr)

# file_name = "../test_human_mouse_rat/Raw_data.2.outgroups/cds_seq/ENSG00000084072.fa"
# inDir     = "../test_human_mouse_rat/Raw_data.2.outgroups/cds_seq/"
# ouFile = "../test_human_mouse_rat/Raw_data.2.outgroups/QC/noStop.txt"
stopCodon_add = function(inDir, ouFile){
  
  File.holder = list.files(inDir, full = TRUE)
  
  noStop = c()
  for (i in File.holder){
    dna = readDNAStringSet(i, format = "fasta")
    wid = width(dna)
    for(j in 1:length(dna)){
      seq  = toString(dna[[j]], width = NULL)
      Stop = str_sub(seq, start = wid[j] - 2, wid[j])
      if(Stop != "TAG" && Stop != "TGA" && Stop != "TAA"){
        noStop = c(noStop, str_extract(basename(i), "[^.]+"))
        break
      }
    }
  }
  cat(noStop, file = ouFile, sep = "\n", append = FALSE)
}

args = commandArgs(trailingOnly = TRUE)
stopCodon_add(args[1], args[2])