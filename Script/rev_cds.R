#Reverse the CDS [for SW power test]

#setwd("~/Dropbox (ASU)/Indel_project/Script")

suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(seqinr))
suppressPackageStartupMessages(library(stringr))

file_name = "../test_human_mouse_rat/Raw_data/cds_seq/ENSG00000000460.fa"
ouDir     = "../test_human_mouse_rat/Rev_seq/"

rev_cds = function(file_name,ouDir){
  
  cds = readDNAStringSet(file_name,format="fasta")
  name = names(cds)
  
  cds1   = str_split(cds, "") 
  cds.r  = list()
  for(i in 1:length(cds)){
    cds.r[[i]] = paste0(rev(cds1[[i]]), collapse="")
  }
  
  write.fasta(sequences=cds.r, names=name, nbchar=80,
              open="w", as.string=TRUE, file.out=paste0(ouDir, basename(file_name)))
}

####################################################
args = commandArgs(trailingOnly = TRUE)
dna_to_aa(args[1], args[2])
