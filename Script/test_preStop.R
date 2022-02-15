# Record any files with early stop codon.
library(Biostrings)
library(stringr)
library(seqinr)
library(stringr)

#setwd("~/Dropbox (ASU)/Indel_project/test_human_mouse_rat/Script")

# inDir  = "../test_human_mouse_rat/Raw_data/cds_seq/"
# ouFile = "../Raw_data.2.outgroups/early_stop.txt"
test_early_stop = function(inDir, ouFile){
  
  File.holder = list.files(inDir, full = TRUE)
  stop.codon  = c("TAG", "TGA", "TAA")
  fileList    = c()
  for( i in File.holder){
    dna = readDNAStringSet(i, format = "fasta")
    len = width(dna)
    for(j in 1:length(dna)){
      pos   = gregexpr(paste(stop.codon, collapse = "|"), dna[[j]])[[1]]
      pos.1 = pos[pos %% 3 == 1]
        if(length(pos.1) > 0 && pos.1[1] < (len[j]  -2) ){ # Always targetting the 1st stop codon
              #The early stop codon can be rare codon
              filename  = str_extract(basename(i), "[^.]+")
              fileList  = c(fileList, filename)
              break
        }
    }
  }
  cat(fileList, file = ouFile, sep = "\n", append = FALSE)
}


args = commandArgs(trailingOnly = TRUE)
test_early_stop(args[1], args[2])


