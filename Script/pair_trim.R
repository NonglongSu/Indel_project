#Trim reference 
suppressPackageStartupMessages(library(Biostrings))

#setwd("~/Dropbox (ASU)/Indel_project/test_human_mouse_rat/Raw_data")

# inD = "mapped_cds_mafft"
# ouD = "pair_mafft/"
pair = function(inD,ouD){
  data = list.files(inD,full.names=T)
  for(i in 1:length(data)){
    dna      = readBStringSet(data[i])
    dna.pair = dna[(length(dna)-1):length(dna)] 
    writeXStringSet(dna.pair, paste0(ouD,basename(data[i])))  
    print(i)
  }
}

args = commandArgs(trailingOnly=TRUE)
pair(args[1],args[2])
