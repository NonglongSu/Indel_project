#Trim reference
suppressPackageStartupMessages(library(Biostrings))

#setwd("~/Dropbox (ASU)/Indel_project/test_human_mouse_rat/Rev_seq")

# inF = "Raw_data/%.mapped.fa"
# ouF = "Raw_data/%.paired.fa"
pair = function(inF,ouF){
  dna      = readBStringSet(inF, format="fasta")
  dna.pair = dna[(length(dna)-1):length(dna)] 
  writeXStringSet(dna.pair, ouF)  
}

args = commandArgs(trailingOnly=TRUE)
pair(args[1],args[2])