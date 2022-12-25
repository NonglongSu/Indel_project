#Reverse the sequence(AA)

#setwd("~/Dropbox (ASU)/Indel_project/test_human_mouse_rat/Script")

library(Biostrings)
library(seqinr)
library(stringr)

# file_name = "../Raw_data/aa_seq/ENSG00000000460.fa"
# ouDir    = "../Raw_data/Tmp/"

dna_to_aa = function(file_name,ouDir){
  
  aa = readAAStringSet(file_name,format="fasta")
  name = names(aa)
  
  tmp     = str_split(aa, "") 
  rev.aa  = list()
  for(i in 1:length(aa)){
    rev.tmp = rev(tmp[[i]])
    rev.aa[[i]] = paste0(rev.tmp,collapse="")
  }
  
  
  write.fasta(sequences=rev.aa, names=name, nbchar=80,
              open="w", as.string=TRUE, file.out=paste0(ouDir,basename(file_name)) )
}

####################################################
args = commandArgs(trailingOnly = TRUE)
dna_to_aa(args[1],args[2])
