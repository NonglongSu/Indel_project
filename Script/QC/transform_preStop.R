#Replace the pre stop codon symbol(*) with X. 

library(readr)
library(Biostrings)
library(stringr)
library(seqinr)
#setwd("~/Dropbox (ASU)/Indel_project/Script")

# file = "../test_human_mouse_rat/Raw_data/QC/preStop.txt"
# dir  = "../test_human_mouse_rat/Raw_data/aa_seq/"

stopCodon_fix = function(file, dir){
  
  File = read_delim(file, "\t", col_names = FALSE)
  
  for (i in 1:nrow(File)) {
    filename = paste0(dir, File[[i, 1]], ".fa")
    aa = readAAStringSet(filename, format = "fasta")
    
    aa.new = list()
    for (j in 1:length(aa)) {
      aa.new[[j]] = str_replace_all(toString(aa[[j]]), "\\*", "X")
    }
    
    write.fasta(aa.new, names(aa), filename,
                open = "w", nbchar = 80, as.string = FALSE)
    
  }
  
  
  
}

args = commandArgs(trailingOnly = TRUE)
stopCodon_fix(args[1], args[2])


