#Record any files that are "not_multiple of 3". 

library(Biostrings)
library(scales)
library(stringr)

#setwd("~/Dropbox (ASU)/Indel_project/Script/90")

# dir    = "../../test_90_species/Raw_data/cds"
# ouD    = "../../test_90_species/Raw_data/QC/Non3/"

find_mul_no_3 = function(f, output){
  
  for(j in f){
    cds   = readDNAStringSet(j, format = "fasta")
    width = width(cds)
    perc  = sum(width %% 3 == 0)/length(cds)
    
    if(perc < 1){
        name = sub(".fasta", "", basename(j))
        cat(name, file = output, sep = "\n", append = TRUE)
    }
  }
}


main = function(dir, ouD){
  Dirs = list.files(dir, full.names = TRUE)
  for (i in 1:length(Dirs)) {
    ouFile = paste0(ouD, basename(Dirs[i]), ".txt")
    Files  = list.files(Dirs[i], full.names = TRUE)
    print(basename(Dirs[i]))
    find_mul_no_3(Files, ouFile)
  }
  
}



args = commandArgs(trailingOnly = TRUE)
main(args[1], args[2])
