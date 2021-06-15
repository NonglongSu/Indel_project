#Translate DNA to Amino acid

library(Biostrings)
library(scales)
library(stringr)
#setwd("~/Dropbox (ASU)/Indel_project/Script/90")


# dir  = "../../test_90_species/Raw_data/cds"
# ouD = "../../test_90_species/Raw_data/aa/"

dna_to_aa = function(f, subD){
 for(i in f){
   cds     = readDNAStringSet(i, format = "fasta")
   protein = Biostrings::translate(cds, genetic.code = GENETIC_CODE, if.fuzzy.codon = "solve")
   
   writeXStringSet(protein, paste0(subD, basename(i)), append = FALSE,
                   compress = FALSE, compression_level = NA, format = "fasta")
 }  
}

main = function(dir, ouD){
  Dirs     = list.files(dir, full.names = TRUE)
  for (i in 1:length(Dirs)) {
    subdir = paste0(ouD, basename(Dirs[i]), "/")
    dir.create(subdir)
    Files  = list.files(Dirs[i], full.names = TRUE)
    print(basename(Dirs[i]))
    dna_to_aa(Files, subdir)
  }
}




#########################################################
args = commandArgs(trailingOnly = TRUE)
main(args[1], args[2])

