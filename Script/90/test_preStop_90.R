# Record any files with early stop codon.
library(Biostrings)
library(stringr)
library(seqinr)
library(stringr)

#setwd("~/Dropbox (ASU)/Indel_project/Script/90")

# dir  = "../../test_90_species/Raw_data/cds"
# ouD = "../../test_90_species/Raw_data/QC/Stop/"

test_early_stop = function(f, output, stp){
  fileList = c()
  for( i in f){
    dna = readDNAStringSet(i, format = "fasta")
    len = width(dna)
    for(j in 1:length(dna)){
      pos   = gregexpr(paste(stp, collapse = "|"), dna[[j]])[[1]]
      pos.1 = pos[pos %% 3 == 1]
        if((length(pos.1) > 0) && (pos.1[1] < (len[j] - 2))){ # Always targetting the 1st stop codon
              #The early stop codon can be rare codon
              filename  = sub(".fasta", "", basename(i))
              fileList  = c(fileList, filename)
              break
        }
     }
  }
  if(!is.null(fileList)){
    cat(fileList, file = output, sep = "\n", append = FALSE)
  }
}

main = function(dir, ouD){
  Dirs        = list.files(dir, full.names = TRUE)
  stop.codon  = c("TAG", "TGA", "TAA")
  for (i in 1:length(Dirs)) {
    ouFile = paste0(ouD, basename(Dirs[i]), ".txt")
    Files  = list.files(Dirs[i], full.names = TRUE)
    print(basename(Dirs[i]))
    test_early_stop(Files, ouFile, stop.codon)
  }
}


########################################
args = commandArgs(trailingOnly = TRUE)
main(args[1], args[2])


