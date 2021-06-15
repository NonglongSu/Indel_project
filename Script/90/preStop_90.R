# Replace the pre stop codon symbol(*) with X. 

library(readr)
library(Biostrings)
library(stringr)
library(seqinr)
#setwd("~/Dropbox (ASU)/Indel_project/Script")

# dir  = "../../test_90_species/Raw_data/QC/Stop"
# ouD  = "../../test_90_species/Raw_data/aa/"



stopCodon_fix = function(filelist, fName, ouD){
  Filelist    = read_delim(filelist, "\t", col_names = FALSE)
  for (i in 1:nrow(Filelist)) {
    fasta.file = paste0(ouD, fName, "/", Filelist[[i, 1]], ".fasta")
    aa         = readAAStringSet(fasta.file, format = "fasta")
    aa.new     = lapply(aa, function(x){ str_replace_all(toString(x), "\\*", "X")})
    write.fasta(aa.new, names(aa), fasta.file,
                open = "w", nbchar = 80, as.string = FALSE)
  }
}

main = function(dir, ouD){
  Files      = list.files(dir, full.names = TRUE)
  for (i in 1:length(Files)) {
    fileName = sub(".txt", "", basename(Files[i]))
    print(paste0(fileName, " fixing"))
    stopCodon_fix(Files[i], fileName, ouD)
  }
}


##############################################
args = commandArgs(trailingOnly = TRUE)
main(args[1], args[2])


