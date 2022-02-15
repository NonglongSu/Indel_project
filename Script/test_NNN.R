# Record any files that have NNNs.

library(Biostrings)
library(stringr)
library(seqinr)
library(stringr)

# setwd("~Dropbox (ASU)/Indel_project/Script/")

# inDir  = "../test_human_mouse_rat/Raw_data/cds_seq/"

test_N = function(inDir, ouFile){
    File.holder = list.files(inDir,full=TRUE)
    
    N.list = c()
    for(i in 1: length(File.holder)){
        dna   = readDNAStringSet(File.holder[i], format = "fasta")
        len   = length(grep('N', dna, ignore.case = TRUE))
        if(len > 0){
           N.list = c(N.list, str_extract(basename(File.holder[i]), "[^.]+"))
        }
    }
    cat(N.list, file = ouFile, sep = "\n", append = FALSE)
}

args = commandArgs(trailingOnly=TRUE)
test_N(args[1], args[2])
