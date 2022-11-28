#Generate a csv/tsv file of 
# - #GC
# - seq length

suppressWarnings(suppressMessages(library(dplyr)))
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(seqinr))


#setwd("~/Dropbox (ASU)/Indel_project/test_90_species")


#########################################################
# inD = "Raw_data/cds/07_yeast_aligned_cds"
# ouF = "Results/GC/07_yeast_aligned_cds.tsv"

main = function(inD,ouF){

  Files = list.files(inD,full.names=T)
  n     = length(Files)
  
  datM           = matrix(0,n,4)
  colnames(datM) = c('C','G','lenA','lenB') 
  
  for (i in 1:length(Files)) {
    dna  = readDNAStringSet(Files[i])
    len  = width(dna)
   
    nucc = oligonucleotideFrequency(dna,width=1)
    nucc = colSums(nucc)
    
    datM[i,] = c(nucc[2],nucc[3],len[1],len[2])
    print(i)
  }
  
  #Output int tsv file
  write.table(datM, file=ouF,sep='\t', append=FALSE, quote=FALSE, row.names=FALSE, col.names=TRUE)
}


########################################
args = commandArgs(trailingOnly=TRUE)
main(args[1], args[2])