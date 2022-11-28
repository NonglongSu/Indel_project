#Record any files that are not divisible by 3, ambiguous N, early stop codon, aa file.
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(seqinr))

#setwd("~/Dropbox (ASU)/Indel_project/Script/90")

testAll = function(f){
  mul3v = c()
  nv    = c()
  stopv = c()
  
  for(j in f){
    dna = readDNAStringSet(j)
    len = width(dna)
    nseq= length(dna)
    
    dnas = str_split(dna,'')
    len1 = length(which(dnas[[1]] %in% DNA_BASES))
    len2 = length(which(dnas[[2]] %in% DNA_BASES))
    
    if(all(c(len1,len2)!=len)){
      nv=c(nv,basename(j))
    }
    
    perc = sum(len%%3 == 0)/nseq
    if(perc!=1){mul3v = c(mul3v,basename(j))}
    
    # len.N = length(grep('N',dna,ignore.case=TRUE))
    # len.q = length(grep('\\?',dna,ignore.case=TRUE))
    # if(len.N+len.q>0){nv=c(nv,basename(j))}
    
    for(k in 1:nseq){
      pos   = gregexpr(paste(stop.codon,collapse="|"), dna[[k]])[[1]]
      pos.1 = pos[pos%%3 == 1]
      if((length(pos.1)>0) && (pos.1[1]<(len[k]-2))){#1st appearing stop codon
        stopv = c(stopv,basename(j))
        break
      }
    }
  }
  res = unique(c(mul3v,nv,stopv))
  if(!is.null(res)){
    return(res)
  }
}


#dir = "../../test_90_species/Raw_data/cds"
#ouD = "../../test_90_species/Raw_data/QC/"

main = function(dir, ouD){
  stop.codon  <<- c("TAG","TGA","TAA")
  Dirs = list.files(dir, full.names=TRUE)
  for (i in 1:length(Dirs)) {
    ouF    = paste0(ouD, basename(Dirs[i]), ".txt")
    Files  = list.files(Dirs[i], full.names=TRUE)
    res    = testAll(Files)
    
    cat(res,file=ouF,sep="\n",append=F)
    print(i)
  }
}


args = commandArgs(trailingOnly = TRUE)
main(args[1], args[2])