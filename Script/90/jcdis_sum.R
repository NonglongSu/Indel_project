#generaet a table of: species_pair,locus,lenA,lenB,match_count,mismatch_count,gapA_count,gapB_count,gapA_len,gapB_len

suppressWarnings(suppressMessages(library(readr)))
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(stats))
suppressPackageStartupMessages(library(dplyr))

#setwd("~/Dropbox (ASU)/Indel_project/test_90_species")

#read files
readF = function(Dir,spec.name){
  Files  = list.files(Dir,full.names=T)
  n      = length(Files)
  res    = matrix(0,n,8)
  locus  = c()
  spname = rep(spec.name,n)
  for(j in 1:n){
    dna       = readBStringSet(Files[j], format="fasta")
    locus[j]  = str_remove(basename(Files[j]),'\\.[^\\.]*$')
    wid       = width(dna)
    seqs      = str_split(dna,'')
    lenA      = nchar(str_remove_all(dna[[1]],'-'))
    lenB      = nchar(str_remove_all(dna[[2]],'-'))
    jc.dis    = jc_distance_sum(seqs)
    id.ct.len = indel_ct_len_sum(seqs)
    res[j,]   = c(lenA,lenB,jc.dis,id.ct.len)
  }
  
  res = as.data.frame(res)
  res = cbind(spname,locus,res)
  return(res)
}

#cal. jukes-cantor corrected distance
jc_distance_sum = function(seqs){
  gap_pos        = seqs[[1]] == "-" | seqs[[2]] == "-"
  match_count    = sum(seqs[[1]] == seqs[[2]] & !gap_pos)
  mismatch_count = sum(seqs[[1]] != seqs[[2]] & !gap_pos)
  return(c(match_count,mismatch_count))
}

#sum the indel length/count of each gene-pair
indel_ct_len_sum = function(seqs){
  countA=countB=lenA=lenB=0
  
  g    = IRangesList(lapply(seqs, function(x){IRanges(x=='-')}))
  countA = length(g[[1]])
  countB = length(g[[2]])
  if((countA!=0)){
    lenA = sum(width(g[[1]]))
  }
  if(countB!=0){
    lenB = sum(width(g[[2]]))
  }
  res = c(countA,countB,lenA,lenB)
  return(res)
}

#####################################################
# inD    = "Raw_data/align_max/01_FcaCaf_aligned_cds.max.tsv"
# ouF    = "Raw_data/JCdis_sum/01_FcaCaf_aligned_cds.max.tsv"
# pat    = "max"
main = function(inD,ouF,pat){
  spec.name     = basename(inD)
  spec_pair_sum = readF(inD,spec.name)
  colnames(spec_pair_sum) = c('species_pair','locus','lenA','lenB','match_count','mismatch_count','gapA_count','gapB_count','gapA_len','gapB_len')
  write.table(spec_pair_sum,ouF,quote=F,col.names=T,row.names=F,sep="\t")
}
#####################################
args=commandArgs(trailingOnly=T)
main(args[1],args[2],args[3])



