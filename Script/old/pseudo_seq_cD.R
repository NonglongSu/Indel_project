#Filter out all files to retain the one with highest score. 
#set up constant distance(24) as the denominator. 

library(Biostrings)
library(BiocGenerics)
library(parallel)
library(S4Vectors)
library(stringr)
library(tidyverse)
library(dplyr)
library(seqinr)
library(stringi)
#setwd("~/Dropbox (ASU)/Indel_project/test_human_mouse_rat/Script")

#Alignment algo (localized global-alignment)
#Hamming distance
Align = function(u,v){
  u = str_split(u,"")
  v = str_split(v,"")
  u = u[[1]]
  v = v[[1]]
  #Cal the identical number of identical bases not including gaps.
  id_pos = length(which(u[u == v] != '-'))
  #Cal the percent sequence identity
  sim = 100*(id_pos)/(constant_wall)
  return (sim)
}

#Action move
#left_sliding
left_slide = function(seq_v,index,wid){
  swap(seq_v[index-1],seq_v[index+wid-1])
  seq   = paste(seq_v,collapse = "")
  return (seq)
}

#right_sliding
right_slide = function(seq_v,index,wid){
  swap(seq_v[index],seq_v[index+wid])
  seq   = paste(seq_v,collapse = "")
  return (seq)
}

#Convert string to vector
str_convert = function(s){
  s1 = str_split(s,"")
  s1 = s1[[1]]
  return(s1)
}

#Generate the pseudo seq with an optimal alignment
Merge = function(idx,wid,x,y){
  
  #assume current state is the optimal state. 
  best_aligned = x 
  
  #Original window seq
  start = idx-(constant_wall/2)
  stop  = idx+(constant_wall/2)+wid-1
  wid.1   = substr(x, start, stop)
  wid_ref = substr(y, start, stop)
  
  #Original alignment score
  sim_ref_public = Align(wid_ref,wid.1) 
  
  #left_aligning
  sim_total.l   = c()
  all_aligned.l = c()
  new_aligned.l = x
  idx.l         = idx
  for(k in 1:wall){
    new_aligned.l = left_slide(str_convert(new_aligned.l),idx.l,wid)
    wid.L=substr(new_aligned.l, start, stop)
    sim_total.l[k]= Align(wid_ref,wid.L) 
    all_aligned.l[k] = new_aligned.l
    idx.l =idx.l-1
  }
  
  #right_aligning
  sim_total.r   = c()
  all_aligned.r = c()
  new_aligned.r = x
  idx.r         = idx
  for(k in 1:wall){
    new_aligned.r = right_slide(str_convert(new_aligned.r),idx.r,wid)
    wid.R=substr(new_aligned.r, start, stop )
    sim_total.r[k] = Align(wid_ref,wid.R) 
    all_aligned.r[k] = new_aligned.r
    idx.r =idx.r+1
  }
  
  #Filtering multiple optimal allignment ---/+++
  sim_all   = c(sim_ref_public,sim_total.l,sim_total.r)
  alg_all   = c(x,all_aligned.l,all_aligned.r)
  index.sim = which(sim_all==max(sim_all))
  if(length(index.sim) > 1){
    subseq(best_aligned,start=idx, end=idx+wid-1) = stri_rand_strings(1,wid,'[+]')
  }else{
    best_aligned = alg_all[index.sim]
  }
  return (best_aligned)
}


#inFile = "Test_pseudo_seq/test.fa"
#inFile = "../Data/updated_cds/ENSG00000000460.fa"

#seq1 = "AAT===AAACAAAGAATGCTTACTGT---ATAAGGCTTACTGTTCTAGCG===ATCACCGCG===TCATGTCTAGTTATGAACGGC------GGTTTAACATTGAATAGCAAGGCACTTCCA---TAATAGGGCCGTC===GTAATTGTCTAATATAG------ATAGTA==="
#seq2 = "TAA------AA===AATTTGATGCTACATTGGATGAGTCTACTTCGAGCGCGCCGCATCGATTGCAAGAGCAGTGTTGCCT===AAGAGCCGTTAGATGCGTCGTTG---ATCGCGTCCGATAATTCGGGAGTTGTGC===CCCAATATTTAATATGATGA===TAGCTATAA"



main = function(inFile,ouDir,num){
  
  #oudir="../Data/mapped_cds_Out/"
  
  dna = readBStringSet(inFile,format = "fasta")
  name = names(dna)
  
  #set up a window (global var)
  num = as.numeric(num)
  wall <<- num
  constant_wall <<- 24
  
  #String mode
  spec.1 = toString(dna[[1]],width=NULL)
  spec.2 = toString(dna[[2]],width=NULL)
  
  #Capture all gaps of different lengths(3,6,9,12)
  dna.1 = str_split(as.character(dna), '')
  g = lapply(dna.1, function(x) { IRanges(x == '-')})
  g = IRangesList(g)
  
  m = g[[1]]
  r = g[[2]]
  
  wid.m = width(m)
  wid.r = width(r)
  
  l.m = length(wid.m)
  l.r = length(wid.r)
  
  pos.m = start(m)
  pos.r = start(r)
  
  
  if(l.m>0 & l.r>0){
    for(i in 1:l.m){
      spec.1 = Merge(pos.m[i],wid.m[i],spec.1,spec.2)
    }
    for(i in 1:l.r){
      spec.2 = Merge(pos.r[i],wid.r[i],spec.2,spec.1)
    }
  }else if(l.m>0 & l.r==0){
    for(i in 1:l.m){
      spec.1 = Merge(pos.m[i],wid.m[i],spec.1,spec.2)
    }
  }else if(l.m==0 & l.r>0){
    for(i in 1:l.r){
      spec.2 = Merge(pos.r[i],wid.r[i],spec.2,spec.1)
    }
  }else{
    quit(save = "default", status = 0, runLast = TRUE)
  }
  
  new_seq = list(spec.1,spec.2)
  
  # Write out to disk
  write.fasta(sequences = new_seq,names = name, nbchar=80,
              open = "w",as.string = TRUE, file.out = paste0(ouDir,basename(inFile)) )
  
  # write.fasta(sequences = new_seq,names = name, nbchar=80,
  #             open = "w",as.string = TRUE, file.out = "Test_pseudo_seq/test_ouput.fa" )
  
}

#compatible with run_test.R
#if(interactive()){
args = commandArgs(trailingOnly = TRUE)
main(args[1],args[2],args[3])
#}

################


