# Slide the gaps and find the best hit. 
suppressWarnings(suppressMessages(library(tidyverse)))
suppressWarnings(suppressMessages(library(dplyr)))
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(BiocGenerics))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(S4Vectors))
suppressPackageStartupMessages(library(seqinr))
suppressPackageStartupMessages(library(stringi))
suppressPackageStartupMessages(library(stringr))

#setwd("~/Dropbox (ASU)/Indel_project/test_human_mouse_rat/Data_6")

#Hamming distance
Align = function(u, v){
  u = str_split(u, "")
  v = str_split(v, "")
  u = u[[1]]
  v = v[[1]]
  
  id_pos    = length(which(u[u == v] != '-'))        # Cal the identical number of identical bases excluding gaps.
  align_pos = Wall * 2
  
  sim = 100 * (id_pos) / (align_pos)
  return (sim)
}

#left_sliding
left_slide = function(seq_v, index, wid){
  swap(seq_v[index - 1], seq_v[index + wid - 1])
  seq   = paste(seq_v, collapse = "")
  return (seq)
}

#right_sliding
right_slide = function(seq_v, index, wid){
  swap(seq_v[index], seq_v[index + wid])
  seq   = paste(seq_v, collapse = "")
  return (seq)
}

#Convert string to vector
str_convert = function(s){
  s1 = str_split(s, "")
  s1 = s1[[1]]
  return(s1)
}

#Generate the pseudo seq with an optimal alignment
Merge = function(idx, wid, x, y){
  best_aligned = x                    
  
  # Original window seq
  start = idx - Wall
  stop  = idx + Wall + wid - 1
  wid.1   = substr(x, start, stop)
  wid_ref = substr(y, start, stop)
  
  #Original alignment score
  sim_ref_public = Align(wid_ref, wid.1) 
  
  #left_aligning
  sim_total.l   = c()
  all_aligned.l = c()
  new_aligned.l = x
  idx.l         = idx
  for(k in 1:Window){
    new_aligned.l = left_slide(str_convert(new_aligned.l), idx.l, wid)
    wid.L         = substr(new_aligned.l, start, stop)
    
    sim_total.l[k]   = Align(wid_ref, wid.L) 
    all_aligned.l[k] = new_aligned.l
    idx.l            = idx.l - 1
  }
  
  #right_aligning
  sim_total.r   = c()
  all_aligned.r = c()
  new_aligned.r = x
  idx.r         = idx
  for(k in 1:Window){
    new_aligned.r = right_slide(str_convert(new_aligned.r), idx.r, wid)
    wid.R         = substr(new_aligned.r, start, stop)
    
    sim_total.r[k]   = Align(wid_ref, wid.R) 
    all_aligned.r[k] = new_aligned.r
    idx.r            = idx.r + 1
  }
  
  #Filtering multiple optimal allignment ---/+++
  sim_all   = c(sim_ref_public, sim_total.l, sim_total.r)
  alg_all   = c(x, all_aligned.l, all_aligned.r)
  index.sim = which(sim_all == max(sim_all))
  if(length(index.sim) > 1){
    subseq(best_aligned,start = idx, end = idx + wid - 1) = stri_rand_strings(1, wid, '[+]')
  }else{
    best_aligned = alg_all[index.sim]
  }
  
  return (best_aligned)
}


# inFile = "../test_human_mouse_rat/Data_6/Mafft/updated_cds/ENSG00000000460.fa"
# oudir  = "../test_human_mouse_rat/Data_6/Data/Mafft/mapped_cds/"

# seq1 = "AAT===AAACAAAGAATGCTTACTGT---ATAAGGCTTACTGTTCTAGCG===ATCACCGCG===TCATGTCTAGTTATGAACGGC------GGTTTAACATTGAATAGCAAGGCACTTCCA---TAATAGGGCCGTC===GTAATTGTCTAATATAG------ATAGTA==="
# seq2 = "TAA------AA===AATTTGATGCTACATTGGATGAGTCTACTTCGAGCGCGCCGCATCGATTGCAAGAGCAGTGTTGCCT===AAGAGCCGTTAGATGCGTCGTTG---ATCGCGTCCGATAATTCGGGAGTTGTGC===CCCAATATTTAATATGATGA===TAGCTATAA"

################################################
main = function(inFile,ouDir,num1,num2){
  dna = readBStringSet(inFile)
  name = names(dna)
  
  #Set up var (global var)
  num1 = as.numeric(num1)
  num2 = as.numeric(num2)
  Window  <<- num1
  Wall    <<- num2
  
  #String mode
  spec.1   = toString(dna[[1]])
  spec.2   = toString(dna[[2]])
  
  #Capture all gaps of different lengths(3,6,9,12)
  dna.1 = str_split(as.character(dna), '')
  g     = lapply(dna.1, function(x) { IRanges(x == '-')})
  g     = IRangesList(g)
  m     = g[[1]]
  r     = g[[2]]
  
  wid.m = width(m)
  wid.r = width(r)
  l.m   = length(wid.m)
  l.r   = length(wid.r)
  pos.m = start(m)
  pos.r = start(r)
  
  
  if((l.m > 0) & (l.r > 0)){
    for(i in 1:l.m){
      spec.1 = Merge(pos.m[i], wid.m[i], spec.1, spec.2)
    }
    for(i in 1:l.r){
      spec.2 = Merge(pos.r[i], wid.r[i], spec.2, spec.1)
    }
  }else if((l.m > 0) & (l.r == 0)){
    for(i in 1:l.m){
      spec.1 = Merge(pos.m[i], wid.m[i], spec.1, spec.2)
    }
  }else if((l.m == 0) & (l.r > 0)){
    for(i in 1:l.r){
      spec.2 = Merge(pos.r[i], wid.r[i], spec.2, spec.1)
    }
  }else{
    quit()
  }
  
  #Filter out the seq with '+++' only.
  flag = unlist(lapply(c(spec.1, spec.2), function(x){grepl('-', x)}))
  if(any(flag)==TRUE){
    new_seq = list(spec.1, spec.2)
    write.fasta(sequences=new_seq, names=name, nbchar=80,
                open="w", as.string = TRUE, file.out=paste0(ouDir, basename(inFile)) )
  }else{
    quit()
  }
}

#######################################
#compatible with run_test.R
args = commandArgs(trailingOnly=TRUE)
 if(interactive()){
  main(args[1], args[2], args[3], args[4])
 }




