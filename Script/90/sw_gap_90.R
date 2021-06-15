# Slide the gaps and find the best hit. 
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(BiocGenerics))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(S4Vectors))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(seqinr))
suppressPackageStartupMessages(library(stringi))
suppressPackageStartupMessages(library(stringr))

# setwd("~/Dropbox (ASU)/Indel_project/Script/90")

# Hamming distance
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

# left_sliding
left_slide = function(seq_v, index, wid){
  swap(seq_v[index - 1], seq_v[index + wid - 1])
  seq   = paste(seq_v, collapse = "")
  return (seq)
}

# right_sliding
right_slide = function(seq_v, index, wid){
  swap(seq_v[index], seq_v[index + wid])
  seq   = paste(seq_v, collapse = "")
  return (seq)
}

# Convert string to vector
str_convert = function(s){
  s1 = str_split(s, "")
  s1 = s1[[1]]
  return(s1)
}

# Generate the pseudo seq with an optimal alignment
Merge = function(idx, wid, x, y){
  #assume current state is the optimal state
  best_aligned = x             
  
  #original window seq
  start = idx - Wall
  stop  = idx + Wall + wid - 1
  wid.1   = substr(x, start, stop)
  wid_ref = substr(y, start, stop)
  
  #original alignment score
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

# Single fasta analysis
sw_update = function(inFile, ouFile){
  if(file.exists(ouFile)){
    next()
  }else{
   print(basename(inFile))
   dna    = readBStringSet(inFile, format = "fasta")
   name   = names(dna)
   spec.1 = toString(dna[[1]])
   spec.2 = toString(dna[[2]])
  
   #capture all gaps of different lengths(3, 6, 9, 12)
   dna.1 = str_split(as.character(dna), '')
   g = lapply(dna.1, function(x) { IRanges(x == '-')})
   g = IRangesList(g)
  
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
  }else{}
  
  #filter out the seq with just '+++'. 
  flag = unlist(lapply(c(spec.1, spec.2), function(x){grepl('-', x)}))
  if(any(flag) == TRUE){
    new_seq = list(spec.1, spec.2)
    write.fasta(sequences = new_seq, names = name, nbchar = 80,
                open = "w", as.string = TRUE, file.out = ouFile)
  }
 }
  
}

#############################################

# Dir   = "../../test_90_species/Data/up_cds"
# ouD   = "../../test_90_species/Data/sw_cds"
# num1  = "6" ; num2 = "12"

main = function(Dir, ouD, num1, num2){

  # Set up a window (global var)
  num1    = as.numeric(num1)
  num2    = as.numeric(num2)
  Window  <<- num1
  Wall    <<- num2
  
  # Read euk- and pro- dirs
  Euk_Pro = list.files(Dir, full.names = TRUE)
  
  for (i in 1:length(Euk_Pro)) {
    dir.create(file.path(ouD, basename(Euk_Pro[i])))
    Each.sp  = list.files(Euk_Pro[i], full.names = TRUE)
    subD = paste0(ouD, "/", basename(Euk_Pro[i]))
    for (j in 1:length(Each.sp)) {
      output = paste0(subD, "/", basename(Each.sp[j]))
      if (file.exists(output)){
        next()
      }else{
        sw_update(Each.sp[j], output)
      }
    }
  }
  
}

#compatible with run_test.R
 #if(interactive()){
 args = commandArgs(trailingOnly = TRUE)
 main(args[1], args[2], args[3], args[4])
 #}

################


