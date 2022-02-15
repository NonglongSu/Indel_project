suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(BiocGenerics))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(S4Vectors))
suppressPackageStartupMessages(library(seqinr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(stringi))
suppressPackageStartupMessages(library(stringr))

# setwd("~/Dropbox (ASU)/Indel_project/Script")

# Find gap range
find_gap_range = function(X){
  X.1   = str_split(as.character(X), '')
  g     = lapply(X.1, function(x) { IRanges(x == '-')})
  g     = IRangesList(g)
  return(g)
}



# Alignment algo (localized- global alignment)
# u,v referencs
# p,q--focal species
# u = "AAT---AAA"
# v = "TAAAAAAAA"
# p = "TAA---TTT"
# q = "TAA---AGA"
# The order of u,v,p,q does not matter!
Align = function(wid.set){
  Tree   = lapply(wid.set, function(x){str_split(x, "")[[1]]})
  cost.T = c()
  for (i in 1:length(Tree[[1]])) {
    #site-by-site comparison
    site = unlist(lapply(Tree, `[[`, i))
    comm = max(as.data.frame(table(site))$Freq)
    diff = length(unique(site))
    
    if((site[1] == site[2]) && (comm == 2) && (diff == 2) ){
      cost = comm - 1 
    }else{
      cost = length(site) - comm
    }
    cost.T = c(cost.T, cost) 
  }
  return(sum(cost.T))
}

#left_sliding
left_slide = function(seq_v, index, wid){
  swap(seq_v[index - 1], seq_v[index + wid - 1])
  seq = paste(seq_v, collapse = "")
  return (seq)
}

#right_sliding
right_slide = function(seq_v, index, wid){
  swap(seq_v[index], seq_v[index + wid])
  seq = paste(seq_v, collapse = "")
  return (seq)
}


#Generate the pseudo seq with an optimal alignment
Merge = function(Arrow, toMatch, idx, wid, seqs){
  
  # Window position
  start = idx - Wall
  stop  = idx + Wall + wid - 1
 
  if(stop > nchar(seqs[1])){
    stop = nchar(seqs[1])
  }
  if(start < 1){
    start = 1
  }
  
  
  #Determine if the gap is insides mouse/rat.  
  if(Arrow == 1){
    seq.res      = seqs[-3]
    best_aligned = seqs[3]
    wid.tar      = substr(seqs[3], start, stop)
  }else{ #Arrow == 0
    seq.res      = seqs[-4]
    best_aligned = seqs[4]
    wid.tar      = substr(seqs[4], start, stop)
  }
  wid.ref = c()
  for (w in 1:length(seq.res)) {
    wid.ref[w] = substr(seq.res[w], start, stop)
  }
  seq.owned = c(seqs[1], seqs[2], best_aligned)
  
  # Prepare for SW methods
  idx.l    = idx.r    = idx
  align.l  = align.r  = best_aligned
  seq.resl = seq.resr = seq.owned
  
  wid.l       = wid.r         = c()
  sim_total.l = sim_total.r   = c()
  wid.set.l = wid.set.r = wid.set.T = list()
  wid.1.l   = wid.1.r   = wid.1.T   = c()    
  
  wid.set.T[[1]] = c(wid.ref[-3], wid.tar)
  
  
  
  # Determine if the event is insertion/deletion
  pattern = c(TRUE, TRUE)    
  
  if(all(toMatch == pattern)){ #insertion.
    for(k in 1:Window){
      for (j in 1:length(seq.resl)){
        seq.resl[j]  = left_slide(str_split(seq.resl[j], "")[[1]], idx.l, wid)
        wid.l[j]     = substr(seq.resl[j], start, stop)
      }
      sim_total.l[k] = Align(c(wid.l, wid.ref[3])) 
      wid.set.l[[k]] = wid.l         
      #all_aligned.l[[k]] = seq.resl
      
      idx.l          = idx.l - 1
    }
    for(k in 1:Window){
      for (j in 1:length(seq.resr)){
        seq.resr[j]  = right_slide(str_split(seq.resr[j], "")[[1]], idx.r, wid)
        wid.r[j]     = substr(seq.resr[j], start, stop)
      }
      sim_total.r[k] = Align(c(wid.r, wid.ref[3])) 
      wid.set.r[[k]] = wid.r   
      #all_aligned.r[[k]] = seq.resr
      
      idx.r          = idx.r + 1
    }
    wid.set.T = append(wid.set.T, c(wid.set.l, wid.set.r))
  }else{ #deletion.
    for(k in 1:Window){
      align.l        = left_slide(str_split(align.l, "")[[1]], idx.l, wid)
      wid.slide      = substr(align.l, start, stop)
      sim_total.l[k] = Align(c(wid.ref, wid.slide)) 
      wid.1.l[k]     = wid.slide
      
      #one_aligned.l[k]  = align.l
      idx.l          = idx.l - 1
    }
    for(k in 1:Window){
      align.r        = right_slide(str_split(align.r, "")[[1]], idx.r, wid)
      wid.slide      = substr(align.r, start, stop )
      sim_total.r[k] = Align(c(wid.ref, wid.slide)) 
      wid.1.r[k]     = wid.slide
      
      #one_aligned.r[k]   = align.r
      idx.r = idx.r + 1
    }
    wid.1.T   = c(wid.tar, wid.1.l, wid.1.r)
  }
  
  # Prepare original alignment score.
  sim_before = Align(c(wid.ref, wid.tar)) 
  sim_all    = c(sim_before, sim_total.l, sim_total.r)
  index.sim  = which(sim_all == min(sim_all))
  
  # Update multiple optimal hits --- ---> +++
  if(all(toMatch == pattern)){
    if(length(index.sim) > 1){
      for (w in 1:length(wid.set.T[[1]])) {
        subseq(seq.owned[w], start = idx, end = idx + wid - 1) = stri_rand_strings(1, wid, '[+]')
      }
    }else{
      best.wid = wid.set.T[[index.sim]]
      for (w in 1:length(best.wid)){
        subseq(seq.owned[w], start = start, end = stop) = best.wid[w]
      }
    }
    new.set = c(seq.owned, seq.res[3])
  }else{
    if(length(index.sim) > 1){
      subseq(best_aligned, start = idx, end = idx + wid - 1) = stri_rand_strings(1, wid, '[+]')
    }else{
      subseq(best_aligned, start = start, end = stop) = wid.1.T[index.sim]
    }
    new.set = c(seq.res[1:2], best_aligned, seq.res[3])
  }
  
  # Keep order the same
  if(Arrow == 0){
    swap(new.set[3], new.set[4]) 
  }
  
  return (new.set)
}


# inFile = "../test_human_mouse_rat/Data_6.2/Mafft/updated_cds/ENSG00000100344.fa"
# ouDir  = "../test_human_mouse_rat/Data_6.2/Mafft/mapped_cds/"

# num1 = 6
# num2 = 12

main = function(inFile, ouDir, num1, num2){
  
  # Read files
  dna = readBStringSet(inFile,format = "fasta")
  name = names(dna)
  
  # Set up vars
  num1 = as.numeric(num1)
  num2 = as.numeric(num2)
  Window  <<- num1
  Wall <<- num2
  
  # String mode
  seq.set = c()
  for (i in 1:length(dna)) {
    seq.set = c(seq.set, toString(dna[[i]]))
  }
  
  
  # Find gap range
  gap = find_gap_range(dna)
  m   = gap[[length(gap) - 1]]          
  r   = gap[[length(gap)]]  
  
  wid.m = width(m)  
  wid.r = width(r)   
  
  l.m   = length(wid.m)
  l.r   = length(wid.r)
  
  pos.m = start(m)
  pos.r = start(r)
  
  h_c       = gap[-((length(gap) - 1):length(gap))]
  toMatch.m = lapply(h_c, function(x){m %in% x})
  toMatch.r = lapply(h_c, function(x){r %in% x})
  
  #seq.set = list(spec.1,spec.2,spec.3,spec.4) 
  
  # Using arrow to indicate where are the gaps from.
  Arrow = 1
  if((l.m > 0) & (l.r > 0)){
    for(i in 1:l.m){
      toMatch.site  = unlist(lapply(toMatch.m, `[[`, i))
      seq.set       = Merge(Arrow, toMatch.site, pos.m[i], wid.m[i], seq.set) 
    }
    Arrow = 0
    for(i in 1:l.r){
      toMatch.site  = unlist(lapply(toMatch.r, `[[`, i))
      seq.set       = Merge(Arrow, toMatch.site, pos.r[i], wid.r[i], seq.set) 
    }
  }else if((l.m > 0) & (l.r == 0)){
    for(i in 1:l.m){
      toMatch.site  = unlist(lapply(toMatch.m, `[[`, i))
      seq.set       = Merge(Arrow, toMatch.site, pos.m[i], wid.m[i], seq.set) 
    }
  }else if((l.m == 0) & (l.r > 0)){
    Arrow = 0
    for(i in 1:l.r){
      toMatch.site  = unlist(lapply(toMatch.r, `[[`, i))
      seq.set       = Merge(Arrow, toMatch.site, pos.r[i], wid.r[i], seq.set) 
    }
  }else{
    quit()
  }
  
  # Filter out the sequence with '+++' only.
  flag = unlist(lapply(seq.set, function(x){grepl('-', x)}))
  if(any(flag) == TRUE){
    Seq.set = as.list(seq.set)
    write.fasta(sequences = Seq.set, names = name, nbchar = 80,
                open = "w", as.string = TRUE, file.out = paste0(ouDir, basename(inFile)) )
  }else{
    quit()
  }  
  
}








# compatible with run_test.R
# if(interactive()){
args = commandArgs(trailingOnly = TRUE)
main(args[1], args[2], args[3], args[4])
#}
