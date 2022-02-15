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
# p,q--focal species
# u,v referencs

# u = "AAT---AAA"
# v = "TAAAAAAAA"
# p = "TAA---TTT"
# q = "TAA---AGA"
###
Align = function(u,v,p,q){
  u = str_split(u,"")[[1]]
  v = str_split(v,"")[[1]]
  p = str_split(p,"")[[1]]
  q = str_split(q,"")[[1]]

  Tree = list(u,v,p,q)
  cost.T = c()
  
  for (i in 1:length(u)) {
    #site-by-site comparison
    site = unlist(lapply(Tree, `[[`, i),use.names = FALSE)
    comm = max(as.data.frame(table(site))$Freq)
    diff = length(unique(site))
    
      if(site[1]==site[2] && comm==2 && diff==2 ){
        cost = comm - 1 
      }else{
        cost = length(site) - comm
      }
    cost.T = c(cost.T,cost) 
  }
  
  return(sum(cost.T))

}

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


#Generate the pseudo seq with an optimal alignment
# z-mouse
# w-rat
###
Merge = function(Arrow,toMatch,idx,wid,seq.set){
  
  x = seq.set[[1]]
  y = seq.set[[2]]
  z = seq.set[[3]]
  w = seq.set[[4]]
  #Window position
  start = idx-wall
  stop  = idx+wall+wid-1
  
  #Determine if the focus is on mouse/rat. 
  if(Arrow == 1){
    wid.ref1 = substr(x, start, stop)
    wid.ref2 = substr(y, start, stop)
    window   = substr(z, start, stop)
    wid.ref3 = substr(w, start, stop)
    best_aligned = z        # mouse
    focalRef     = w
  }else{
    wid.ref1 = substr(x, start, stop)
    wid.ref2 = substr(y, start, stop)
    wid.ref3 = substr(z, start, stop)
    window   = substr(w, start, stop)
    best_aligned = w        # rat
    focalRef     = z
  }
  
  #Prepare for SW methods
  idx.l   = idx.r   = idx
  wid.L   = wid.R   = c()
  
  sim_total.l   = sim_total.r   = c()
  all_aligned.l = all_aligned.r = alg_all = list()
  one_aligned.l = one_aligned.r = alg_one = c()
  
  align.l = align.r = best_aligned
  seq.resl = seq.resr = alg_all[[1]] = c(x,y,best_aligned)
  
  #Determine if the event is insertion/deletion
  pattern = c(TRUE,TRUE)    
  if(all(toMatch==pattern)){
    for(k in 1:wall){
      for (j in 1:length(seq.resl)){
        seq.resl[j] = left_slide(str_split(seq.resl[j],"")[[1]],idx.l,wid)
        wid.L[j]    = substr(seq.resl[j], start, stop)
      }
      wid.L = c(wid.L,wid.ref3)
      sim_total.l[k]     = Align(wid.L[1],wid.L[2],wid.L[3],wid.L[4]) 
      all_aligned.l[[k]] = seq.resl
      idx.l = idx.l-1
    }
    for(k in 1:wall){
      for (j in 1:length(seq.resr)){
        seq.resr[j] = right_slide(str_split(seq.resr[j],"")[[1]],idx.r,wid)
        wid.R[j]   = substr(seq.resr[j], start, stop)
      }
      wid.R = c(wid.R,wid.ref3)
      sim_total.r[k]     = Align(wid.R[1],wid.R[2],wid.R[3],wid.R[4]) 
      all_aligned.r[[k]] = seq.resr
      idx.r = idx.r+1
    }
    alg_all   = append(alg_all,c(all_aligned.l,all_aligned.r))
  }else{
    for(k in 1:wall){
      align.l = left_slide(str_split(align.l,"")[[1]],idx.l,wid)
      wid.tmp = substr(align.l, start, stop)
      sim_total.l[k]     = Align(wid.ref1,wid.ref2,wid.ref3,wid.tmp) 
      one_aligned.l[k]   = align.l
      idx.l =idx.l-1
    }
    for(k in 1:wall){
      align.r = right_slide(str_split(align.r,"")[[1]],idx.r,wid)
      wid.tmp = substr(align.r, start, stop )
      sim_total.r[k]     = Align(wid.ref1,wid.ref2,wid.ref3,wid.tmp) 
      one_aligned.r[k]   = align.r
      idx.r =idx.r+1
    }
    alg_one   = c(best_aligned,one_aligned.l,one_aligned.r)
  }
  
  #Prepare original alignment score.
  sim_ref_public = Align(window,wid.ref1,wid.ref2,wid.ref3) 
  sim_all        = c(sim_ref_public,sim_total.l,sim_total.r)
  index.sim      = which(sim_all==min(sim_all))
  #Update multiple optimal allignment ---/+++
  if(all(toMatch==pattern)){
    if(length(index.sim) > 1){
      for (i in 1:length(alg_all[[1]])) {
        subseq(alg_all[[1]][i],start=idx, end=idx+wid-1) = stri_rand_strings(1,wid,'[+]')
      }
      new.set = c(alg_all[[1]],focalRef)
    }else{
      best_aligned = alg_all[[index.sim]]
      new.set = c(best_aligned,focalRef)
    }
  }else{
    if(length(index.sim) > 1){
      subseq(best_aligned,start=idx, end=idx+wid-1) = stri_rand_strings(1,wid,'[+]')
      new.set = c(x,y,best_aligned,focalRef)
    }else{
      best_aligned = alg_one[index.sim]
      new.set = c(x,y,best_aligned,focalRef)
    }
  }
  #Keep the same order.
  if(Arrow ==0){
    swap(new.set[3],new.set[4]) 
  }
  
  return (new.set)
}


#inFile = "../Raw_data.2.outgroups/test.fa"

main = function(inFile,ouDir,num){
  
  #Read file
  dna = readBStringSet(inFile,format = "fasta")
  name = names(dna)
  
  #set up a window  
  num = as.numeric(num)
  wall <<- num
  
  #String mode
  spec.1 = toString(dna[[1]],width=NULL)                ## human
  spec.2 = toString(dna[[2]],width=NULL)                ## hamster
  spec.3 = toString(dna[[3]],width=NULL)                ## mouse
  spec.4 = toString(dna[[4]],width=NULL)                ## rat
  #Find all gaps
  dna.1 = str_split(as.character(dna), '')
  g = lapply(dna.1, function(x) { IRanges(x == '-')})
  g = IRangesList(g)
  
  m = g[[3]]  
  r = g[[4]]  
  
  wid.m = width(m)
  wid.r = width(r)
  
  l.m = length(wid.m)
  l.r = length(wid.r)
  
  pos.m = start(m)
  pos.r = start(r)
  
  #Prepare for indentifying the ins/del
  h   = g[[1]]
  c   = g[[2]]
  h_c = list(g[[1]],g[[2]])
  
  toMatch.m = lapply(h_c, function(x){m %in% x})
  toMatch.r = lapply(h_c, function(x){r %in% x})
  
  seq.set = list(spec.1,spec.2,spec.3,spec.4)
  Arrow   = 1
  #Update the species set.
  if(l.m>0 & l.r>0){
    for(i in 1:l.m){
      toMatch.site  = unlist(lapply(toMatch.m, `[[`, i), use.names=FALSE)
      seq.set       = Merge(Arrow,toMatch.site,pos.m[i],wid.m[i],seq.set) 
    }
    for(i in 1:l.r){
      Arrow = 0
      toMatch.site  = unlist(lapply(toMatch.r, `[[`, i), use.names=FALSE)
      seq.set       = Merge(Arrow,toMatch.site,pos.r[i],wid.r[i],seq.set) 
    }
  }else if(l.m>0 & l.r==0){
    for(i in 1:l.m){
      toMatch.site  = unlist(lapply(toMatch.m, `[[`, i), use.names=FALSE)
      seq.set       = Merge(Arrow,toMatch.site,pos.m[i],wid.m[i],seq.set) 
    }
  }else if(l.m==0 & l.r>0){
    for(i in 1:l.r){
      Arrow = 0
      toMatch.site  = unlist(lapply(toMatch.r, `[[`, i), use.names=FALSE)
      seq.set       = Merge(Arrow,toMatch.site,pos.r[i],wid.r[i],seq.set) 
    }
  }else{
    quit(save = "default", status = 0, runLast = TRUE)
  }
  
  #Output the fasta files
 #ouDir="../Tmp/"
  newSeq.set = as.list(seq.set)
  write.fasta(sequences = newSeq.set, names = name, nbchar=80,
              open = "w", as.string = TRUE, file.out = paste0(ouDir,basename(inFile)) )
  
  
}

#compatible with run_test.R
#if(interactive()){
args = commandArgs(trailingOnly = TRUE)
main(args[1],args[2],args[3])
#}
