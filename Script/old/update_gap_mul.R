#Update the fasta files especially for muting the disqualified gaps. 

library(tidyverse)
library(Biostrings)
library(BiocGenerics)
library(parallel)
library(S4Vectors)
library(seqinr)
library(dplyr)
library(stringi)

#setwd("~/Dropbox (ASU)/Indel_project/test_human_mouse_rat/Script")

#Initial word size (Window) as 33(15+3+15)
#convert ---/---       to ===/===
#convert ---AAA---     to ===AAA===

#modify any large gaps
Modify_large = function(x,X){
  st  = start(x)
  en  = end(x)
  wid = width(x) 
  for(i in 1:length(x)){
    subseq(X,start=st[i], end=en[i]) = stri_rand_strings(1,wid[i],'[=]')
  }
  return(X)
}

####
# seq   = "AATATATTTAAA---TAACAAGTAATA---ATGCTTACTGTATAG---"
# seq.1 = str_split(as.character(seq), '')
# y = lapply(seq.1, function(x) { IRanges(x == '-')})
# y = IRangesList(y)
# y = y[[1]]
####

#Adjust close gap-gap in the same species
ad_gap_test = function(seq,y){
  
  if(length(y)==0){
    Res = list(seq,y)
  }else{
   len = length(y)
   pos.start = start(y)
   pos.end   = end(y)
   wid.all   = width(y) 
   #at least 3 extra nts to separate windows.
   filter = Wall + 3
   Record = c()
   
   for(i in 1:len){
     pos = pos.start[i]
     wid = wid.all[i]
     en  = pos.end[i]
     left  = substr(seq, start = pos-filter, stop = pos-1)
     right = substr(seq, start = pos+wid, stop = pos+wid+(filter-1))
     l1 = grep("[-=]", left)
     r1 = grep("[-=]", right)
    
     if(length(l1)!=0 || length(r1)!=0){
       subseq(seq,start=pos, end=en) = stri_rand_strings(1,wid,'[=]')
       Record = c(Record,i)
     }
   }
   if(length(Record)>0){y = y[-Record]}        #update the g.
   Res = list(seq,y)
  }
  return(Res)
}

######
# tmp="../Raw_data.2.outgroups/test.fa"
# dna.tmp   = readDNAStringSet(tmp,format = "fasta")
# dna.tmp.1 = str_split(as.character(dna.tmp), '')
# g.tmp = lapply(dna.tmp.1, function(x) { IRanges(x == '-')})
# g.tmp = IRangesList(g.tmp)
# 
# tmp.1 = toString(dna.tmp[[1]],width=NULL)
# tmp.2 = toString(dna.tmp[[2]],width=NULL)
# tmp.3 = toString(dna.tmp[[3]],width=NULL)
# tmp.4 = toString(dna.tmp[[4]],width=NULL)
# 
# Z = list(tmp.1,tmp.2,tmp.3,tmp.4)
# G = g.tmp
######

#Deal with special scenarios (---
#                             ------)
overlap_gap_test = function(Z,G){
  
  M.sRest = list(start(G[[1]]),start(G[[2]]),start(G[[4]]))
  M.eRest = list(end(G[[1]]),end(G[[2]]),end(G[[4]]))
  
  R.sRest = list(start(G[[1]]),start(G[[2]]),start(G[[3]]))
  R.eRest = list(end(G[[1]]),end(G[[2]]),end(G[[3]]))
  
  startMatch.M = lapply(M.sRest, function(x){start(G[[3]]) %in% x})
  endMatch.M   = lapply(M.eRest, function(x){end(G[[3]]) %in% x})
  
  startMatch.R = lapply(R.sRest, function(x){start(G[[4]]) %in% x})
  endMatch.R   = lapply(R.eRest, function(x){end(G[[4]]) %in% x})
  #Record the order of "bad" gaps
  Record.m = c()
  Record.r = c()
  if(length(G[[3]])>0){
    for (i in 1:length(M.sRest)) {
      for (j in 1:length(G[[3]])) {
        if(startMatch.M[[i]][j] != endMatch.M[[i]][j]){
          subseq(Z[[3]],start=start(G[[3]])[j], end=end(G[[3]])[j]) = stri_rand_strings(1,width(G[[3]])[j],'[=]')
          Record.m = c(Record.m,j)
        }
      }
    }
  }
  if(length(G[[4]])>0){
  for (i in 1:length(R.sRest)) {
    for (j in 1:length(G[[4]])) {
      if(startMatch.R[[i]][j] != endMatch.R[[i]][j]){
        subseq(Z[[4]],start=start(G[[4]])[j], end=end(G[[4]])[j]) = stri_rand_strings(1,width(G[[4]])[j],'[=]')
        Record.r = c(Record.r,j)
       }
      }
    }
  }
  if(length(Record.m)>0){G[[3]] = G[[3]][-Record.m]}
  if(length(Record.r)>0){G[[4]] = G[[4]][-Record.r]}
  
  Res = list(Z,G)
  return(Res)
}

#Deal with gaps based on maximum parsimony
phylo_gap_test = function(Z,G){
  
  H.start = start(G[[1]])
  C.start = start(G[[2]])
  M.start = start(G[[3]])
  R.start = start(G[[4]])
  H.end = end(G[[1]])
  C.end = end(G[[2]])
  M.end = end(G[[3]])
  R.end = end(G[[4]])
  H.wid = width(G[[1]])
  C.wid = width(G[[2]])
  M.wid = width(G[[3]])
  R.wid = width(G[[4]])
  
  M.rest = list(G[[1]],G[[2]],G[[4]])
  R.rest = list(G[[1]],G[[2]],G[[3]])
  #Use focal-species gaps to match the rest of group.
  pattern.1 = c(TRUE,TRUE,FALSE)
  pattern.2 = c(FALSE,FALSE,FALSE)
  #match gaps
  toMatch.1 = lapply(M.rest, function(x){G[[3]] %in% x})
  toMatch.2 = lapply(R.rest, function(x){G[[4]] %in% x})
  spec.M.rest = Z[-3]
  spec.R.rest = Z[-4]

  #Compare mouse with other species
  if(length(G[[3]])>0){
    for (i in 1:length(G[[3]])) {
      tmp.ls.1   = lapply(toMatch.1, `[[`, i)
      tmp.vec.1  = unlist(tmp.ls.1, use.names=FALSE)
      #Check if the phylogeny tree is correct.
      if(all(tmp.vec.1==pattern.1) || all(tmp.vec.1==pattern.2)){
        for (j in 1:length(spec.M.rest)) {#Check all rest groups
          lf = substr(spec.M.rest[j], start = M.start[i]-Wall-3,stop = M.start[i]-1)
          rt = substr(spec.M.rest[j], start = M.start[i]+M.wid[i],stop = M.start[i]+M.wid[i]+(Wall+2))
          l2 = grep("[-=]", lf)
          r2 = grep("[-=]", rt)
          #Check if they are independent gaps.
          if(length(l2)!=0 || length(r2)!=0){
            for (w in 1:length(Z)) {#Do not use subseq()=stri_rand_strings() because that will replace everything even there is no gaps.
              test.gap = substr(Z[[w]],M.start[i],M.end[i])
              test.gap = gsub('-','=',test.gap)
              substr(Z[[w]],M.start[i],M.end[i]) = test.gap
            }
            break
          }
        }
      }else{#Bad phylogeny tree also be updated.
        for (w in 1:length(Z)) {
          test.gap = substr(Z[[w]],M.start[i],M.end[i])
          test.gap = gsub('-','=',test.gap)
          substr(Z[[w]],M.start[i],M.end[i]) = test.gap
        }
      }
    }
  }
  #Compare rat with other species
  if(length(G[[4]])>0){
    for (i in 1:length(G[[4]])) {
      tmp.ls.2   = lapply(toMatch.2, `[[`, i)
      tmp.vec.2  = unlist(tmp.ls.2, use.names=FALSE)
      #Check if the phylogeny tree is correct.
      if(all(tmp.vec.2==pattern.1) || all(tmp.vec.2==pattern.2)){
        for (j in 1:length(spec.R.rest)) {
          lf = substr(spec.R.rest[j], start = R.start[i]-Wall-3,stop = R.start[i]-1)
          rt = substr(spec.R.rest[j], start = R.start[i]+R.wid[i],stop = R.start[i]+R.wid[i]+(Wall+2))
          l2 = grep("[-=]", lf)
          r2 = grep("[-=]", rt)
          #Check if they are independent gaps.
          if(length(l2)!=0 | length(r2)!=0){
            for (w in 1:length(Z)) {
              test.gap = substr(Z[[w]],R.start[i],R.end[i])
              test.gap = gsub('-','=',test.gap)
              substr(Z[[w]],R.start[i],R.end[i]) = test.gap
            }
            break
          }
        }
      }else{
        for (w in 1:length(Z)) {
          test.gap = substr(Z[[w]],R.start[i],R.end[i])
          test.gap = gsub('-','=',test.gap)
          substr(Z[[w]],R.start[i],R.end[i]) = test.gap
        }
      }
    }
  }
#Update the remained free-gaps within human/hamster species 
H_C.rest  = list(G[[3]],G[[4]])
 toMatch.3 = lapply(H_C.rest, function(x){G[[1]] %in% x})
 toMatch.4 = lapply(H_C.rest, function(x){G[[2]] %in% x})
  
  if (length(G[[1]]) > 0){
    for (m in 1:length(G[[1]])) {
      tmp.ls.3   = lapply(toMatch.3, `[[`, m)
      tmp.vec.3  = unlist(tmp.ls.3, use.names=FALSE)
      if(all(tmp.vec.3==FALSE)){#update all outgroup gaps. 
        subseq(Z[[1]],start=H.start[m], end=H.end[m]) = stri_rand_strings(1,H.wid[m],'[=]')
      }
    }
  }
  if(length(G[[2]])>0){
    for (n in 1:length(G[[2]])){
      tmp.ls.4   = lapply(toMatch.4, `[[`, n)
      tmp.vec.4  = unlist(tmp.ls.4, use.names=FALSE)
      if(all(tmp.vec.4==FALSE)){
        subseq(Z[[2]],start=C.start[n], end=C.end[n]) = stri_rand_strings(1,C.wid[n],'[=]')
      }
    }
  }

  return(Z)
}


#inFile = "../Raw_data.2.outgroups/mapped_cds_mafft_Tmp/ENSG00000117335.fa"
#nFile = "../Tmp1/version3.fa"
main = function(inFile,ouDir,num){
  
  #oudir="../Data/updated_cds/"
  dna = readBStringSet(inFile,format = "fasta")
  name = names(dna)
  dna.l = width(dna)
  #Variables setting
  num = as.numeric(num)
  Window  <<- num
  Wall <<- 2*Window
  
  #String mode
  spec.1 = toString(dna[[1]],width=NULL)
  spec.2 = toString(dna[[2]],width=NULL)
  spec.3 = toString(dna[[3]],width=NULL)
  spec.4 = toString(dna[[4]],width=NULL)
  #Find gap range
  dna.1 = str_split(as.character(dna), '')
  g = lapply(dna.1, function(x) { IRanges(x == '-')})
  g = IRangesList(g)
  
  h = g[[1]]
  c = g[[2]]
  m = g[[3]]
  r = g[[4]]
  wid.h = width(h)
  wid.c = width(c)
  wid.m = width(m)
  wid.r = width(r)
  #Define gap range
  gap = c(3,6,9,12)
  
  h.null = h[!(wid.h %in% gap)]
  c.null = c[!(wid.c %in% gap)]
  m.null = m[!(wid.m %in% gap)]
  r.null = r[!(wid.r %in% gap)]
  h.null.l = length(h.null)
  c.null.l = length(c.null)
  m.null.l = length(m.null)
  r.null.l = length(r.null)
  
  #Remove large gaps
  if(h.null.l>0){
    spec.1 = Modify_large(h.null,spec.1)
  }
  if(c.null.l>0){
    spec.2 = Modify_large(c.null,spec.2)  
  }
  if(m.null.l>0){
    spec.3 = Modify_large(m.null,spec.3) 
  }
  if(r.null.l>0){
    spec.4 = Modify_large(r.null,spec.4) 
  }
  
  #Keep the legit gaps.
  g.1 = g[width(g) %in% gap] 
  
  #clean-up of start-end positioned gaps
  spec  =c(spec.1,spec.2,spec.3,spec.4)
  for (i in 1:length(g.1)) {
    if(length(g.1[[i]]) > 0){
      Start = start(g.1[[i]])
      Stop  = end(g.1[[i]])
      Width = width(g.1[[i]])
      for(j in 1:length(g.1[[i]])){
        if(Start[j] < Window || Stop[j] > (dna.l[i]-Window+1) ){
          subseq(spec[i],start=Start[j], end=Stop[j]) = stri_rand_strings(1,Width[j],'[=]')
          g.1[[i]] = g.1[[i]][-j]
        }
      }
    }
  }
  
  #clean-up of close gaps in the same species
  spec.list = c()
  g.2 = list()
  for (k in 1:length(spec)) {
    spec.list[k] = ad_gap_test(spec[k],g.1[[k]])[1]
    g.2[k]       = ad_gap_test(spec[k],g.1[[k]])[2]
  }
  
  #clean-up of the overlap gaps in different species
  spec.list.1 = overlap_gap_test(spec.list,g.2)[[1]]
  g.3 = overlap_gap_test(spec.list,g.2)[[2]]
  
  #clean-up of phylogeny-aware indels
  spec.list.2 = phylo_gap_test(spec.list.1,g.3) 
  #Write to file
  write.fasta(sequences = spec.list.2,names = name, nbchar=80,
              open = "w",as.string = TRUE, file.out = paste0(ouDir,basename(inFile)) )
  
}

#compatible with run_test.R
#if(interactive()){
args = commandArgs(trailingOnly = TRUE)
main(args[1],args[2],args[3])
#}

