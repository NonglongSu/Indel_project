# Create a table of alignment score for all of INDELs. 

library(Biostrings)
library(BiocGenerics)
library(parallel)
library(S4Vectors)
library(seqinr)
library(stringr)

#setwd("~/Dropbox (ASU)/Indel_project/test_human_mouse_rat/Script")

Index = function(x){
  res = c()
  posAll = gregexpr(text=x,'[^-]---[^-]')
  if(posAll[[1]][1] != -1){
    res = posAll[[1]][posAll[[1]]>lwall]
  }
  return(res)
}

Align = function(u,v){
  u = str_split(u,"")
  v = str_split(v,"")
  u = u[[1]]
  v = v[[1]]
  id_pos = length(which(u[u == v] != '-'))
  align_pos = wall*2
  sim = 100*(id_pos)/(align_pos)
  return (sim)
}

#left_sliding
left_slide = function(seq,seq_v,index){
  str_sub(seq,start=index,end=index) = '-'
  str_sub(seq,start=(index+3),end=(index+3)) = seq_v[index]
  return (seq)
}

#right_sliding
right_slide = function(seq,seq_v,index){
  str_sub(seq,start=index+1,end=index+1) = seq_v[index+4]
  str_sub(seq,start=(index+4),end=(index+4)) = '-'
  return (seq)
}

#Convert string to vector
str_convert = function(s){
  s1 = str_split(s,"")
  s1 = s1[[1]]
  return(s1)
}

Merge = function(idx,x,y){
  
  #Original window seq
  start = idx-lwall 
  stop  = idx+rwall
  wid.1   = substr(x, start, stop)
  wid_ref = substr(y, start, stop)
  
  #Original alignment score
  sim_ref_public = Align(wid_ref,wid.1) 
  
  #left_aligning
  sim_total.l   = c()
  new_aligned.l = x
  idx.l         = idx
  for(k in 1:wall){
    new_aligned.l = left_slide(new_aligned.l,str_convert(new_aligned.l),idx.l)
    wid.l=substr(x = new_aligned.l, start, stop)
    sim_total.l[k]= Align(wid_ref,wid.l) 
    idx.l =idx.l-1
  }
  
  #right_aligning
  sim_total.r   = c()
  new_aligned.r = x
  idx.r         = idx
  for(k in 1:wall){
    new_aligned.r = right_slide(new_aligned.r,str_convert(new_aligned.r),idx.r)
    wid.r=substr(x = new_aligned.r, start, stop )
    sim_total.r[k] = Align(wid_ref,wid.r) 
    idx.r =idx.r+1
  }
  
  #Recording multiple optimal allignment
  #alg_all= c(all_aligned.l,all_aligned.r)
  
  sim_all = c(sim_ref_public,sim_total.l,sim_total.r)
  
  return(sim_all)
}

scoring = function(idx,x,y){
  ST = c()
  for(i in 1:length(idx)){
    Score = Merge(idx[i],x,y)
    Score = as.double(Score)
    ST = rbind(ST,Score)
  }
  
  return(ST)
}


#inFile="../Data/tmp/test_pseudo_seq.fa"
#inFile = "../Data/updated_cds/ENSG00000000460.fa"

# seq1 = "AAT===AAACAAAGAATGCTTACTGT---ATAAGGCTTACTGTTCTAGCG===ATCACCGCG===TCATGTCTAGTTATGAACGGC------GGTTTAACATTGAATAGCAAGGCACTTCCA---TAATAGGGCCGTC===GTAATTGTCTAATATAG------ATAGTA==="
# seq2 = "TAA------AA===AATTTGATGCTACATTGGATGAGTCTACTTCGAGCGCGCCGCATCGATTGCAAGAGCAGTGTTGCCT===AAGAGCCGTTAGATGCGTCGTTG---ATCGCGTCCGATAATTCGGGAGTTGTGC===CCCAATATTTAATATGATGA===TAGCTATAA"

main = function(dir){
  
  inFile = list.files(dir, full.names = TRUE)
  
  #global var
  window_size <<- 27
  wall <<- (window_size-3)/4
  lwall <<- wall -1
  rwall <<- wall + 3
  
  Tab = c()
  
  for(i in 1:length(inFile)){
    dna   = readBStringSet(inFile[i],format = "fasta")
    mouse = toString(dna[[2]],width=NULL)
    rat   = toString(dna[[3]],width=NULL)
    idx.m = Index(mouse)
    idx.r = Index(rat)
    
    if(is.null(idx.m) && is.null(idx.r)){
        next
    }else if(!is.null(idx.m) && !is.null(idx.r) ){
        s_m = scoring(idx.m,mouse,rat)
        s_r = scoring(idx.r,rat,mouse)
        S   = rbind(s_m,s_r)
    }else if(!is.null(idx.m) && is.null(idx.r) ){
        S = scoring(idx.m,mouse,rat)
    }else{
        S = scoring(idx.r,rat,mouse)
    }
    
    Tab = rbind(Tab,S)
  }
  
  colnames(Tab) <- c(0,-1,-2,-3,-4,-5,-6,1,2,3,4,5,6)
  rownames(Tab) <- c()
  TAB = round(Tab,2)

  write.table(TAB,file = "../Data/aln_socre_report.txt",
              sep = "\t", append = FALSE,col.names = TRUE,row.names = FALSE)
}

#if(interactive()){
args = commandArgs(trailingOnly = TRUE)
main(args[1])
#}

