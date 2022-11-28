# Update all gaps from 90 species
suppressWarnings(suppressMessages(library(tidyverse)))
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(BiocGenerics))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(S4Vectors))
suppressPackageStartupMessages(library(seqinr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(stringi))



#setwd("~/Dropbox (ASU)/Indel_project/Script/90")

# Initial word size (Window) as 27(6+3+6). 
# convert ---AAA/AAA--- to ===AAA/AAA===
# convert ---/---       to ===/===
# convert ---AAA---     to ===AAA===

#Remove any large gaps
filter_Long = function(x, X){
  if(length(x) == 0){
    return(X)
  }else{
    st  = start(x)
    en  = end(x)
    wid = width(x) 
    for(i in 1:length(x)){
      subseq(X, start = st[i], end = en[i]) = stri_rand_strings(1, wid[i], '[=]')
    }
  }
  return(X)
}

#Clean up of start/end position of the seq
start_stop_test = function(X){
  Start = substr(X, 1, Window)
  IsGap.1 = grepl('-', Start)
  
  Stop = substr(X, (nchar(X) - (Window - 1)), nchar(X))
  IsGap.2 = grepl('-', Stop) 
  
  if(IsGap.1 == TRUE){
    Start = gsub('-', '=', Start)
    substr(X, 1, Window) = Start
  }
  
  if(IsGap.2 == TRUE){
    Stop = gsub('-', '=', Stop)
    substr(X, (nchar(X) - (Window - 1)), nchar(X)) = Stop
  }
  return(X)
}

#Adjust gap distance
ad_gap_dis = function(y, seq, ref){
  
  if(length(y) == 0){
    return(c(seq, ref))
    }
  
  pos.all = start(y)
  wid.all = width(y) 
  pos.en  = end(y)
  
  for(i in 1:length(pos.all)){
    pos = pos.all[i]
    wid = wid.all[i]
    en  = pos.en[i]
    left   = substr(seq, start = pos - Wall - 3, stop = pos - 1)              ## extra three is to make sure each gap is independent.
    right  = substr(seq, start = pos + wid, stop = pos + wid + (Wall + 2))
    window = substr(ref, start = pos - Wall - 3, stop = pos + wid + (Wall + 2))
    
    l1  = grep("[-=]", left)
    r1  = grep("[-=]", right)
    Wid = grep("[-=]", window)
    
    if(length(l1) != 0 | length(r1) != 0 | length(Wid) != 0){
      subseq(seq, start = pos, end = en) = stri_rand_strings(1, wid, '[=]')
    }
  }
  res = c(seq, ref)
  return(res)
}

#Single fasta analysis
gap_update = function(inFile,ouFile){
    dna   = readDNAStringSet(inFile, format="fasta")
    name  = names(dna)
    spec.1= toString(dna[[1]])
    spec.2= toString(dna[[2]])
  
     #Find gap range
    dna.1 = str_split(as.character(dna),'')
    g     = lapply(dna.1, function(x) { IRanges(x=='-')})
    g     = IRangesList(g)
  
    m     = g[[1]]
    r     = g[[2]]
    wid.m = width(m)
    wid.r = width(r)
  
    #Test gap range
    gap = c(3,6,9,12)
    m.null = m[!(wid.m %in% gap)]
    r.null = r[!(wid.r %in% gap)]
  
    M = filter_Long(m.null,spec.1)
    R = filter_Long(r.null,spec.2)  
  
    #Clean up of start-end region of the seq
    M.1 = start_stop_test(M)
    R.1 = start_stop_test(R)
  
    #Update old gap length 
    M.2 = str_split(as.character(M.1),'')
    R.2 = str_split(as.character(R.1),'')
  
    g.m = lapply(M.2, function(x) { IRanges(x=='-')})
    g.r = lapply(R.2, function(x) { IRanges(x=='-')})
  
    g.m = IRangesList(g.m)[[1]]
    g.r = IRangesList(g.r)[[1]]
  
    #Update the gap via swapping reference
    M_R = ad_gap_dis(g.m, M.1, R.1)
    R_M = ad_gap_dis(g.r,M_R[2],M_R[1])
  
    #Update the focal-seq
    new_seq = list(R_M[2],R_M[1])
    write.fasta(sequences=new_seq, names=name, nbchar=100,
              open="w", as.string=TRUE, file.out=ouFile)
}

####################################################
# Dir  = "../../test_90_species/Raw_data/mapped_cds/07_yeast_aligned_cds"
# ouD  = "../../test_90_species/Data/up_cds/07_yeast_aligned_cds"
# num1 = "6" 
# num2 = "12"

main = function(Dir, ouD, num1, num2){
 
  num1    = as.numeric(num1)
  num2    = as.numeric(num2)
  Window  <<- num1
  Wall    <<- num2
  
  # Read euk- and pro- dirs
  Euk_Pro = list.files(Dir,full.names=TRUE)
  
  for (i in 1:length(Euk_Pro)) {
    ouF = paste0(ouD,"/",basename(Euk_Pro[i]))
    if (file.exists(ouF)){
      next()
    }else{
      gap_update(Euk_Pro[i],ouF)
      print(i)
    }
  }
}


############################################
#compatible with run_test.R
#if(interactive()){
 args = commandArgs(trailingOnly=TRUE)
 main(args[1], args[2], args[3], args[4])
#}


