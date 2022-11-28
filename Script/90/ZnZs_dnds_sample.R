#Generate a csv/tsv file of 
# - typeN
# - typeS
# - Nd
# - N
# - Sd
# - S
# - length of seqA
# - length of seqB 
suppressWarnings(suppressMessages(library(tidyverse)))
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(seqinr))
suppressPackageStartupMessages(library(jsonlite))
suppressPackageStartupMessages(library(profvis))
suppressPackageStartupMessages(library(stringi))


#setwd("~/Dropbox (ASU)/Indel_project/test_90_species")

#Determine number of non-syn/syn sites
cal_sites = function(x,tag){
  unit.x  = substr(x, tag, tag+2)
  s.site  = 0
  
  if(grepl(pat,unit.x)){
    return(c(0, 0))
  }else{
    section = syn[unit.x][[1]]
    for (k in 1:3) {
      Base     = DNA_BASES[substr(unit.x, k, k) != DNA_BASES]
      unit.sub = unit.x
      for (w in 1:length(Base)) {
        substr(unit.sub, k, k) = Base[w]
        if(unit.sub %in% section){
          s.site = s.site + 1/3
        }
      }
    }
  }
  n.site = 3 - s.site
  res = c(n.site, s.site)
  return(res)
}

#Determine number of non-syn/syn substitutions
cal_subs = function(x,y,rdnas,tag){
  unit.x  = substr(x,tag,tag+2)
  unit.y  = substr(y,tag,tag+2)
  sub.s   = 0   
  sub.n   = 0  
  
  if(any(grepl(pat,c(unit.x,unit.y)))){
    return(c(0,0))
  }else{
    section = syn[unit.x][[1]]
    unix    = rdnas[[1]][tag:(tag+2)]
    uniy    = rdnas[[2]][tag:(tag+2)]
    mut.pos = which(unix != uniy) 
    mut.num = length(mut.pos)
    
    if(mut.num == 1){#ATG/TTG
      if(unit.y %in% section){
        sub.s = sub.s + 1
      }else{
        sub.n = sub.n + 1
      }
    }else if(mut.num == 2){#ATG/TCG
      unit.list = collect_allSubs(unit.x, unit.y, mut.pos)
      sub.rate  = cal_subsub(unit.list, sub.s, sub.n, section)
      sub.s = sub.rate[1] / 2
      sub.n = sub.rate[2] / 2
    }else if(mut.num == 3){#ATG/CCT
      unit.list = collect_allSubs(unit.x, unit.y, mut.pos)
      sub.rate  = cal_subsub(unit.list, sub.s, sub.n, section)
      sub.s = sub.rate[1] / 6
      sub.n = sub.rate[2] / 6
    }else{#no substitution ATG/ATG
      return(c(0,0))
    }
  }
  return(c(sub.s,sub.n))
}

# Find all possible substitution pathways.        
# unit.x = "TTT"
# unit.y = "TAC"
# unit.y = "GAC"
# mut.pos = 2:3
# mut.pos = 1:3
collect_allSubs = function(unit.x,unit.y,mut.pos){
  unit.set  = c()
  for (i in 1:length(mut.pos)) {
    unit.sub = unit.x
    substr(unit.sub, mut.pos[i], mut.pos[i]) = substr(unit.y, mut.pos[i], mut.pos[i])
    unit.set = c(unit.set, unit.sub)
    if (length(mut.pos) == 3) {
      for (j in mut.pos[mut.pos != mut.pos[i]]) {
        unit.sb               = unit.sub
        substr(unit.sb, j, j) = substr(unit.y, j, j)
        unit.set              = c(unit.set, unit.sb)
      }
    }
  }
  unit.list = list()
  if(length(unit.set) == 2){
    for (k in 1:length(unit.set)) {
      unit.list = c(unit.list, list(c(unit.x, unit.set[k], unit.y)))
    }
  }else{
    for (k in 1:length(unit.set)) {
      if (k %% 3 == 1) {
        for (w in (k + 1):(k + 2)) {
          unit.list = c(unit.list, list(c(unit.x, unit.set[k], unit.set[w], unit.y)))
        }
      }
    }
  }
  return(unit.list)
}

# Sum the number of substitutions of all possible pathways
# section = c("TTT","TTC")
# 12 , 6
cal_subsub = function(unit.list,sub.s,sub.n,sec){
  for (m in 1:length(unit.list)) {
    sec.var = sec
    for (n in 2:(length(unit.list[[m]])) ) {
      if(unit.list[[m]][n] %in% sec.var){
        sub.s = sub.s + 1
      }else{
        sub.n = sub.n + 1
      }
      sec.var = syn[unit.list[[m]][n]][[1]]
    }
  }
  return(c(sub.s, sub.n))
}


#Remove all gap-positioned string and return a sub-flag
subflag = function(aln2,ngroups){
  flag = rep(0,ngroups)
  R1   = c()
  R2   = c()
  for(k in 1:ngroups){
    dna   = DNAStringSet(c(aln2$Seq1[k],aln2$Seq2[k]))
    g     = IRangesList(lapply(str_split(dna,''), function(x){IRanges(x=='-')}))
    g     = unlist(g)
    r.dna = stri_sub_replace_all(dna,from=sort(start(g)),to=sort(end(g)),replacement='')
    if(r.dna[1]!=r.dna[2]){
      flag[k]=1
    }
    R1[k] = r.dna[1]
    R2[k] = r.dna[2]
  }
  res = list(flag,R1,R2)
  return(res)
}


###############################################
#get the number of typeN, typeS indels (<=12nts)
get_typeNS_sam = function(dna){
  dna.str  = str_split(dna,'')
  seqA     = dna.str[[1]]
  seqB     = dna.str[[2]]
  g        = IRangesList(lapply(dna.str, function(x){IRanges(x=='-')}))
  
  g1 = g[[1]][which(width(g[[1]])<=12)] 
  g2 = g[[2]][which(width(g[[2]])<=12)] 
  
  Sm = Nm = Sr = Nr = 0
  
  if(length(g1)+length(g2) == 0){
    res = c(0,0)
  }else{
    if(length(g1)>0){#gaps in mouse
      subs = cal_zed_subs(seqA,seqB,g1) 
      Nm   = Nm + subs[1]
      Sm   = Sm + subs[2]
      
    }
    if(length(g2)>0){#gaps in rat
      subs = cal_zed_subs(seqB,seqA,g2)
      Nr   = Nr + subs[1]
      Sr   = Sr + subs[2]
      
    }
    s.indel = Sm + Sr 
    n.indel = Nm + Nr
    res     = c(s.indel,n.indel)
  }
  
  return(res)
}

#Number of observed mutational indel events
#sa-gap
#sb-ref
cal_zed_subs = function(sa,sb,gz){
  imat = matrix(0,3,2)
  ps1  = start(gz)
  ps2  = end(gz)
  rem  = ps1 %% 3
  
  #1A-- -AA AAA  2A-- -AA AAA   4AAA AA- --A  3AAA AAA AAA  >phase1
  #1AAA A-- -AA  2AAA AA- --A   4A-- -AA AAA  3A-- -A- --A
  
  #5AA- --A AAA   >phase2
  #5AAA AA- --A  
  
  #6--- NNN NNN NNN --- END  7--- NNN  8NNN NNN  9NNN --- end  10NNN NNN end >phase0    
  #6NNN --- NNN --- NNN END  7NNN NNN  8--- NNN  9NNN NNN end  10NNN --- end             
  
  for (i in 1:length(rem)) {
    if(rem[i] == 2){#phase 1 
      if(sa[ps2[i]+2]=='-'){#3
        tmp.0=sb[ps1[i]-1]
        tmp.1=sb[ps2[i]+1]
        j = ps2[i]+2
        repeat{
          j     = j+1
          if(sa[j]!='-'){
            break
          }
        }
        tmp.2  = sb[j]
        ref    = paste0(tmp.0,tmp.1,tmp.2,            collapse="")
        sub1   = paste0(tmp.0,sb[ps1[i]],sb[ps1[i]+1],collapse="")
        sub2   = paste0(sb[ps2[i]],tmp.1,sb[ps2[i]+2],collapse="")
        sub3   = paste0(sb[j-2],sb[j-1],tmp.2        ,collapse="")
        sec    = syn[[ref]]
        if( (sub1 %in% sec) || (sub2 %in% sec) || (sub3 %in% sec)){ 
          imat[2,2]=imat[2,2]+1 
        }else{
          imat[2,1]=imat[2,1]+1 
        }
      }else{#no-3
        if(sb[ps2[i]+1]=='-'){#1
          j = ps2[i]+1
          repeat{
            j      = j+1
            tmp.1  = sb[j]
            if(tmp.1!='-'){
              break
            }
          }
          tmp.0=sb[ps1[i]-1]
          tmp.2=sb[j+1]
        }else if(sb[ps2[i]+2]=='-'){#2,4
          tmp.0=sb[ps1[i]-1]
          tmp.1=sb[ps2[i]+1]
          j = ps2[i]+2
          repeat{
            j     = j+1
            tmp.2 = sb[j]
            if(tmp.2!='-'){
              break
            }
          }
        }else if(sb[ps1[i]-1]=='-'){#1
          j = ps1[i]-1
          repeat{
            j      = j-1
            tmp.0  = sb[j]
            if(tmp.0!='-'){
              break
            }
          }
          tmp.1 = sb[ps2[i]+1]
          tmp.2 = sb[ps2[i]+2]
        }else{#normal
          tmp.0 = sb[ps1[i]-1]
          tmp.1 = sb[ps2[i]+1]
          tmp.2 = sb[ps2[i]+2]
        }
        ref    = paste0(tmp.0,tmp.1,tmp.2,            collapse="")
        sub1   = paste0(tmp.0,sb[ps1[i]],sb[ps1[i]+1],collapse="")
        sub2   = paste0(sb[ps2[i]],tmp.1,tmp.2,      collapse="")
        sec    = syn[[ref]]
        if(sub1 %in% sec || sub2 %in% sec){ 
          imat[2,2]=imat[2,2]+1 
        }else{
          imat[2,1]=imat[2,1]+1 
        }
      }
    }else if(rem[i] == 0){#phase2
      if(sa[ps1[i]-2]=='-'){#3
        j = ps1[i]-2
        repeat{
          j = j-1
          if(sa[j]!='-'){
            break
          }
        }
        tmp.0  = sb[j]
        tmp.1  = sb[ps1[i]-1]
        tmp.2  = sb[ps2[i]+1]
        ref    = paste0(tmp.0,tmp.1,tmp.2,            collapse="")
        sub1   = paste0(tmp.0,sb[j+1],sb[j+2],        collapse="")
        sub2   = paste0(sb[ps1[i]-2],tmp.1,sb[ps1[i]],collapse="")
        sub3   = paste0(sb[ps2[i]-1],sb[ps2[i]],tmp.2,collapse="")
        sec    = syn[[ref]]
        if((sub1 %in% sec) || (sub2 %in% sec) || (sub3 %in% sec)){ 
          imat[3,2]=imat[3,2]+1 
        }else{
          imat[3,1]=imat[3,1]+1 
        }
      }else{#no-3
        if(sb[ps1[i]-2]=='-'){#2,4
          j = ps1[i]-2
          repeat{
            j      = j-1
            tmp.0  = sb[j]
            if(tmp.0!='-'){
              break
            }
          }
          tmp.1 = sb[ps1[i]-1]
          tmp.2 = sb[ps2[i]+1]
        }else if(sb[ps2[i]+1]=='-'){#5
          j = ps2[i]+1
          repeat{
            j      = j+1
            tmp.2  = sb[j]
            if(tmp.2!='-'){
              break
            }
          }
          tmp.0=sb[ps1[i]-2]
          tmp.1=sb[ps1[i]-1]
        }else if(sb[ps1[i]-1]=='-'){#5
          j = ps1[i]-1
          repeat{
            j      = j-1
            tmp.1  = sb[j]
            if(tmp.1!='-'){
              break
            }
          }
          tmp.0 = sb[j-1]
          tmp.1 = sb[j]
          tmp.2 = sb[ps2[i]+1]
        }else{
          tmp.0 = sb[ps1[i]-2]
          tmp.1 = sb[ps1[i]-1]
          tmp.2 = sb[ps2[i]+1]
        }
        ref    = paste0(tmp.0,tmp.1,tmp.2,            collapse="")
        sub1   = paste0(tmp.0,tmp.1,sb[ps1[i]],       collapse="")
        sub2   = paste0(sb[ps2[i]-1],sb[ps2[i]],tmp.2,collapse="")
        sec    = syn[[ref]]
        if(sub1 %in% sec || sub2 %in% sec){ 
          imat[3,2]=imat[3,2]+1  
        }else{
          imat[3,1]=imat[3,1]+1  
        }
      }
    }else{#phase0
      imat[1,2]=imat[1,2]+1   
    }
  }
  res = colSums(imat)
  return(res)
}


#inD = "Raw_data/align_sampling/01_FcaCaf_aligned_cds"
#ouF = "Results/ZD_sum/01_FcaCaf_aligned_cds.sample.tsv"
main = function(inD,ouF){
  
  source(paste0("../Script/sources/codon_call.R"))
  #construct codons and its degeneracy
  co.res  = codon_call()
  syn     <<- co.res[[3]]
  
  #pattern to be ignored
  toMatch =  c('-','\\+','=','N')
  pat     <<-paste0(toMatch, collapse="|")
  
  #read files
  Files = list.files(inD, full.names=TRUE)
  n     = length(Files)
  datM  = matrix(0,n,8)
  colnames(datM) = c('Zn','Zs','Nd','N','Sd','S','lenA','lenB') 
  
  for (i in 1:n) {
    ##init 0
    n.m  = 0 #number of non-synonymous sites (mouse)
    s.m  = 0 #number of synonymous sites (mouse)
    n.r  = 0 #number of non-synonymous sites (rat)
    s.r  = 0 #number of synonymous sites (rat)
    
    ##read json
    dna.json = fromJSON(Files[i])
    
    ##cal the expected syn-/non-subs only once
    dna0  = DNAStringSet(c(dna.json$aln[2,1],dna.json$aln[2,2]))
    len0  = width(dna0)[1]
    M0    = toupper(toString(dna0[1])) 
    R0    = toupper(toString(dna0[2])) 
    
    lenA  = nchar(str_remove_all(M0,'-'))
    lenB  = nchar(str_remove_all(R0,'-'))
    
    j0 = 1                
    while (j0<len0){
      nSite_M  = cal_sites(M0,j0)
      n.m      = n.m + nSite_M[1]
      s.m      = s.m + nSite_M[2]
      
      nSite_R  = cal_sites(R0,j0) 
      n.r      = n.r + nSite_R[1]
      s.r      = s.r + nSite_R[2]
      j0 = j0+3                   
    }
    non.num = 0.5*(n.m+n.r)
    syn.num = 0.5*(s.m+s.r)
    
    ##cal the observed syn-/non-subs 100 times
    Dat.tmp  = dna.json %>% tidyr::unpack(aln)
    colnames(Dat.tmp)[1:2] = c('Seq1','Seq2')
    aln      = Dat.tmp %>% dplyr::group_by(Seq1,Seq2)
    aln2     = aln %>% dplyr::group_keys(Seq1,Seq2,weight,log_weight)
    gpsize   = group_size(aln)
    ngroups  = n_groups(aln)
    
    ##normalized weight
    minw     = min(aln2$log_weight)
    norm.wei = exp(aln2$log_weight-minw)/sum(exp(aln2$log_weight-minw)*gpsize)
    #print(sum(norm.wei*gpsize))
    
    ##sub-flag
    Flags = subflag(aln2,ngroups)
    Flag  = Flags[[1]]
    rseq1 = Flags[[2]]
    rseq2 = Flags[[3]]
    aln2  = cbind(aln2,rseq1,rseq2)
    
    Sm    = rep(0,ngroups) 
    Nm    = rep(0,ngroups)
    tyn   = rep(0,ngroups) 
    tys   = rep(0,ngroups)
    
    #profvis({
    for(k in 1:ngroups){
      dna   = DNAStringSet(c(aln2$Seq1[k],aln2$Seq2[k]))
      len   = width(dna)[1]
      
      #typeN/typeS
      NS.value  = get_typeNS_sam(dna)  
      tys[k]    = NS.value[1]
      tyn[k]    = NS.value[2]
      
      if(Flag[k]==0){#no subs
        next
      }else{
        rdna  = DNAStringSet(c(aln2$rseq1[k],aln2$rseq2[k]))
        rlen  = width(rdna)[1]
        rdnas = str_split(rdna,'')
        M     = toupper(toString(rdna[1])) 
        R     = toupper(toString(rdna[2])) 
        j = 1                 
        while (j<rlen) {
          Subs = cal_subs(M,R,rdnas,j) 
          Sm[k]= Sm[k] + Subs[1]
          Nm[k]= Nm[k] + Subs[2]
          j = j+3                   
        }
      }
    }
    #})
    
    #add Nd,Sd
    Sd = sum(Sm*gpsize*norm.wei)
    Nd = sum(Nm*gpsize*norm.wei)
    
    #add typeN, typeS
    tyS = sum(tys*gpsize*norm.wei)
    tyN = sum(tyn*gpsize*norm.wei)
    
    datM[i,] = c(tyN,tyS,Nd,non.num,Sd,syn.num,lenA,lenB)
    print(i)
  }
  
  #Output int tsv file
  write.table(datM, file=ouF,sep='\t', append=FALSE, quote=FALSE, row.names=FALSE, col.names=TRUE)
}


###########################################
args = commandArgs(trailingOnly=TRUE)
main(args[1],args[2])