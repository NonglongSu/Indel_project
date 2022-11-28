#Generate a csv/tsv file of 
# - typeN
# - typeS
# - Nd
# - N
# - Sd
# - S
# - length of mouse seq
# - length of rat seq 
suppressWarnings(suppressMessages(library(dplyr)))
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(seqinr))


#setwd("~/Dropbox (ASU)/Indel_project/test_90_species")

########################################################
#Determine number of non-syn/syn sites
cal_sites = function(x,tag){
  
  unit.x  = substr(x, tag, tag+2)
  unit.x  = toupper(unit.x)
  toMatch = c('-','\\+','=','N')
  s.site  = 0
  
  if(grepl(paste0(toMatch, collapse="|"), unit.x)){#independence 
    return(c(0, 0))
  }else{
    section = codon[[which(sapply(codon, function(X){unit.x %in% X}))]]
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
cal_subs = function(x, y, tag){
  
  unit.x  = substr(x, tag, tag + 2)
  unit.y  = substr(y, tag, tag + 2)
  sub.s   = 0   
  sub.n   = 0  
  toMatch = c('-', '\\+', '=', 'N')
  
  if(grepl(paste0(toMatch, collapse = "|"), unit.x) || grepl(paste0(toMatch, collapse = "|"), unit.y)){
    return(c(sub.s, sub.n))
  }else{
    section = codon[[which(sapply(codon, function(X){unit.x %in% X}))]]
    mut.num = mapply(function(X, Y) sum(X != Y), strsplit(unit.x, ""), strsplit(unit.y, ""))
    mut.pos = mapply(function(X, Y) which(X != Y), strsplit(unit.x, ""), strsplit(unit.y, ""))
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
      return(c(sub.s, sub.n))
    }
  }
  return(c(sub.s, sub.n))
}

# Find all possible substitution pathways.        
# unit.x = "TTT"
# unit.y = "TAC"
# unit.y = "GAC"
# mut.pos = 2:3
# mut.pos = 1:3
collect_allSubs = function(unit.x, unit.y, mut.pos){
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

#Sum the number of substitutions of all possible pathways
cal_subsub = function(unit.list, sub.s, sub.n, sec){
  for (m in 1:length(unit.list)) {
    sec.var = sec
    for (n in 2:(length(unit.list[[m]])) ) {
      if(unit.list[[m]][n] %in% sec.var){
        sub.s = sub.s + 1
      }else{
        sub.n = sub.n + 1
      }
      sec.var = codon[[which(sapply(codon, function(X){unit.list[[m]][n] %in% X}))]]
    }
  }
  return(c(sub.s, sub.n))
}


##############################
#get the number of typeN, typeS indels (<=12nts)
get_typeNS_sw = function(dna,M,R){
  dna.str  = str_split(dna,'')
  g        = IRangesList(lapply(dna.str, function(x){IRanges(x=='-')}))
  
  g1    = g[[1]][which(width(g[[1]])<=12)] 
  g2    = g[[2]][which(width(g[[2]])<=12)] 
  wid.m = width(g1)
  wid.r = width(g2)
  pos.m = start(g1)
  pos.r = start(g2)
  l.m   = length(wid.m)
  l.r   = length(wid.r)
  
  eff.1 = 0
  non.1 = 0
  eff.2 = 0
  non.2 = 0
  
  if(l.m+l.r == 0){
    res = c(0,0)
  }else{
    if(l.m>0){#gaps in mouse
      eff_M = eff_phase(R,pos.m,wid.m,l.m)    
      eff.1 = eff.1 + eff_M[1]
      non.1 = non.1 + eff_M[2]
    }
    if(l.r>0){#gaps in rat
      eff_R = eff_phase(M,pos.r,wid.r,l.r)
      eff.2 = eff.2 + eff_R[1] 
      non.2 = non.2 + eff_R[2]
    }
    typeN  = eff.1 + eff.2
    typeS  = non.1 + non.2
  }
  return(c(typeN,typeS))
}

#Number of effective/non-effecitive indels
eff_phase = function(seq1,pos,wid,len){#seq1:ref
  seq1.char = str_split(seq1,"")[[1]]
  eff.sub = 0
  non.sub = 0
  for(j in 1: len){
    if(pos[j]%%3 == 1){#phase-0
      non.sub = non.sub+1
    }else if(pos[j]%%3 == 0){#phase-2
      pos.ori = pos[j]-2 
      unit.1  = substr(seq1, pos.ori, pos.ori+2)
      unit.2  = substr(seq1, pos.ori+wid[j], pos.ori+wid[j]+2)
      sub     = paste0(c(seq1.char[pos[j]-2], seq1.char[pos[j]-1], seq1.char[pos[j]+wid[j]]), collapse="")
      sec     = codon[[which(sapply(codon, function(X){sub %in% X}))]]
      if(unit.1%in%sec || unit.2%in%sec){
        non.sub = non.sub+1
      }else{
        eff.sub = eff.sub+1
      }
    }else{#phase-1
      pos.ori = pos-1 
      unit.1  = substr(seq1, pos.ori, pos.ori+2)
      unit.2  = substr(seq1, pos.ori + wid[j], pos.ori+wid[j]+2)
      sub     = paste0(c(seq1.char[pos[j]-1], seq1.char[pos[j]+wid[j]], seq1.char[pos[j]+wid[j]+1]), collapse="")
      sec     = codon[[which(sapply(codon, function(X){sub %in% X}))]]
      if(unit.1%in%sec || unit.2%in%sec){
        non.sub = non.sub+1
      }else{
        eff.sub = eff.sub+1
      }
    }
  }  
  res = c(eff.sub,non.sub) 
  return(res)
}




######################################################################
# inD = "Data/sw_cds/01_FcaCaf_aligned_cds"
# ouF = "Results/ZD_sum/01_FcaCaf_aligned_cds.sw.tsv"

main = function(inD,ouF){
  #source the script
  source("../Script/sources/codon_call.R")
  
  #construct codons and its degeneracy
  co.res = codon_call()
  syn   <<- co.res[[3]]
  
  codon <<- list (c("TTT","TTC"),
                  c("TTA","TTG","CTT","CTC","CTA","CTG"),
                  c("ATT","ATC","ATA"),
                  c("ATG"),
                  c("GTT","GTC","GTA","GTG"),
                  c("TCT","TCC","TCA","TCG","AGT","AGC"),
                  c("CCT","CCC","CCA","CCG"),
                  c("ACT","ACC","ACA","ACG"),
                  c("GCT","GCC","GCA","GCG"),
                  c("TAT","TAC"),
                  c("CAT","CAC"),
                  c("CAA","CAG"),
                  c("AAT","AAC"),
                  c("AAA","AAG"),
                  c("GAT","GAC"),
                  c("GAA","GAG"),
                  c("TGT","TGC"),
                  c("TGG"),
                  c("CGT","CGC","CGA","CGG","AGA","AGG"),
                  c("GGT","GGC","GGA","GGG"),
                  c("TAA","TGA","TAG")
  )
  
  
  Files = list.files(inD,full.names=T)
  n     = length(Files)
  
  datM           = matrix(0,n,8)
  colnames(datM) = c('Zn','Zs','Nd','N','Sd','S','lenA','lenB') 
  g.type         = c('-', '\\+', '=')
  g.pat          = paste0(g.type, collapse="|")
  
  for (i in 1:length(Files)) {
    dna  = readBStringSet(Files[i])
    len  = width(dna)[1]
    M    = toupper(toString(dna[1])) 
    R    = toupper(toString(dna[2]))
    
    
    lenA   = nchar(str_remove_all(dna[[1]],g.pat)) 
    lenB   = nchar(str_remove_all(dna[[2]],g.pat))
    
    #count typeN,typeS
    tyS = tyN = 0
    NS.value  = get_typeNS_sw(dna,M,R)
    tyN       = tyN + NS.value[1]
    tyS       = tyS + NS.value[2]
    
    
    #Cal. the dN/dS of each pair 
    n.m   = 0 #number of non-synonymous sites 
    s.m   = 0 #number of synonymous sites
    n.r   = 0
    s.r   = 0
    Nd    = 0 #number of non-synonymous mutations
    Sd    = 0 #number of synonymouse mutations
   
    j = 1   
    while (j < len) {
      Subs    = cal_subs(M,R,j)
      Sd   = Sd + Subs[1]
      Nd   = Nd + Subs[2]
      
      nSite_M  = cal_sites(M,j)
      n.m      = n.m + nSite_M[1]
      s.m      = s.m + nSite_M[2]
      
      nSite_R  = cal_sites(R,j) 
      n.r      = n.r + nSite_R[1]
      s.r      = s.r + nSite_R[2]
      
      j = j + 3  
    }
    
    N = n.m + n.r
    S = s.m + s.r
    
    datM[i,] = c(tyN,tyS,Nd,N,Sd,S,lenA,lenB)
    print(i)
  }
  
  #Output int tsv file
  write.table(datM, file=ouF,sep='\t', append=FALSE, quote=FALSE, row.names=FALSE, col.names=TRUE)
  
}



args = commandArgs(trailingOnly=TRUE)
main(args[1], args[2])
