#Generate a csv/tsv file of (<=12nts)
# - ZN
# - Zs
# - ZN
# - ZS
# - Nd
# - N
# - Sd
# - S
# - lenA
# - lenB

suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(seqinr))
suppressPackageStartupMessages(library(dplyr))

#setwd("~/Dropbox (ASU)/Indel_project/test_human_mouse_rat")
########################################################


#Determine number of non-syn/syn sites
cal_sites = function(x,tag){
  
  unit.x  = substr(x, tag, tag+2)
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



#Number of non-syn/syn indel sites (3mer-gap)
#pick one as ancestor
cal_zed_sites = function(anc.seq,wid,len){
  D_mer    = len-wid+1
  anc.char = str_split(anc.seq,"")[[1]] 
  
  s.d = 0
  n.d = 0
  i   = 1
  while (i <= D_mer) {
    if(i %% 3 == 1){#Phase-0
      unit.1 = substr(anc.seq, i, i+2)
      unit.2 = substr(anc.seq, i+wid, i+wid+2)
      if(grepl('N',unit.1) || grepl('N',unit.2)){
        i = i + 3
      }else{
        s.d = s.d + 1
        i   = i + 1
      }
    }else{
      if(i %% 3 == 0){#phase-2
        sub = paste0(c(anc.char[i-2], anc.char[i-1], anc.char[i+wid]), collapse="")
        sec = codon[[which(sapply(codon, function(X){sub%in%X}))]]
        if(unit.1%in%sec || unit.2%in%sec){
          s.d = s.d + 1
        }else{
          n.d = n.d + 1  
        }
      }else{#phase-1
        sub = paste0(c(anc.char[i-1], anc.char[i+wid], anc.char[i+wid+1]), collapse="")
        sec = codon[[which(sapply(codon, function(X){sub%in%X}))]]
        if(unit.1%in%sec || unit.2%in%sec){
          s.d = s.d  + 1
        }else{
          n.d = n.d + 1
        }
      }
      i = i+1
    }
  }
  res = c(s.d,n.d) 
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

#obtain the znzs pars (<=12nts)
get_znzs_par = function(dna){
  Sm = Nm = Sr = Nr = 0
  
  M        = toString(dna[1]) 
  R        = toString(dna[2]) 
  toRemove = paste0(c('-','\\+','='), collapse="|")
  M1       = str_remove_all(M,toRemove)
  R1       = str_remove_all(R,toRemove)
  len.M1   = nchar(M1)
  len.R1   = nchar(R1)
  
  dna.str  = str_split(dna,'')
  seqA     = dna.str[[1]]
  seqB     = dna.str[[2]]
  g        = IRangesList(lapply(dna.str, function(x){IRanges(x == '-')}))
  
  g1 = g[[1]][which(width(g[[1]])<=12)] 
  g2 = g[[2]][which(width(g[[2]])<=12)] 
  
  nSite.M  = cal_zed_sites(M1,3,len.M1)    
  nSite.R  = cal_zed_sites(R1,3,len.R1)
  s.Site   = (nSite.M[1] + nSite.R[1])/2
  n.Site   = (nSite.M[2] + nSite.R[2])/2
  
  if(length(g1)+length(g2) == 0){
    res = c(s.Site,n.Site,0,0)
  }else{
    if(length(g1)>0){#gaps in mouse
      subs = cal_zed_subs(seqA,seqB,g1) 
      Nm   = Nm + subs[1]
      Sm   = Sm + subs[2]
      
    }
    if(length(g2)>0){#gaps in rat
      subs = cal_zed_subs(seqB, seqA, g2)
      Nr   = Nr + subs[1]
      Sr   = Sr + subs[2]
      
    }
    s.indel = Sm + Sr 
    n.indel = Nm + Nr
    res     = c(s.Site,n.Site,s.indel,n.indel)
  }
  
  return(res)
}



######################################################################
inD   = "Raw_data/coati_align"
ouF   = "Results/dnds_znzs_sum.tsv"

main = function(inD, ouF){
  
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
  
  #source the script
  source("../Script/sources/codon_call.R")
  co.res = codon_call()
  syn   <<- co.res[[3]]
  
  #read file
  Files = list.files(inD,full.names=T)
  n     = length(Files)
  datM           = matrix(0,n,10)
  colnames(datM) = c('Zn','Zs','ZN','ZS','Nd','Sd','N','S','lenA','lenB') 
  
  
  for (i in 1:length(Files)) {
    dna = readBStringSet(Files[i])
    len = width(dna)[1]
    
    M   = toupper(toString(dna[1]))
    R   = toupper(toString(dna[2]))
    
    ##seq len
    lenA = nchar(str_remove_all(dna[[1]],'-'))
    lenB = nchar(str_remove_all(dna[[2]],'-'))
    
    
    
    ##Zn,Zs,ZN,ZS
    ZS = ZN = Zs = Zn = 0
    zed.value = get_znzs_par(dna)
    ZS        = ZS + zed.value[1]
    ZN        = ZN + zed.value[2]
    Zs        = Zs + zed.value[3]
    Zn        = Zn + zed.value[4]
    
    
    #N,S,Nd,Sd
    n.m   = 0 #number of non-synonymous sites (relative to species)
    s.m   = 0 #number of synonymous sites (relative to species)
    n.r   = 0
    s.r   = 0
    N.sub = 0 #number of non-synonymous mutations
    S.sub = 0 #number of synonymouse mutations
    
    j = 1   
    while (j < len) {
      Subs    = cal_subs(M, R, j)
      S.sub   = S.sub + Subs[1]
      N.sub   = N.sub + Subs[2]
      
      nSite_M  = cal_sites(M, j)
      n.m      = n.m + nSite_M[1]
      s.m      = s.m + nSite_M[2]
      
      nSite_R  = cal_sites(R, j) 
      n.r      = n.r + nSite_R[1]
      s.r      = s.r + nSite_R[2]
      
      j = j + 3  
    }
    N = (n.m + n.r)/2
    S = (s.m + s.r)/2
    
    datM[i,] = c(Zn,Zs,ZN,ZS,N.sub,S.sub,N,S,lenA,lenB)
    
    print(i)
  }
  
  
  #Output a tsv file
  write.table(datM, file=ouF,sep='\t', append=FALSE, quote=FALSE, row.names=FALSE, col.names=TRUE)
}



###########################################
args = commandArgs(trailingOnly=T)
main(args[1], args[2])