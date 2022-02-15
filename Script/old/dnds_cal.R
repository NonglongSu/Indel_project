library(Biostrings)
library(stringr)
library(seqinr)


#setwd("~/Dropbox (ASU)/Indel_project/Script")

#Number of non-syn/syn sites
cal_sites = function(x,y,j){
  base  = c('A',"T","C","G")
  scoreboard = 3
  score = -1/3
  
  unit.x = substr(x,j,j+2)
  unit.y = substr(y,j,j+2)
  
  if( grepl('-|N',unit.x) || grepl('-|N',unit.y) ){
    return(0)
  }else{
    section = codon[[which(sapply(codon, function(X){unit.x %in% X}))]]
    for (k in 1:3) {
      Base = base[substr(unit.x,k,k) != base]
      unit.sub = unit.x
      for (w in 1:length(Base)) {
        substr(unit.sub,k,k) = Base[w]
        if(unit.sub %in% section){
          scoreboard = scoreboard + score
        }
      }
    }
  }
  return(scoreboard)  
}

#Number of non-syn/syn substitutions
cal_subs = function(x,y,j){
  unit.x = substr(x,j,j+2)
  unit.y = substr(y,j,j+2)
  
  sub.s = 0
  sub.n = 0
  
  if( grepl('-|N',unit.x) || grepl('-|N',unit.y) ){
    return(c(sub.s,sub.n))
  }else{
    section = codon[[which(sapply(codon, function(X){unit.x %in% X}))]]
    mut.num = mapply(function(X,Y) sum(X!=Y),strsplit(unit.x,""),strsplit(unit.y,""))
    mut.pos = mapply(function(X,Y) which(X!=Y),strsplit(unit.x,""),strsplit(unit.y,""))
    if(mut.num == 1){
      if(unit.y %in% section){
        sub.s = sub.s + 1
      }else{
        sub.n = sub.n + 1  
      }
    }else if(mut.num == 2){
      unit.list = collect_allSubs(unit.x, unit.y, mut.pos)
      sub.rate  = cal_subsub(unit.list,sub.s,sub.n,section)
      sub.s = sub.rate[1]/2
      sub.n = sub.rate[2]/2
    }else if(mut.num == 3){
      unit.list = collect_allSubs(unit.x, unit.y, mut.pos)
      sub.rate  = cal_subsub(unit.list,sub.s,sub.n,section)
      sub.s = sub.rate[1]/6
      sub.n = sub.rate[2]/6
    }else{#no substitution
      return(c(sub.s,sub.n))
    }  
  }
  return(c(sub.s,sub.n))
}

#Find all possible substitution pathways.
collect_allSubs = function(unit.x,unit.y,mut.pos){
  unit.set  = c()
  for (i in 1:length(mut.pos)) {
    unit.sub = unit.x 
    substr(unit.sub,mut.pos[i],mut.pos[i]) = substr(unit.y,mut.pos[i],mut.pos[i])
    unit.set = c(unit.set,unit.sub)
    if (length(mut.pos) == 3) {
      for (j in mut.pos[mut.pos!=mut.pos[i]] ) {
        unit.sb = unit.sub 
        substr(unit.sb,j,j) = substr(unit.y,j,j)
        unit.set = c(unit.set,unit.sb)
      }
    }
  }
  unit.list = list()
  if(length(unit.set)==2){
    for (k in 1:length(unit.set)) {
          unit.list = c(unit.list,list(c(unit.x,unit.set[k],unit.y)))
        }
  }else{
    for (k in 1:length(unit.set)) {
      if (k %% 3 == 1) {
        for (w in (k+1):(k+2)) {
          unit.list = c(unit.list,list(c(unit.x,unit.set[k],unit.set[w],unit.y)))
        }
      }
    }
  }
  return(unit.list)
}

#Sum the number of substitutions of all possible pathways 
cal_subsub = function(unit.list,sub.s,sub.n,sec){
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
  return(c(sub.s,sub.n))
}



# inDir  = "../test_human_mouse_rat/Data_6.2/Mafft/mapped_anc"
# ouFile = "../test_human_mouse_rat/Data_6.2/Results/dNdS.txt"

main = function(inDir, ouFile){
  
  # Create a codon table by hand
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
  
  Files = list.files(inDir, full.names = TRUE)
  
  DF = c()
  for (i in 1:length(Files)) {
    dna   = readDNAStringSet(Files[i] ,format = "fasta")
    len   = width(dna)[1]
    
    A = toString(dna[[1]])
    M = toString(dna[[2]])
    R = toString(dna[[3]])
    
    #A1 = str_split(A,"")
    #Cal. the dN/dS with mouse
    n.m  = 0 # number of non-synonymous sites (relative to species)
    Nd.m = 0 # number of non-synonymous mutations
    Sd.m = 0 # number of synonymouse mutations
    j = 1
    while (j < len) {
      nSite  = cal_sites(A,M,j) 
      Nd   = cal_subs(A,M,j)
      n.m  = n.m + nSite
      Sd.m = Sd.m + Nd[1]
      Nd.m = Nd.m + Nd[2]
      j = j+3
    } 
    s.m = len - n.m  # S = (3r-N)
    P.n_m = Nd.m/n.m
    P.s_m = Sd.m/s.m
    
    #Cal. the dN/dS with rat
    n.r = 0
    Nd.r = 0
    Sd.r = 0
    j = 1
    while (j < len) {
      nSite  = cal_sites(A,R,j)
      Nd   = cal_subs(A,R,j)
      n.r  = n.r + nSite
      Sd.r = Sd.r + Nd[1]
      Nd.r = Nd.r + Nd[2]
      j = j+3
    } 
    s.r = len - n.r  
    P.n_r = Nd.r/n.r
    P.s_r = Sd.r/s.r
    
  #Jukes Cantor formula (1969)
    dNdS_M = log(1-4*P.n_m/3,exp(1)) / log(1-4*P.s_m/3,exp(1))
    dNdS_R = log(1-4*P.n_r/3,exp(1)) / log(1-4*P.s_r/3,exp(1))
    
  #Generate a database  
    df = data.frame("File.name"=basename(Files[i]), "M_dNdS"=dNdS_M, "R_dNdS"=dNdS_R)
    DF = rbind(DF,df)
  }
  
  
  write.table(DF, file = ouFile, 
              sep = "\t", append = FALSE, quote = FALSE, row.names = FALSE, col.names = TRUE )
}

args = commandArgs(trailingOnly = TRUE)
main(args[1],args[2])


