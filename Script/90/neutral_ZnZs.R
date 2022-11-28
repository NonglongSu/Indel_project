#neutral Zn/Zn+Zs
#Simulate the Prob of typeN and typeS in phase1 and phase 2 via 64*4 codon space.
#Assume symmetry between insertion and deletion, hence I only simulate insertion here.

suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(seqinr))

#setwd("~/Dropbox (ASU)/Indel_project/test_90_species")

#count the nucleotide freq
count_nuc_freq = function(input){
  nuc.count = 0
  for (i in input){
    dna       = readDNAStringSet(i)
    nuc.count = nuc.count + oligonucleotideFrequency(dna,width=1)
  }
  nuc.freq  = colSums(nuc.count)/sum(nuc.count)
  if(sum(nuc.freq)==1){
    return(nuc.freq)
  }else{
    print("The nucleotide frequency sums up not 1!")
  }
}

#count the codon freq
count_codon_freq = function(input){
  cod.count = 0
  for (i in input){
    dna       = readDNAStringSet(i)
    cod.count = cod.count + trinucleotideFrequency(dna)
  }
  cod.freq  = colSums(cod.count)/sum(cod.count)
  if(sum(cod.freq)==1){
    return(cod.freq)
  }else{
    print("The codon frequency sums up not 1!")
  }
}

#simulated 3-mer del in a 6-mer seq
#ty-1: Zs; ty-0:Zn
sim_del = function(s1,s2,pos){
  
  c1  = str_split(s1, "")[[1]]
  c2  = str_split(s2, "")[[1]]
  
  ty = 0
  if(pos == 1){
    sub = paste0(c1[1],c2[2],c2[3],collapse="")
    sec = syn[[sub]]
    if(s1 %in% sec || s2 %in% sec){ 
      ty = 1  
    }
  }else{
    sub = paste0(c1[1],c1[2],c2[3],collapse="")
    sec = syn[[sub]]
    if(s1 %in% sec || s2 %in% sec){ 
      ty = 1 
    }
  }
  return(ty)
}


#########################################################MAIN
#inD="Raw_data/cds"
#ouF="Results/ZnZs/neutral_ZnZs.txt"
main = function(inD,ouF){
  
  #construct codons and its degeneracy
  codons = cbind(rep(DNA_BASES, each=16),
                 rep(DNA_BASES, times=4, each=4),
                 rep(DNA_BASES, 16))
  codonstrs  = apply(codons, 1, stringr::str_c, collapse="")
  syn        = syncodons(codonstrs)
  names(syn) = toupper(names(syn))
  syn        <<- lapply(syn, toupper)
  
  n = length(codonstrs)
 
  Dirs  = list.files(inD,full.names=TRUE)
  N     = length(Dirs)
  nucf  = matrix(0,N,4)
  codf  = matrix(0,N,64)
  for (i in 1:length(Dirs)) {
    Files = list.files(Dirs[i],full.names=TRUE)
    nucf[i,] = count_nuc_freq(Files)
    codf[i,] = count_codon_freq(Files)
    print(i)
  }
  print(rowSums(nucf))
  print(rowSums(codf))
  
  
  
  set.seed(8088)
  Zrat = matrix(0,n,n)
  for (i in 1:n) {
    s1 = codonstrs[i]
    for(j in 1:n) {
      s2      = codonstrs[j]
      sam.pos = sample(0:2,10000,replace=T)
      ph0     = length(which(sam.pos==0))
      ph1     = length(which(sam.pos==1))
      ph2     = length(which(sam.pos==2))
      
      typ1    = sim_del(s1, s2, 1)
      typ2    = sim_del(s1, s2, 2)
      
      Zn        = ph1*as.integer(typ1==0) + ph2*as.integer(typ2==0)
      Zs        = ph0 + ph1*as.integer(typ1==1) + ph2*as.integer(typ2==1)
      Zrat[i,j] = Zn/(Zn+Zs)
    }
  }
   
  
  neu  = c()
  neu1 = c()
  for(i in 1:N){
    f0  = nucf[i,]
    cf0 = sapply(seq(n), function(x){prod(f0[match(codons[x,],DNA_BASES)])})
    
    cf1 = codf[i,]
    
    dmat = matrix(0,n,n)
    dmat1= matrix(0,n,n)
    for(j in 1:n){
      for(k in 1:n){
        dmat[j,k]  = cf0[j]*cf0[k]*Zrat[j,k]
        dmat1[j,k] = cf1[j]*cf1[k]*Zrat[j,k]
      }
    }
    neu[i]  = sum(dmat)
    neu1[i] = sum(dmat1)
  }
  
  
  specname = basename(Dirs)
  df       = data.frame('species'=specname,'neutral_nuc_ratio'=neu, 'neutral_cod_ratio'=neu1)
  write.table(df,ouF,col.names=T,row.names=F,sep='\t',quote=F)
}

#####################################
args=commandArgs(trailingOnly=T)
main(args[1],args[2])

