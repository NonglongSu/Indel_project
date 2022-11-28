#Cal. the mean znzs value of each bin groups
#Correlate omega vs znzs

suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(seqinr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(readr))

#setwd("~/Dropbox (ASU)/Indel_project/test_human_mouse_rat")
########################################################

#Number of non-syn/syn indel sites
#dna=DNAStringSet(c("A---TCGATCGT","AAAA---TACAA"))
#site 5 2   subs 1 0
#     6 1        0 1 
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


##################################################################
# inD   = "Raw_data/coati_align"
# file1 = "Results/omega_bin/omega_gap.txt"
# file2 = "Results/omega_bin/omega_bin.txt"
# ouFig = "Figure/omega_bin/omega_znzs.pdf"

main = function(inD, file1, file2, ouFig){
  
  #source the script
  source("../Script/sources/codon_call.R")
  co.res = codon_call()
  syn   <<- co.res[[3]]
  
  
  codon <<-  list(c("TTT","TTC"),
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
  
  #Extract data from files
  Files   = list.files(inD,full.names=T)
  DF.rank = read_delim(file1, "\t", col_names=TRUE)
  geneId  = DF.rank[[1]]
  dnds    = DF.rank[[2]]
  index   = read_delim(file2, "\t", col_names=FALSE)[[1]]
  
  
  #take the average of each bin group of dnds.
  dnds.avg=c()
  for(i in 1:length(index)){
    if(i==1){
      dnds.avg[i] = mean(dnds[i:index[i]])
    }else{
      dnds.avg[i] = mean(dnds[(index[i-1]+1):index[i]])
    }
  }
  
  #cal. the average of each bin group of znzs.
  znzs.avg = c()
  k        = 1
  for (i in index) {
    S.site = N.site = S.indels = N.indels = 0
    while (k<=i) {
      dna       = readBStringSet(Files[geneId[k]])
      zed.value = get_znzs_par(dna)
      S.site    = S.site + zed.value[1]
      N.site    = N.site + zed.value[2]
      S.indels  = S.indels + zed.value[3]
      N.indels  = N.indels + zed.value[4]
      k         = k + 1
    }
    znzs     = (N.indels/N.site)/(S.indels/S.site)
    znzs.avg = c(znzs.avg, znzs)
    print(i)
  }
  
  W.zed = data.frame("znzs"=znzs.avg, "dnds"=dnds.avg)
  
  
  # #Plot
  pdf(ouFig)
  plot(dnds.avg, znzs.avg,
       main = "znzs vs dnds",
       xlab = "dnds", ylab = "znzs", pch=20, col = "#88CCEE")
  dev.off()
  
}



args = commandArgs(trailingOnly = TRUE)
main(args[1], args[2], args[3], args[4])
