#Data_6/mafft_edge data
#cateogrize the bias data to ECM (https://academic.oup.com/mbe/article/24/7/1464/986344)

suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(stringr))

#setwd("~/Dropbox (ASU)/Indel_project/test_human_mouse_rat")

#1-same cat, 0-otherwise
aa_cat = function(x,mod){
  sv = rep(0,nrow(x))
  for(i in 1:nrow(x)){
    aaA = str_split(x$amino_acid.1[i],'')
    aaB = str_split(x$amino_acid.2[i],'')
    aa  = sapply(list(aaA,aaB), "[[", 1)
    
    for(j in 1:nrow(aa)){
      sec = mod[[aa[j,1]]]
      if(aa[j,2] %in% sec){
        sv[i] = sv[i] + 1
      }
    }
  }
  sv
}


######################################################

inF1 = "Data_6/Figure/mafft_edge/align.ori.l.txt"
inF2 = "Data_6/Figure/mafft_edge/align.better.l.txt"
inF3 = "Data_6/Figure/mafft_edge/align.ori.r.txt"
inF4 = "Data_6/Figure/mafft_edge/align.better.r.txt"

#make a list according to the ECM paper
ecm = list(  W='W',Y=c('Y','F'),F=c('Y','F'),
             L=c('L','M','I','V'),M=c('L','M','I','V'),I=c('L','M','I','V'),V=c('L','M','I','V'),
             C='C',
             T=c('T','S','A','E','D','N','Q','K','R','H'),S=c('T','S','A','E','D','N','Q','K','R','H'),A=c('T','S','A','E','D','N','Q','K','R','H'),
             E=c('T','S','A','E','D','N','Q','K','R','H'),D=c('T','S','A','E','D','N','Q','K','R','H'),N=c('T','S','A','E','D','N','Q','K','R','H'),
             Q=c('T','S','A','E','D','N','Q','K','R','H'),K=c('T','S','A','E','D','N','Q','K','R','H'),R=c('T','S','A','E','D','N','Q','K','R','H'),
             H=c('T','S','A','E','D','N','Q','K','R','H'),
             G='G',
             P='P')


mec = list(  W='W',F=c('L','F'),L=c('L','F'),
             M=c('M','I','V','E','D','Q','K'),I=c('M','I','V','E','D','Q','K'),V=c('M','I','V','E','D','Q','K'),
             E=c('M','I','V','E','D','Q','K'),D=c('M','I','V','E','D','Q','K'),Q=c('M','I','V','E','D','Q','K'),
             K=c('M','I','V','E','D','Q','K'),
             C=c('C','T','S','A','N','R','G'),T=c('C','T','S','A','N','R','G'),S=c('C','T','S','A','N','R','G'),
             A=c('C','T','S','A','N','R','G'),N=c('C','T','S','A','N','R','G'),R=c('C','T','S','A','N','R','G'),
             G=c('C','T','S','A','N','R','G'),
             R=c('R','H','Y'),H=c('R','H','Y'),Y=c('R','H','Y'),
             P=c('P','S'),S=c('P','S'))

dat1 = read.table(inF1,header=T)
dat2 = read.table(inF2,header=T)
dat3 = read.table(inF3,header=T)
dat4 = read.table(inF4,header=T)

Dat1 = rbind(dat1,dat3)
Dat2 = rbind(dat2,dat4)

score1 = aa_cat(Dat1,ecm)
score2 = aa_cat(Dat2,ecm)

slen = length(which(score1 != score2))

cat(sprintf("Percentage of original alignment in ECM: %.3f\n",length(which(score1>score2))/slen))  #0.845






