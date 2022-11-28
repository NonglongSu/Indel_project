#Calculate the proportion of phase-0/1/2 via coati-alignpair method. 

suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(seqinr))
suppressPackageStartupMessages(library(ggplot2))  
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(tidyverse))


#setwd("~/Dropbox (ASU)/Indel_project/Script")

#only grap the positions of gaps of 3,6,9,12.
phasing_3 = function(f){
  dna   = readBStringSet(f)
  dna.1 = str_split(dna,'')
  g     = lapply(dna.1, function(x) {IRanges(x == '-')})
  g     = IRangesList(g)
  ug    = unlist(g)
  ug    = ug[which(width(ug)<=12)]
  
  pos = c()
  if(length(ug)>0){
    pos = start(ug)
  }
  return(pos)
}

Record = function(PST){
  rem      = PST %% 3
  phase.0  = length(which(rem == 1))
  phase.1  = length(which(rem == 2)) 
  phase.2  = length(which(rem == 0)) 
  df.phase = c(phase.0,phase.1,phase.2)
  return(df.phase)
}  


##############################################################
# dir   = "../test_human_mouse_rat/Raw_data/coati_align"
# ouF   = "../test_human_mouse_rat/Results/phase_coati.txt"

main = function(dir, ouF){
  
  Files  = list.files(dir,full.names=TRUE)
  p.name = c("phase-0","phase-1","phase-2")
  
  pst = c()
  for(i in 1:length(Files)){
    phase_score = phasing_3(Files[i])
    pst         = c(pst,phase_score)
    print(i)
  }
  
  phaseM = Record(pst)
  PhaseM = phaseM/(sum(phaseM))
  
  #output
  names(PhaseM) = p.name
  capture.output(print(PhaseM, print.gap=3),file=ouF)
  
}

#######################################
args = commandArgs(trailingOnly=TRUE)
main(args[1], args[2])