#Calculate the proportion of phase-0/1/2 via coati-max method. 
suppressWarnings(suppressMessages(library(tidyverse)))
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(seqinr))
suppressPackageStartupMessages(library(ggplot2))  
suppressPackageStartupMessages(library(ggrepel))


#setwd("~/Dropbox (ASU)/Indel_project/test_90_species")

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

# dir   = "Raw_data/align_max"
# ouF   = "Results/znzs/phase_coati.txt"

main = function(dir, ouF){
  
  Files  = list.files(dir,full.names=TRUE)
  p.name = c("phase-0","phase-1","phase-2")
  
  PST = list()
  for(i in 1:length(Files)){
    files = list.files(Files[i],full.names=TRUE)
    pst   = c()
    for(j in files){
      phase_score = phasing_3(j)
      pst         = c(pst,phase_score)
    }
    PST[[i]] = pst
    print(i)
  }
  
  
  phaseM = matrix(0,90,3)
  for (k in 1:length(PST)) {
    phaseM[k,] = Record(PST[[k]])
  }
  PhaseM = phaseM/(rowSums(phaseM))
  
  #Output a table
  Pdf.out   = PhaseM
  spec.name = unlist(lapply(Files,function(x){basename(x)}))
  Pdf.out   = data.frame(Pdf.out)
  colnames(Pdf.out) = p.name
  Pdf.out   = cbind(spec.name,Pdf.out)
  write.table(Pdf.out,ouF,row.names=F,sep='\t',quote=F)
  #capture.output(print(Pdf.out, print.gap=3),file=ouF)
  
}

#######################################
args = commandArgs(trailingOnly=TRUE)
main(args[1], args[2])



# test_3 = function(f){
#   dna   = readBStringSet(f)
#   dna.1 = str_split(dna,'')
#   g     = lapply(dna.1, function(x) {IRanges(x == '-')})
#   g     = IRangesList(g)
#   ug    = unlist(g)
#   ug    = ug[which(width(ug)%%3 != 0)]
#   
#   if(length(ug)>0){
#     return(f)
#   }
# }
# 
# Test = list()
# for(i in 1:length(Files)){
#   files = list.files(Files[i],full.names=TRUE)
#   tes   = c()
#   for(j in files){
#     outli = test_3(j)
#     tes         = c(tes,outli)
#     if(length(tes)!=0){
#       print(j)
#       break}
#   }
#   Test[[i]] = tes
#   print(i)
# }



