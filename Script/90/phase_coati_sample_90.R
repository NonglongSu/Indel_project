#Calculate the proportion of phase-0/1/2 via coati-sampling method. 
suppressWarnings(suppressMessages(library(tidyverse)))
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(seqinr))
suppressPackageStartupMessages(library(jsonlite))
suppressPackageStartupMessages(library(ggplot2))  
suppressPackageStartupMessages(library(ggrepel))



#setwd("~/Dropbox (ASU)/Indel_project/Script/90")

#only grap the positions of gaps of 3,6,9,12.
phasing_sample = function(fi){
  dna.json = fromJSON(fi)
  
  minw     = min(dna.json$log_weight)
  norm.wei = exp(dna.json$log_weight-minw)/sum(exp(dna.json$log_weight-minw))
  #print(sum(norm.wei))
  
  phaseV = matrix(0,nrow(dna.json),3)
  for(i in 1:nrow(dna.json)){
    dna   = DNAStringSet(c(dna.json$aln[i,1],dna.json$aln[i,2]))
    dna.1 = str_split(dna,'')
    g     = lapply(dna.1, function(x) {IRanges(x=='-')})
    g     = IRangesList(g)
    ug    = unlist(g)
    ug    = ug[which(width(ug)<=12)]
    
    if(length(ug)>0){
      pos       = start(ug)
      phaseV[i,]= Record(pos)
    }
  }
  norm.phaseV = norm.wei*phaseV
  
  #print(colSums(norm.phaseV))
  if(all(is.nan(colSums(norm.phaseV)))){stop("NAN error")}
  
  return(colSums(norm.phaseV))
}

Record = function(PST){
  rem      = PST %% 3
  phase.0  = length(which(rem == 1))
  phase.1  = length(which(rem == 2)) 
  phase.2  = length(which(rem == 0)) 
  phaseV   = c(phase.0,phase.1,phase.2)
  return(phaseV)
}  


##############################################################

# dir   = "../../test_90_species/Raw_data/align_sampling/01_FcaCaf_aligned_cds"
# ouF   = "../../test_90_species/Results/phase_sample/01_FcaCaf_aligned_cds.txt"

main = function(dir,ouF){
  Files = list.files(dir,full.names=TRUE)
  pst   = c(0,0,0)
  for(i in 1:length(Files)){
    phase_score = phasing_sample(Files[i])
    pst         = pst + phase_score   
    print(i)
  }
  
  pst = pst/sum(pst)
  
  #Output a line
  write(pst,ouF,sep='\t')
}

#######################################
args = commandArgs(trailingOnly=TRUE)
main(args[1], args[2])



