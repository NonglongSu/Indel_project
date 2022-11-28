#check if a single gene gives many indels (default>=10)
#usage:Rscript --vanilla ../Script/check_more_gap.R Data_6/Mafft/mapped_cds Data_6/one_gene_more_gap.txt 10

suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(stringr))
#setwd("~/Dropbox (ASU)/Indel_project/Script")

# inD = "../test_human_mouse_rat/Data_6/Mafft/mapped_cds"
# ouF = "../test_human_mouse_rat/Data_6/one_gene_more_gap.txt"
# num = '10'
main = function(inD,ouF,num){
  flag = as.numeric(num)
  Files= list.files(inD,full.names=T)
  
  for(i in 1:length(Files)){
    dna = readBStringSet(Files[i])
    l   = length(dna)
    dna = dna[(l-1):l] 
    
    dnas  = str_split(dna,'')
    g     = IRangesList(lapply(dnas, function(x){IRanges(x %in% c("-", "+","="))}))
    num.g = length(width(unlist(g)))
    if(num.g>flag){
      cat(basename(Files[i]),file=ouF,append=T,sep='\n')
      print(i)
    }
  }
  
}


args = commandArgs(trailingOnly=T)
main(args[1],args[2],args[3])