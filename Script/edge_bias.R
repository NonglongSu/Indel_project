#Generate the Edge bias source file(txt)

library(Biostrings)
library(BiocGenerics)
library(parallel)
library(S4Vectors)
library(seqinr)
library(stringr)
library(tidyr)

# setwd("~/Dropbox (ASU)/Indel_project/Script")

phasing = function(file){
  dna   = readBStringSet(file, format = "fasta")
  l     = length(dna)
  dna.1 = str_split(as.character(dna), '')
  
  g1 = lapply(dna.1, function(x) { IRanges(x %in% c("-", "+"))})
  g2 = lapply(dna.1, function(x) { IRanges(x == "+")})
  g1 = IRangesList(g1)
  g2 = IRangesList(g2)
  
  mr     = g1[(l-1):l]
  wid.mr = unlist(width(mr))
  pos.mr = unlist(start(mr))
  
  flag  = rep(0, length(wid.mr))
  flag1 = unlist(start(g2[(l-1):l]))
  flag[which(pos.mr %in% flag1)] = rep(1, length(flag1))
  
  files = rep(file, length(pos.mr))
  df.mr = data.frame("pos" = pos.mr, "wid" = wid.mr, "file" = files, "flag" = flag)
  
  return(df.mr)
}

Extract = function(file, idx, wid){
  dna   = readBStringSet(file, format = "fasta")
  l     = length(dna)
  x = toString(dna[[l-1]])
  y = toString(dna[[l]])
  
  start = idx - Window
  stop  = idx + Window + wid - 1
  wid.1 = substr(x, start, stop)
  wid.2 = substr(y, start, stop)
  
  gap   = grep('-', wid.1)
  gap.l = length(gap)
  
  if(gap.l > 0){ #keep wid.1 as reference
    swap(wid.1, wid.2)
  }
  df = data.frame("wid.1" = wid.1, "wid.2" = wid.2)
  return(df)
}


table_magic = function(x, k, y, z){
  pos.edge  = x$pos[which(y == z)]
  wid.edge  = x$wid[which(y == z)]
  
  file.list = as.character(k$file[which(y == z)])
    
  windows = c()
  for(i in 1:length(file.list)){
    val     = Extract(file.list[i], pos.edge[i], wid.edge[i]) 
    windows = rbind(windows, val)
  }
  amino_acid.1 = Biostrings::translate(DNAStringSet(windows$wid.1), 
                                       genetic.code = GENETIC_CODE, if.fuzzy.codon = "solve")
  windows$amino_acid.1 = as.character(amino_acid.1)
  
  wid2 = str_split(windows$wid.2, "")
  amino_acid.2 = list()
  for (i in 1:length(wid2)){
    wid2.char = translate(wid2[[i]], NAstring = "-", ambiguous = FALSE)
    wid2.str  = paste0(wid2.char, collapse = "")
    amino_acid.2[i] = wid2.str
  }
  windows$amino_acid.2 = unlist(amino_acid.2)
  
  return(windows)
}

# dir1   = "../test_human_mouse_rat/Data_6/Mafft/updated_cds"
# dir2   = "../test_human_mouse_rat/Data_6/Mafft/mapped_cds"

main = function(dir1, dir2, ouDir, num){
  Files = list.files(dir2, full.names = FALSE)
  
  PST.1 = c()
  PST.2 = c()
  for(i in 1:length(Files)){
    file1 = paste0(dir1, "/", Files[i])
    file2 = paste0(dir2, "/", Files[i])
    
    pScore.1 = phasing(file1)
    pScore.2 = phasing(file2)
    PST.1    = rbind(PST.1, pScore.1)
    PST.2    = rbind(PST.2, pScore.2)
  }
  
  # Filter the "+++'
  PST.2f = PST.2[PST.2$flag == 0, ]
  PST.1f = PST.1[which(PST.2$flag == 0), ]   
  

  # Edge value
  num    = as.numeric(num) 
  Window <<- num 
  Edge   = c(-Window, Window)                           
  diff.T = PST.2f$pos - PST.1f$pos 
  
  
  # Left bias
  align_ori.l    = table_magic(PST.1f, PST.1f, diff.T, Edge[1]) 
  align_better.l = table_magic(PST.1f, PST.2f, diff.T, Edge[1])  
  
  # Right bias
  align_ori.r    = table_magic(PST.1f, PST.1f, diff.T, Edge[2]) 
  align_better.r = table_magic(PST.1f, PST.2f, diff.T, Edge[2])  
  

  #Write out
  ouFile.1 = "align.ori.l.txt"
  ouFile.2 = "align.better.l.txt"
  ouFile.3 = "align.ori.r.txt"
  ouFile.4 = "align.better.r.txt"
 
  align.edge = list(align_ori.l, align_better.l, align_ori.r, align_better.r)
  ouFile     = c(ouFile.1, ouFile.2, ouFile.3, ouFile.4)
  
  for(i in 1:length(ouFile)){
    write.table(align.edge[i], file = paste0(ouDir, basename(ouFile[i])), 
                sep = "\t", append = FALSE, quote = FALSE, row.names = FALSE, col.names = TRUE)
  }
  
}

#######################################
args = commandArgs(trailingOnly = TRUE)
main(args[1], args[2], args[3], args[4])



#wind.l = unite(windows.l,"window",c(amino_acid.1,wid.1,wid.2,amino_acid.2),sep='\n',remove = TRUE)
#windows.l$wid.12 = paste(windows.l$wid.1,windows.l$wid.2, sep = '\n')
