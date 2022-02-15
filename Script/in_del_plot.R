library(Biostrings)
library(stringr)

phasing = function(file){
  
  dna = readBStringSet(file, format = "fasta")
  l   = length(dna)
  
  dna.1 = str_split(as.character(dna), '')
  g1 = lapply(dna.1, function(x) { IRanges(x %in% c("-", "+"))})
  g2 = lapply(dna.1, function(x) { IRanges(x == "+")})
  g1 = IRangesList(g1)
  g2 = IRangesList(g2)
  
  mr     = g1[(l-1):l]
  wid.mr = unlist(width(mr))
  pos.mr = unlist(start(mr))
  
  # record "+++"
  flag  = rep(0, length(wid.mr))
  flag1 = unlist(start(g2[(l-1):l]))
  flag[which(pos.mr %in% flag1)] = rep(1, length(flag1))
  
  # record file
  files = rep(basename(file), length(pos.mr))
  
  # record which species
  spec = c(rep('m', length(mr[[1]])), rep('r', length(mr[[2]])))
  
  df.mr = data.frame("pos" = pos.mr, "wid" = wid.mr, "flag" = flag, "file" = files, "spec" = spec)
  return(df.mr)
}  


indel_watch = function(f3, psj) {
  dna     =  readBStringSet(f3, format = "fasta")
  dna.ref = dna[[1]]
  pos = psj$pos
  wid = psj$wid
  tag = psj$spec
  
  pat   = substr(dna.ref, pos, pos + wid - 1)
  ins.m = ins.r = del.m = del.r = 0
  
 if(tag == 'm'){
   if(grepl('-', pat) == TRUE){
     ins.r = 1
   }else{
     del.m = 1
   }
 }else{
   if(grepl('-', pat) == TRUE){
     ins.m = 1
   }else{
     del.r = 1
   }
 }
 
  ind.m = ins.m + del.m
  ind.r = ins.r + del.r
  res   = c(ind.m, ind.r)
  return(res)
}


Record = function(PST, gapName){
  pos.3  = PST[PST$wid == 3, ]$pos
  pos.6  = PST[PST$wid == 6, ]$pos
  pos.9  = PST[PST$wid == 9, ]$pos
  pos.12 = PST[PST$wid == 12,]$pos
  
  pos.lst  = list(pos.3, pos.6, pos.9, pos.12)
  DF.phase = data.frame() 
  for (i in 1:length(pos.lst)) {
      rem      = pos.lst[[i]] %% 3
      phase.0  = length(which(rem == 1))
      phase.1  = length(which(rem == 2)) 
      phase.2  = length(which(rem == 0)) 
      df.phase = data.frame(Phase_0 = phase.0, Phase_1 = phase.1, Phase_2 = phase.2)
      DF.phase = rbind(DF.phase, df.phase)
  }
  DF.phase = as.table(noquote(t(DF.phase)))
  colnames(DF.phase) = gapName
  return(DF.phase)
}  

# setwd("~/Dropbox (ASU)/Indel_project/Script")

# file.1 = "../test_human_mouse_rat/Data_6/Mafft/mapped_cds/ENSG00000000460.fa"
# file.2 = "../test_human_mouse_rat/Data_6/Mafft/updated_cds/ENSG00000000460.fa"

# dir1  = "../test_human_mouse_rat/Data_6/Mafft/updated_cds"
# dir2  = "../test_human_mouse_rat/Data_6/Mafft/mapped_cds"
# ouFig = "../test_human_mouse_rat/Data_6/Figure/Phase.mafft.sep.pdf"
# Name  = c('mouse', 'rat')
# num   = '6'

main = function(dir1, dir2, ouFig, Name, num){

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
  
  # Filter out the "+++"
  PST.2f = PST.2[PST.2$flag == 0, ]
  PST.1f = PST.1[which(PST.2$flag == 0), ]  

  # Catch the seq state of ref to seperate indels
  Ind.m = c()
  Ind.r = c()
  for (j in 1:nrow(PST.1f)) {
    file3   = paste0(dir1, "/", PST.1f$file[[j]])
    ind.sep = indel_watch(file3, PST.1f[j, ])
    Ind.m   = c(Ind.m, ind.sep[1])
    Ind.r   = c(Ind.r, ind.sep[2])
  }
  PST.tag = cbind(PST.2f, Ind.m, Ind.r)
  
  
  ### PART II
  PST.m = PST.tag[PST.tag$Ind.m == 1, ]
  PST.r = PST.tag[PST.tag$Ind.r == 1, ]
  
  gap.len  = c('3', '6', '9', '12')
  Phase.m  = Record(PST.m, gap.len)
  Phase.r  = Record(PST.r, gap.len)
  
  # Plot
  pdf(ouFig)
  par(mfrow = c(1, 2), oma = c(2, 0, 2, 0))
  
  safe_colorblind_palette = c("#88CCEE", "#CC6677", "#DDCC77")
  barplot(Phase.m, xlab = Name[1], ylab = "Count", col = safe_colorblind_palette ) 
  barplot(Phase.r, xlab = Name[2], ylab = "Count", col = safe_colorblind_palette ) 
  
  Window = as.numeric(num)
  mtext("The proportion of phased-indels across length of gaps", outer = TRUE, cex = 1.5)
  mtext(paste0("Gap length (window size = ", Window, ")"), outer = TRUE, cex = 1.5, side = 1)
  legend(2.5, 800, bg ="transparent", bty = "n", legend = rownames(Phase.m), fill = safe_colorblind_palette, cex = 0.75)
  dev.off()
  
  
}
  
args = commandArgs(trailingOnly = TRUE)
main(args[1], args[2], args[3], eval(parse(text = args[4])), args[5])
  

