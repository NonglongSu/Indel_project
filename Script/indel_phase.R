# Calculate the proportion of phase 0 / phase 1 / phase 2 indels 
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(stringr))

#setwd("~/Dropbox (ASU)/Indel_project/test_human_mouse_rat/Data_6")

#Generate dataframe: pos | width
phasing = function(file){
  dna = readBStringSet(file, format = "fasta")
  l   = length(dna)
  
  dna.1 = str_split(as.character(dna), '')
  g = lapply(dna.1, function(x) { IRanges(x == '-')})
  g = IRangesList(g)
  
  m = g[[l-1]]
  r = g[[l]]
  
  wid.m = width(m)
  wid.r = width(r)
  
  l.m = length(wid.m)
  l.r = length(wid.r)
  
  pos.m = start(m)
  pos.r = start(r)
  
  
  df.m = data.frame("pos" = pos.m, "wid" = wid.m)
  df.r = data.frame("pos" = pos.r, "wid" = wid.r)
  
  if(l.m > 0 & l.r == 0){
    return(df.m)
  }else if(l.m == 0 & l.r > 0){
    return(df.r)
  }else if(l.m > 0 & l.r > 0){
    df = merge(df.m, df.r, all = TRUE)
    return(df)
  }else{
    return(NULL)
  }
  
}

#Distinguish phases
Record = function(PST){
  if(is.null(PST)){
    return(PST)
  }
  rem   = PST %% 3
  phase.0 = length(which(rem == 1))
  phase.1 = length(which(rem == 2)) 
  phase.2 = length(which(rem == 0)) 
  
  df.phase = data.frame(Phase_0 = phase.0, Phase_1 = phase.1, Phase_2 = phase.2)
  return(df.phase)
}  

#Generate the phase dataframe
phase_gen = function(inD){
  File.total = list.files(inD, full.names=TRUE)
  PST = c()
  for(i in 1:length(File.total)){
    phase_score = phasing(File.total[i])
    PST         = rbind(PST, phase_score)
    print(i)
  }
  
  pos.3  = PST[PST$wid == 3, ]$pos
  pos.6  = PST[PST$wid == 6, ]$pos
  pos.9  = PST[PST$wid == 9, ]$pos
  pos.12 = PST[PST$wid == 12,]$pos
  
  phase.3  = Record(pos.3)
  phase.6  = Record(pos.6)
  phase.9  = Record(pos.9)
  phase.12 = Record(pos.12)
  
  Phase.df  = rbind(phase.3, phase.6, phase.9, phase.12)
  Phase.df
}

#ggplot of phase prop
gg_phase = function(dat,win){
  df      = as.table(noquote(t(dat)))
  colnames(df) = c('3','6','9','12')
  safe_colorblind_palette = c("#0072B2","#009E73","#AA4499")
  
  barplot(df, main=paste0("Window size of ",win), xlab="Gap length", ylab="Count", 
          col=safe_colorblind_palette) 
  legend("topright", bg="transparent", bty="n", rownames(df), fill=safe_colorblind_palette, cex=1)
}

###############################################
#ouFig  = "Figure/Phase.mafft.pdf"
main = function(ouFig){
  
  dir1  = "../Data_3/Mafft/mapped_cds"
  dir2  = "Mafft/mapped_cds"
  dir3  = "../Data_9/Mafft/mapped_cds"
  dir4  = "../Data_12/Mafft/mapped_cds"
  
  pha1 = phase_gen(dir1)
  pha2 = phase_gen(dir2)
  pha3 = phase_gen(dir3)
  pha4 = phase_gen(dir4)
  
  
  #ggplot
  pdf(ouFig,onefile=T)
  par(mfrow=c(2,2))
  gg_phase(pha1,3)
  gg_phase(pha2,6)
  gg_phase(pha3,9)
  gg_phase(pha4,12)
  
  dev.off()   
  
}

#######################################
args = commandArgs(trailingOnly=TRUE)
main(args[1])
