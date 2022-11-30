# Calculate the positions  of phase 0 / phase 1 / phase 2 indels 
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggpubr))

#Generate dataframe: pos | width
phasing = function(file){
  dna = readBStringSet(file, format="fasta")
  l   = length(dna)
  
  dna.1 = str_split(as.character(dna), '')
  g = lapply(dna.1, function(x) { IRanges(x=='-')})
  g = IRangesList(g)
  
  m = g[[l-1]]
  r = g[[l]]
  
  wid.m = width(m)
  wid.r = width(r)
  
  l.m = length(wid.m)
  l.r = length(wid.r)
  
  pos.m = start(m)
  pos.r = start(r)
  
  if(l.m > 0 & l.r == 0){
    df = data.frame("pos"=pos.m, "wid"=wid.m, "seq.len"=width(dna)[1])
  }else if(l.m == 0 & l.r > 0){
    df = data.frame("pos"=pos.r, "wid"=wid.r, "seq.len"=width(dna)[1])
  }else if(l.m > 0 & l.r > 0){
    dfm = data.frame("pos"=pos.m, "wid"=wid.m, "seq.len"=width(dna)[1])
    dfr = data.frame("pos"=pos.r, "wid"=wid.r, "seq.len"=width(dna)[1])
    df  = merge(dfm, dfr, all=TRUE)
  }else{
    return(NULL)
  }
  return(df)
}

#Distinguish phases
Record = function(PST){
  pos   = PST$pos
  phase = c()
  for (j in 1:length(pos)) {
    rem   = pos[j] %% 3
    if(rem==1){
      phase[j]=0
    }else if(rem==2){
      phase[j]=1
    }else{
      phase[j]=2
    }
  }
  
  PST.new = cbind(PST,phase)
  PST.new
}  

#Generate the phase dataframe
phase_gen = function(files){
  PST = c()
  for(i in 1:length(files)){
    phase_score = phasing(files[i])
    PST         = rbind(PST, phase_score)
    print(i)
  }
  Record(PST)
}


######################################
#Generate a type-N (0) /typeS (1) vector
eff_gen = function(files){
  typeNS = c()
  
  for(i in 1:length(files)){
    dna = readBStringSet(files[i])
    len = length(dna)
    M   = toString(dna[[len-1]])
    R   = toString(dna[[len]])
    
    dna.1 = str_split(as.character(dna),'')
    g     = lapply(dna.1, function(x) {IRanges(x=='-')})
    g     = IRangesList(g)
    
    m     = g[[len-1]]
    r     = g[[len]]
    wid.m = width(m)
    wid.r = width(r)
    pos.m = start(m)
    pos.r = start(r)
    l.m   = length(wid.m)
    l.r   = length(wid.r)
    
    if(l.m > 0){#mouse
      eff_M  = eff_phase(R,pos.m,wid.m,l.m)    
      typeNS = c(typeNS,eff_M) 
    }
    if(l.r > 0){#rat
      eff_R  = eff_phase(M,pos.r,wid.r,l.r)
      typeNS = c(typeNS,eff_R) 
    }
    print(i)
  }
  typeNS
}

#check if the gap is typeN (0) /typeS (1).
eff_phase = function(seq1, pos, wid, len){#seq1--reference
  seq1.char = str_split(seq1, "")[[1]]
  tyns = c()
  for(j in 1: len){
    if(pos[j] %% 3 == 1){#phase-0
      tyns[j] =  1 
    }else if(pos[j] %% 3 == 0){#phase-2
      pos.ori = pos[j]-2 
      unit.1  = substr(seq1, pos.ori, pos.ori+2)
      unit.2  = substr(seq1, pos.ori+wid[j], pos.ori+wid[j]+2)
      sub     = paste0(c(seq1.char[pos[j]-2], seq1.char[pos[j]-1], seq1.char[pos[j]+wid[j]]), collapse = "")
      sec     = codon[[which(sapply(codon, function(X){sub %in% X}))]]
      if(unit.1 %in% sec || unit.2 %in% sec){
        tyns[j] = 1
      }else{
        tyns[j] = 0
      }
    }else{#phase-1
      pos.ori = pos-1 
      unit.1  = substr(seq1, pos.ori, pos.ori+2)
      unit.2  = substr(seq1, pos.ori + wid[j], pos.ori+wid[j]+2)
      sub     = paste0(c(seq1.char[pos[j]-1], seq1.char[pos[j]+wid[j]], seq1.char[pos[j]+wid[j]+1]), collapse = "")
      sec     = codon[[which(sapply(codon, function(X){sub %in% X}))]]
      if(unit.1 %in% sec || unit.2 %in% sec){
        tyns[j] = 1
      }else{
        tyns[j] = 0
      }
    }
  }  
  tyns
}

#combine phase with typeN/S into a single table
sum_phase_type = function(inD){
  File.total = list.files(inD, full.names=TRUE)
  pha        = phase_gen(File.total)
  typeNS     = eff_gen(File.total)
  cbind(pha,typeNS)
}

##ggplot
gg_plot = function(df,wind){
  Pos    = round(df$pos/df$seq.len,3)
  newtab = cbind(Pos,df)
  g1  = ggplot(newtab,aes(x=factor(phase), y=Pos, color=factor(phase))) + 
    geom_dotplot(dotsize=.2, drop=T, binwidth=1/40, fill=NA, binaxis='y', stackdir='center') + 
    labs(x="Phase",y="Position",) + theme(axis.title=element_text(size=8))
  gg1 = g1 + coord_flip()
  
  g2  = ggplot(newtab,aes(x=factor(typeNS), y=Pos, color=factor(typeNS))) + 
    geom_dotplot(dotsize=.2, drop=T, binwidth=1/40, fill=NA, binaxis='y', stackdir='center') + 
    labs(x="typeN(0) vs typeS(1)",y="Position") + theme(axis.title=element_text(size=8))
  gg2 = g2 + coord_flip()
  
  gg    = ggarrange(gg1, gg2, labels=c("A","B"), ncol=1,nrow=2)
  gg    = annotate_figure(gg,top=text_grob(paste0("Window size = ",wind),face="bold",size=10))
  print(gg)
}






################################
#setwd("~/Dropbox (ASU)/Indel_project/test_human_mouse_rat")
#ouFig = "Results/Phase.mafft.pos.pdf"
main = function(ouFig){
  #Create a codon table
  codon <<- list (c("TTT","TTC"),
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
  
  dir1  = "Data_3/Mafft/mapped_cds"
  dir2  = "Data_6/Mafft/mapped_cds"
  dir3  = "Data_9/Mafft/mapped_cds"
  dir4  = "Data_12/Mafft/mapped_cds"
  
  sumtab1 = sum_phase_type(dir1)
  sumtab2 = sum_phase_type(dir2)
  sumtab3 = sum_phase_type(dir3)
  sumtab4 = sum_phase_type(dir4)
  
  ##PART II
  wind1=3; wind2=6; wind3=9; wind4=12
  
  
  pdf(ouFig,onefile=T, paper="letter", width=6.375, height=8.875)
  gg_plot(sumtab1,wind1)
  gg_plot(sumtab2,wind2)
  gg_plot(sumtab3,wind3)
  gg_plot(sumtab4,wind4)
  dev.off()
  
}
  
######################################
args = commandArgs(trailingOnly = TRUE)
main(args[1])