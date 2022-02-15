#Calculate the genetic distance between homologs (Hamming distance)(Quality control)

library(Biostrings)
library(BiocGenerics)
library(parallel)
library(S4Vectors)


#setwd("~/Dropbox (ASU)/Indel_project/Script")
#filename="../Raw_data.2.outgroups/Tmp/ENSG00000158796.fa"

perc = function(filename){
  y   = readBStringSet(filename,format = "fasta")
  wid = width(y)
  l   = length(wid)
  
  spec.1  = as.vector(y[[1]],mode = "any")
  spec.2  = as.vector(y[[2]],mode = "any")
  spec.3  = as.vector(y[[3]],mode = "any")
  spec.4  = as.vector(y[[4]],mode = "any")
  
  #Cal the identical number of bases between any two species(gaps excluded)
  identical_pos1 = length(which(spec.1[spec.1 == spec.2] != '-'))
  identical_pos2 = length(which(spec.1[spec.1 == spec.3] != '-'))
  identical_pos3 = length(which(spec.1[spec.1 == spec.4] != '-'))
  identical_pos4 = length(which(spec.2[spec.2 == spec.3] != '-'))
  identical_pos5 = length(which(spec.2[spec.2 == spec.4] != '-'))
  identical_pos6 = length(which(spec.3[spec.3 == spec.4] != '-'))
  
  id_pos = c(identical_pos1,identical_pos2,identical_pos3,identical_pos4,identical_pos5,identical_pos6)
  
  #Cal shared gaps between any two species
  shared_gap1 = length(which(spec.1[spec.1 == spec.2] == '-'))
  shared_gap2 = length(which(spec.1[spec.1 == spec.3] == '-'))
  shared_gap3 = length(which(spec.1[spec.1 == spec.4] == '-'))
  shared_gap4 = length(which(spec.2[spec.2 == spec.3] == '-'))
  shared_gap5 = length(which(spec.2[spec.2 == spec.4] == '-'))
  shared_gap6 = length(which(spec.3[spec.3 == spec.4] == '-'))
  
  #Cal total gaps of each species
  l1 = length(which(spec.1 =='-'))
  l2 = length(which(spec.2 =='-'))
  l3 = length(which(spec.3 =='-'))
  l4 = length(which(spec.4 =='-'))
  
  #Cal the similarity of two seqs
  combinations = factorial(l)/(factorial(l-2)*factorial(2))
  comb.l = rep(wid[1],combinations)
  
  aligned.all = comb.l - c((l1+l2-shared_gap1),(l1+l3-shared_gap2),(l1+l4-shared_gap3),
                           (l2+l3-shared_gap4),(l2+l4-shared_gap5),(l3+l4-shared_gap6))
  similarity = 100*(id_pos)/(aligned.all)
  return (similarity)
}

distance_plot= function(dir,ouFile,ouFig){
  
  inFile = list.files(dir, full.names = TRUE)
  
  spec.1_2 = c()
  spec.1_3 = c()
  spec.1_4 = c()
  spec.2_3 = c()
  spec.2_4 = c()
  spec.3_4 = c()
  
  for(i in 1:length(inFile)){
    sim = perc(inFile[i])
    spec.1_2[i] = sim[1]
    spec.1_3[i] = sim[2]
    spec.1_4[i] = sim[3]
    spec.2_3[i] = sim[4]
    spec.2_4[i] = sim[5]
    spec.3_4[i] = sim[6]
  }
  
  # Check if the smallest distance always exist between mouse and rat.
  ab_spec = spec.1_2[spec.1_2 >spec.3_4 | spec.1_3>spec.3_4 | spec.1_4>spec.3_4 | 
                     spec.2_3>spec.3_4  | spec.2_4>spec.3_4 ]
  
  
  # match file orders
  v1 = match(ab_spec,spec.1_2)
  for(i in 1:length(v1)){
    cat(basename(inFile[v1[i]]),file = ouFile,append = TRUE,sep = "\n")
  }
  
  
  # Plot the similiary distribution
  #outFig="../Figure/distance_of_orthologs_cds.pdf"
  pdf(ouFig)
  color = c("red","green","blue","black","yellow","gray")
  
  plot(density(spec.3_4),xlab="percent sequence identity",ylab="freq",
       main = " Kernal density plot of human-hamster-mouse-rat similary",col=color[1])
  lines(density(spec.1_2),col=color[2])
  lines(density(spec.1_3),col=color[3])
  lines(density(spec.1_4),col=color[4])
  lines(density(spec.2_3),col=color[5])
  lines(density(spec.2_4),col=color[6])
  
  legend("topleft",legend = c("mouse-rat","human-hamster","human-mouse","human-rat","hamster-mouse","hamster-rat"),
         col=color,lty=1,bty="n",cex=1.2)
  
  invisible(dev.off())
  embedFonts(ouFig)
  
}

args = commandArgs(trailingOnly = TRUE)
distance_plot(args[1],args[2],args[3])
