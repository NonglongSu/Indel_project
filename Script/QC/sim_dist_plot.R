#Calculate the genetic distance between homologs (similarity %)

suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(readr))

#setwd("~/Dropbox (ASU)/Indel_project/test_human_mouse_rat/Raw_data")

Read_file = function(Files){
  sim   = list()
  sim2  = list()
  
  for(i in 1:length(Files)){
    sim[[i]] = perc(Files[i])
  }
  
  for (j in 1:length(sim[[1]])) {
    sim2[[j]] = unlist(lapply(sim, `[[`, j))
  }
  
  return(sim2)
}

#cal. similiarity
perc = function(filename){
  Seq = readBStringSet(filename, format="fasta")
  wid = width(Seq)
  
  specs = str_split(Seq,'')
  sim_vec = c()  
  for (i in 1:(length(specs)-1)) {
    for (j in (i + 1):length(specs)) {
      identical_pos = length(which(specs[[i]][specs[[i]] == specs[[j]]] != '-'))
      shared_gaps   = length(which(specs[[i]][specs[[i]] == specs[[j]]] == '-'))
      total_gaps    = length(which(specs[[i]] =='-')) + length(which(specs[[j]] =='-'))
      sim.score     = 100 * (identical_pos) / (wid[1] - (total_gaps-shared_gaps))
      sim_vec       = c(sim_vec, sim.score)
    }
  }
  return(sim_vec)
}

#plot 
plot_Data = function(color, data, name.comb, x_note){
  
  # hist(data[[1]], xlab=x_note,col=color[1],main=NULL,breaks=40,freq=F)
  # hist(data[[2]],col=color[2],add=TRUE,breaks=40,freq=F)
  # hist(data[[3]],col=color[3],add=TRUE,breaks=40,freq=F)
  
  plot(density(data[[1]]),  col=color[1], lwd=3, main="",ylim=c(0,0.2),
       xlab=x_note)
  lines(density(data[[2]]), col=color[2], lwd=3)
  lines(density(data[[3]]), col=color[3], lwd=3)
  
  legend("topleft", legend=name.comb, col=color, lty=1, bty="n", cex=1)
}


##########################################
# dir1   = "Alignments/mafft_out"
# dir2   = "mapped_cds_mafft"
# inFile = "../Input_species/homo3.1.txt"
# ouFile = "QC/low_sim.txt"
main = function(dir1, dir2, inFile, ouFig){
  
  File1 = list.files(dir1, full.names=TRUE)
  File2 = list.files(dir2, full.names=TRUE)
  data1 = Read_file(File1)
  data2 = Read_file(File2)
  
  #Remove data with missing values.
  data1.1 = lapply(data1, function(x){x[!is.na(x)]})
  data2.1 = lapply(data2, function(x){x[!is.na(x)]})
  
  #adjust legend
  name  = read_delim(inFile,"\t",col_names=FALSE,show_col_types=F)
  name2 = name[[3]]
  l     = length(name2)
  
  name.comb = c()
  for (i in 1:(l-1)) {
    for (j in (i+1):l) {
      name.comb = c(name.comb, paste0(c(name2[i],name2[j]),collapse = "-"))
    }
  }
  
  #Plot the density distribution
  color   = c("#56B4E9","#009E73","#CC79A7")
  x_note1 = "Sequence identity AA(%)"
  x_note2 = "Sequence identity CDS(%)"
  
  pdf(ouFig)
  par(mfrow = c(1,2))
  
  plot_Data(color, data1.1, name.comb, x_note1)
  plot_Data(color, data2.1, name.comb, x_note2)
  
  invisible(dev.off())
  embedFonts(ouFig)
}



##########################################
args = commandArgs(trailingOnly=TRUE)
main(args[1], args[2], args[3], args[4])


