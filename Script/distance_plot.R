#Calculate the genetic distance between homologs (Hamming distance)(Quality control)

library(Biostrings)
library(BiocGenerics)
library(readr)
library(parallel)
library(S4Vectors)

#setwd("~/Dropbox (ASU)/Indel_project/Script")

perc = function(filename){
  
  Seq = readBStringSet(filename, format = "fasta")
  wid = width(Seq)
  l   = length(wid)
  
  specs = lapply(Seq, function(x){as.vector(x)})
  sim_vec = c()  
  for (i in 1:(length(specs) - 1)) {
    for (j in (i + 1):length(specs)) {
      identical_pos = length(which(specs[[i]][specs[[i]] == specs[[j]]] != '-'))
      shared_gaps   = length(which(specs[[i]][specs[[i]] == specs[[j]]] == '-'))
      total_gaps    = length(which(specs[[i]] =='-')) + length(which(specs[[j]] =='-'))
      similarity    = 100 * (identical_pos) / (wid[1] - (total_gaps - shared_gaps))
      sim_vec       = c(sim_vec, similarity)
    }
  }
  return(sim_vec)
}

Read_file = function(dir){
  
  Files = list.files(dir, full.names = TRUE)
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

plot_Data = function(color, data, title, comb, x_note){
  plot(density(data[[length(data)]]), xlab = x_note, ylab = "Density",
       main = paste0("Kernal density plot of ", title, " similarity"), col = color[length(color)])
  
  for (i in 1: (length(comb) - 1)) {
    lines(density(data[[i]]), col = color[i])
  }
  legend("topleft", legend = comb,
         col = color, lty = 1, bty = "n", cex = 0.8)
}


# dir1 = "../test_human_mouse_rat/Raw_data_2/Alignments/mafft_out"
# dir2 = "../test_human_mouse_rat/Raw_data_2/mapped_cds_mafft"
# inFile = "../Input/homo4.1.txt"
distance_plot = function(dir1, dir2, inFile, ouFig){
  
  data1 = Read_file(dir1)
  data2 = Read_file(dir2)
  
  #Remove data with missing values.
  data1.1 = lapply(data1, function(x){x[!is.na(x)]})
  data2.1 = lapply(data2, function(x){x[!is.na(x)]})
  
  #Working on legend
  name  = read_delim(inFile, "\t", col_names = FALSE )
  name2 = name[, 3][[1]]
  l     = length(name2)
  
  name.comb = c()
  for (i in 1:(l-1)) {
    for (j in (i+1):l) {
      name.comb = c(name.comb, paste0(c(name2[i], name2[j]), collapse = "-"))
    }
  }
  name.title = paste0(name2, collapse = '-')
  
  # # Check if the smallest distance always exist between mouse and rat.
  # ab_spec = spec.1_2[spec.1_2 > spec.2_3 | spec.1_3 > spec.2_3]
  # 
  # # match file orders
  # v1 = match(ab_spec, spec.1_2)
  # for(i in 1:length(v1)){
  #   cat(basename(Files[v1[i]]), file = ouFile, append = TRUE, sep = "\n")
  # }
  
  
  # Plot the density distribution
  pdf(ouFig)
  par(mfrow = c(2, 1))
  
  color   = rainbow(length(name.comb))
  x_note1 = "Sequence identity AA(%)"
  x_note2 = "Sequence identity CDS(%)"
  
  plot_Data(color, data1.1, name.title, name.comb, x_note1)
  plot_Data(color, data2.1, name.title, name.comb, x_note2)
  
  invisible(dev.off())
  embedFonts(ouFig)
}


args = commandArgs(trailingOnly = TRUE)
distance_plot(args[1], args[2], args[3], args[4])


