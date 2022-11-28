#Add stopcodon (TAG) to the end of the sequence. 

#setwd("~/Dropbox (ASU)/Indel_project/Script")
library(readr)
library(Biostrings)
library(stringr)
library(seqinr)


#file = "../test_human_mouse_rat/Raw_data/QC/noStop.txt"
#dir  = "../test_human_mouse_rat/Raw_data/cds_seq/"

stopCodon_add = function(file, dir){
  
  File = read_delim(file, "\t", col_names = FALSE)

  for (i in 1:nrow(File)) {
    filename = paste0(dir, File[[i, 1]], ".fa")
    cds = readDNAStringSet(filename, format = "fasta")
    wid = width(cds)
    new_cds = list()
    for(j in 1:length(cds)){
      dna = toString(cds[[j]]) 
      Stop = str_sub(cds[[j]], start = wid[j] - 2, end = wid[j])
      if(Stop != "TAG" & Stop != "TGA" & Stop != "TAA"){
        dna = gsub('^(.*)$', '\\1TAG\\2', cds[[j]])
      }
      new_cds[[j]] = dna
    }
    name = names(cds)
    write.fasta(sequences = new_cds, names = name, nbchar = 80,
                open = "w",as.string = TRUE, file.out = filename)
  }
  
}

args = commandArgs(trailingOnly = TRUE)
stopCodon_add(args[1], args[2])


############################
#file_name = "../Raw_data/haha.fa"
# cds <- readDNAStringSet(inFile,format="fasta")


