#Generate pair species only. 

library(seqinr)

pair = function(inFile){
  
  #inFile = file.path("..","Raw_data",inFile)
  #inFile = file.path("..","Raw_data.2.outgroups",inFile)
  
  data = read.fasta(inFile,seqtype = "DNA",as.string = FALSE,
                     set.attributes = FALSE,forceDNAtolower = FALSE)
  
  data.new = tail(data,2)
  name = names(data.new)
  
  write.fasta(sequences = data.new, names = name, nbchar=80,
              open = "w",as.string = FALSE, file.out = inFile)
}

args = commandArgs(trailingOnly = TRUE)
pair(args[1])
