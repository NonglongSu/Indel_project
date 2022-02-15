#Extract the CDS from ensembl API.

#setwd("~/Dropbox (ASU)/Indel_project/Script")
library(methods)
library(httr)
library(jsonlite)
library(xml2)
library(readr)
library(Biostrings)
library(stringr)


# file = "../test_human_mouse_rat/Raw_data_2/geneId.txt"
# n    = "ENSG00000128815"
# output   = "../test_human_mouse_rat/Raw_data/cds_seq/"    #make cds
# output.1 = "../Raw_data.2.outgroups/cds_seq/Tmp1/"        #make cds.1
# output.2 = "../Raw_data.2.outgroups/cds_seq/Tmp2/"        #make cds.2

#name="ENSG00000158796"
main = function(n, file, output){
  
  server = "http://rest.ensembl.org"
  geneId = read_delim(file, "\t", col_names = FALSE)  
  
  n          = str_extract(basename(n), "[^.]+")
  gene.stem  = geneId[which(geneId[[1]] == n), ]
  
  cds.list = list()
  for (i in 1:length(gene.stem)) {
      name = gene.stem[[i]]
      ext  = paste("/lookup/id/", name, "?expand=1",sep = "")
      r    = GET(paste(server, ext, sep = ""), content_type("application/json"), timeout(30))
      stop_for_status(r)
      id   = fromJSON(toJSON(content(r)))
      id_  = id$Transcript[["id"]][id$Transcript$is_canonical==1]
      
      #quality control
      if(is.null(id_)){
        quit()
      }
      
      ext_1 = paste("/sequence/id/",id_,"?type=cds",sep = "")
      r_1   = GET(paste(server, ext_1, sep = ""), content_type("text/x-fasta"))
      stop_for_status(r_1)
      cds_  = (content(r_1))
      
      #quality control
      if(is.null(cds_)){
        quit()
      }
      
      cds.list = c(cds.list, cds_)
  }

  write.table(cds.list, paste0(output, n, ".fa"), quote = FALSE,
              row.names = FALSE, col.names = FALSE, append = FALSE, sep = "\n")
  
}

args = commandArgs(trailingOnly = TRUE)
main(args[1], args[2], args[3])
  


