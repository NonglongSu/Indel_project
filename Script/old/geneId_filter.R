#Extract the CDS from ensembl API

# setwd("~/Dropbox (ASU)/Indel_project/test_human_mouse_rat/Script")

library(methods)
library(httr)
library(jsonlite)
library(xml2)
library(readr)
library(Biostrings)

# file.1 = "~/Downloads/trans_ID_NF.txt"
# file.2 = "../Raw_data/geneId.txt"
# ouFile="../Raw_data/geneId_new.txt"

main = function(file.1,file.2,ouFile){


  server = "https://rest.ensembl.org"

  gene.filter = read_tsv(file.1,col_names = "cds_NF")
  gene.list   = read_tsv(file.2,col_names = TRUE,quote = "\"")
  h.gene = gene.list[,1]

  canonical.lst = list()
  for (i in 1:nrow(h.gene)) {
    name = gsub('\"',"",h.gene[[1]][i],fixed=TRUE)
    ext  = paste("/lookup/id/",name,"?expand=1",sep="")
    r    = GET(paste(server,ext,sep=""),content_type("application/json"))
    stop_for_status(r)
  
    id  = fromJSON(toJSON(content(r)))
    id_ = id$Transcript[["id"]][id$Transcript$is_canonical==1]
  
    canonical.lst[i] = id_
}
    
  #canonical.lst = read_tsv("~/Downloads/gene.list.txt",col_names = FALSE)
  
  gene.lst = data.frame(matrix(unlist(gene.list),nrow=nrow(gene.list),byrow = FALSE))
  gene.lst["h.trans"] = canonical.lst
  
  
  gene.filter = unlist(gene.filter, use.names = FALSE)
  gene.lt     = gene.lst[!(gene.lst$h.trans %in% gene.filter),]
  gene.l      = gene.lt[-4]
  
  write.table(gene.l,ouFile,quote = FALSE,
              row.names=FALSE,col.names = FALSE,append = FALSE,sep = "\t")
}
args = commandArgs(trailingOnly=TRUE)
main(args[1],args[2],args[3])
