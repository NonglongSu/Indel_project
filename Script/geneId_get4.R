#extract the gene ID from ensembl database

# rm(list=ls())
# invisible(dev.off())
######################

library(BiocGenerics)
library(S4Vectors)
library(Biostrings)
library(biomaRt)
library(readr)
library(dplyr)

# attributes = c("ensembl_gene_id",
#                "cchok1gshd_homolog_ensembl_gene",
#                "cchok1gshd_homolog_orthology_type",
#                "mmusculus_homolog_ensembl_gene",
#                "mmusculus_homolog_orthology_type" ,
#                "rnorvegicus_homolog_ensembl_gene",
#                "rnorvegicus_homolog_orthology_type"
# )
# orth_filtered = orth %>% filter(`mmusculus_homolog_orthology_type`=="ortholog_one2one" & 
#                                   `rnorvegicus_homolog_orthology_type`=="ortholog_one2one" &
#                                   `cchok1gshd_homolog_orthology_type`=="ortholog_one2one")

# input  = "../Input/homo4.1.txt"
# output = "../test_human_mouse_rat/Raw_data.2/geneId.txt"
# setwd("~/Dropbox (ASU)/Indel_project/Script")

main = function(input, output){
  #Read input
  name = read_delim(input, "\t", col_names = FALSE)
  
  #Extract 
  ensembl = useMart("ensembl")
  data    = useDataset(paste0(toString(name[1, 1]), "_gene_ensembl"), mart = ensembl)
  
  # b = listAttributes(human)
  # b$name
  # a = listFilters(human)
  # a$name
  
  attributes = c("ensembl_gene_id")  
  for (i in 2:nrow(name)) {
    att1 = paste0(toString(name[i, 1]), "_homolog_ensembl_gene")
    att2 = paste0(toString(name[i, 1]), "_homolog_orthology_type")
    attributes = c(attributes, att1, att2)
  }
  
  filters = "with_ccds"
  orth    = getBM(attributes, filters, values = TRUE,
                  mart = data, uniqueRows = TRUE)
  
  #  element-wise comparison single &. 
  orth_filtered = orth[(orth[3] == "ortholog_one2one" & orth[5] == "ortholog_one2one"  & orth[7] == "ortholog_one2one"), ]
  orth_com      = orth_filtered[, c(1, 2, 4, 6)]
  
  write.table(orth_com, output, sep = '\t', row.names = F, col.names = F, quote = FALSE)
  
}

args = commandArgs(trailingOnly = TRUE)
main(args[1], args[2])





