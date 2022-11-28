#extract the gene ID from ensembl API
######################
# .libPaths()
# [1] "/home/zzhu49/R/x86_64-pc-linux-gnu-library/3.4"
# [2] "/usr/local/lib/R/site-library"                 
# [3] "/usr/lib/R/site-library"                       
# [4] "/usr/lib/R/library"  

library(BiocGenerics)
library(S4Vectors)
library(Biostrings)
library(biomaRt)
library(readr)
library(dplyr)
# setwd("Dropbox (ASU)/Indel_project/Script")

# human = useDataset("hsapiens_gene_ensembl", mart = ensembl)
# attributes = c("ensembl_gene_id",
#                "mmusculus_homolog_ensembl_gene",
#                "mmusculus_homolog_orthology_type" ,
#                "rnorvegicus_homolog_ensembl_gene",
#                "rnorvegicus_homolog_orthology_type"
# )
# orth_filtered = orth %>% filter(`mmusculus_homolog_orthology_type` == "ortholog_one2one" &
#                                   `rnorvegicus_homolog_orthology_type` == "ortholog_one2one" )


# input  = "../Input/homo3.1.txt"
# output = "../test_human_mouse_rat/Raw_data/geneId.txt"

main =function(input, output){
  #Read input
  name = read_delim(input, "\t", col_names = FALSE)
  
  #Extract 
  ensembl = useMart("ensembl")
  data    = useDataset(paste0(toString(name[1, 1]), "_gene_ensembl"), mart = ensembl)l
  
  # b = listAttributes(human)
  # b$name
  # a = listFilters(human)
  # a$name
    
  attributes = c("ensembl_gene_id")  
  for (i in 2:nrow(name)) {
    att1 = paste0(toString(name[i, 1]), "_homolog_ensembl_gene")
    att2 = paste0(toString(name[i, 1]), "_homolog_orthology_type")
    attributes = c(attributes ,att1, att2)
  }

  filters = "with_ccds"
  orth    = getBM(attributes, filters, values = TRUE,
                   mart = data, uniqueRows = TRUE)
  
  # Element-wise comparison with &. 
  orth_filtered = orth[(orth[3] == "ortholog_one2one" & orth[5] == "ortholog_one2one"), ]
  orth_com      = orth_filtered[, c(1, 2, 4)]
  
  # Write data
  write.table(orth_com, output, sep = '\t', row.names = F, col.names = F, quote = FALSE)
  
}

args = commandArgs(trailingOnly = TRUE)
main(args[1], args[2])





