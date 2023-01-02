#Extract all CDS sequence from 1:1:1 orthologous ID of dropophila species

#rm(list = ls())
#setwd("~/Dropbox (ASU)/Indel_project/test_drosophila/Script")

library(Biostrings)
library(readr)
library(dplyr)

file ="../Raw_data/FlyB_Id.txt"
file1="../Raw_data/Dyak.fasta"
file2="../Raw_data/Dsim.fasta"
file3="../Raw_data/Dmel.fasta"

main = function(file,file1,file2,file3){
  
  #Pull out all orthologous ID
  data_T = read_delim(file,"\t", col_names = TRUE)
  
  Id_Dyak_toMatch = data_T$Dyak_ID
  Id_Dsim_toMatch = data_T$Dsim_ID
  Id_Dmel_toMatch = data_T$Dmel_ID
  
  #Pull out whole IDs from separate fasta files.
  cds_Dyak = readDNAStringSet(file1,format="fasta")
  cds_Dsim = readDNAStringSet(file2,format="fasta")
  cds_Dmel = readDNAStringSet(file3,format="fasta")
  
  name_Dyak = names(cds_Dyak)
  name_Dsim = names(cds_Dsim)
  name_Dmel = names(cds_Dmel)
  
  
  #grab the pattern between.
  name_Dyak_Id = sub(".*parent= *(.*?) *,.*", "\\1", name_Dyak)
  name_Dsim_Id = sub(".*parent= *(.*?) *,.*", "\\1", name_Dsim)
  name_Dmel_Id = sub(".*parent= *(.*?) *,.*", "\\1", name_Dmel)
  
  #Pattern match test
  #Remove duplicate Id in fasta files
  #Id_Dyak_matched = unique (grep(paste(Id_Dyak_toMatch[1:10000],collapse="|"), name_Dyak_Id, value=TRUE))
  
  Id_Dyak_matched = Id_Dyak_toMatch[Id_Dyak_toMatch %in% name_Dyak_Id]
  Id_Dsim_matched = Id_Dsim_toMatch[Id_Dsim_toMatch %in% name_Dsim_Id]
  Id_Dmel_matched = Id_Dmel_toMatch[Id_Dmel_toMatch %in% name_Dmel_Id]
  #All of them are matched "Perfecto except for Id_Dmel. "
  
  #We have to createa new dataset
  dataT = data_T %>% filter(data_T$Dmel_ID %in% Id_Dmel_matched)
  
  
  #Resonstruction
  names(cds_Dyak) = dataT$Dyak_ID
  names(cds_Dsim) = dataT$Dsim_ID
  names(cds_Dmel) = dataT$Dmel_ID
  
  cds_Dyak_new = cds_Dyak[names(cds_Dyak) %in% dataT$Dyak_ID]
  cds_Dsim_new = cds_Dsim[names(cds_Dsim) %in% dataT$Dsim_ID]
  cds_Dmel_new = cds_Dmel[names(cds_Dmel) %in% dataT$Dmel_ID]
  
  
  nameId = names(cds_Dyak_new)
  output = "../Raw_data/cds_seq/"
  
  for(i in 1:length(cds_Dyak_new)){
    seq = c(cds_Dyak_new[i],cds_Dsim_new[i],cds_Dmel_new[i])
    writeXStringSet(seq,paste0(output,nameId[i],".fa"), append=FALSE,
                    compress=FALSE, compression_level=NA, format="fasta")
  }
}



args = commandArgs(trailingOnly = TRUE)
if (length(args)< 4) {
  stop("At least four arguments must be supplied.\n", call.=FALSE)
}
main(args[1],args[2],args[3],args[4])

