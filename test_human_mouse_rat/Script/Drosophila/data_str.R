#Get all 1:1:1 orthologous IDs of Dropophila species

#rm(list = ls())
#setwd("~/Dropbox (ASU)/Indel_project/test_drosophila/Script")
#filename="../Raw_data/Table_Dmel.txt"

library(dplyr)
library(purrr)
library(tidyverse)
library(readr)

feat=function(x){
  feature = c("Dsim","Dyak")
  all(feature %in% x)
}

main = function(filename){
  
  output = "../Raw_data/FlyB_Id.txt"
  data = read_delim(filename,"\t", col_names = TRUE)
  
  #1:1:1 ortholog filtering.
  dataA = data %>%
    group_by(`## FBgn_ID`) %>%
    filter(Ortholog_GeneSymbol=="Dsim" | Ortholog_GeneSymbol=="Dyak")
  
  #Filter based on occurrence.
  dataB = dataA %>% group_by(`## FBgn_ID`) %>%
    filter(n() == 2)
  
  #Ensure if only if Dyak,Dsim exist in each group.
 dataC = dataB %>%
         group_by(`## FBgn_ID`) %>%
         summarize(Feat = feat(Ortholog_GeneSymbol)) %>%
         filter(Feat == TRUE)
 
 dataD = dataB[dataB$`## FBgn_ID` %in% dataC$`## FBgn_ID`,]
 
 #Reframe the data
 Dsim = dataD %>% group_by(`## FBgn_ID`) %>% filter(Ortholog_GeneSymbol=="Dsim")
 Dyak = dataD %>% group_by(`## FBgn_ID`) %>% filter(Ortholog_GeneSymbol=="Dyak")
 
 
 new_data = data.frame("Dyak_ID" = Dyak$Ortholog_FBgn_ID, "Dsim_ID" = Dsim$Ortholog_FBgn_ID, "Dmel_ID" = Dsim$`## FBgn_ID`)
 
 #Remove the duplicated Id
 n = new_data[!duplicated(new_data$Dyak_ID),]
 
 write.table(n,file=output,
             sep='\t',quote = FALSE, row.names=F,col.names=TRUE)

}

args = commandArgs(trailingOnly = TRUE)
main(args[1])


