library(Biostrings)
library(stringr)

#setwd("~/Dropbox (ASU)/Indel_project/test_human_mouse_rat/Script")

# file.1 = "../Data/mapped_cds_Out/ENSG00000090905.fa"
# file.2 = "../Data/mapped_cds_plus/ENSG00000090905.fa"
# file.3 = "../Raw_data/mapped_cds_mafft_Tmp/ENSG00000090905.fa"

phasing = function(file.1,file.2,file.3){
  
  dna.1   = readBStringSet(file.1,format = "fasta")
  dna.1   = str_split(as.character(dna.1), '')
  
  dna.2   = readBStringSet(file.2,format = "fasta")
  dna.2   = str_split(as.character(dna.2), '')
  
  dna.3   = readBStringSet(file.3,format = "fasta")
  dna.3   = str_split(as.character(dna.3), '')
  
  g.1 = lapply(dna.1, function(x) { IRanges(x == '-')})
  g.1 = IRangesList(g.1)
  
  g.2 = lapply(dna.2, function(x) { IRanges(x == '-')})
  g.2 = IRangesList(g.2)
  
  m = g.1[[1]]
  r = g.1[[2]]
  
  wid.m = width(m)
  wid.r = width(r)
  
  l.m = length(wid.m)
  l.r = length(wid.r)
  
  pos.m = start(m)
  pos.r = start(r)
  
  df.m = data.frame("pos"=pos.m,"wid"=wid.m)
  df.r = data.frame("pos"=pos.r,"wid"=wid.r)
  
  col.name = c("insert_1","insert_2","delete_1","delete_2")
  
  if(l.m > 0 & l.r == 0){
   iw.1 = indel_watch_1(dna.3,g.2,l.m)
   df.m[col.name] = c(iw.1[[1]],iw.1[[2]],iw.1[[3]],iw.1[[4]])
   return(df.m)
  }else if(l.m == 0 & l.r > 0){
    iw.2 = indel_watch_2(dna.3,g.2,l.r)
    df.r[col.name] = c(iw.2[[1]],iw.2[[2]],iw.2[[3]],iw.2[[4]])
    return(df.r)
  }else if(l.m > 0 & l.r > 0){
    iw.1 = indel_watch_1(dna.3,g.2,l.m)
    iw.2 = indel_watch_2(dna.3,g.2,l.r)
    df.m[col.name] = c(iw.1[[1]],iw.1[[2]],iw.1[[3]],iw.1[[4]])
    df.r[col.name] = c(iw.2[[1]],iw.2[[2]],iw.2[[3]],iw.2[[4]])
    df = merge(df.m,df.r,all=TRUE)
    return(df)
  }else{
    return(NULL)
  }

}
  

indel_watch_1 = function(dna.3,g.2,l.m) {
  m     = g.2[[1]] 
  pos.m = start(m) 
  
  insert_1 = c()
  insert_2 = c()
  delete_1 = c()
  delete_2 = c()
  
  for(i in 1:l.m) {
    insert.1 = delete.1 = insert.2 = delete.2 = 0
    id = pos.m[i]
      if(dna.3[[1]][id] == '-'){
        insert.2 = 1
       }else{
        delete.1 = 1
      }
      insert_1[i] = insert.1
      insert_2[i] = insert.2
      delete_1[i] = delete.1
      delete_2[i] = delete.2
    }
  indel_total = list(insert_1,insert_2,delete_1,delete_2)
  return(indel_total)
}

indel_watch_2 = function(dna.3,g.2,l.r) {
  r     = g.2[[2]]
  pos.r = start(r)
  
  insert_1 = c()
  insert_2 = c()
  delete_1 = c()
  delete_2 = c()

  for(i in 1:l.r) {
    insert.1 = delete.1 = insert.2 = delete.2 = 0
    id = pos.r[i]
      if(dna.3[[1]][id] == '-'){
        insert.1 = 1
       }else{
        delete.2 = 1
      }
      insert_1[i] = insert.1
      insert_2[i] = insert.2
      delete_1[i] = delete.1
      delete_2[i] = delete.2
    }
  indel_total = list(insert_1,insert_2,delete_1,delete_2)
  return(indel_total)
  
}


# dir1 = "/home/zzhu49/Dropbox (ASU)/Indel_project/test_human_mouse_rat/Data/mapped_cds_Out"
# dir2 = "/home/zzhu49/Dropbox (ASU)/Indel_project/test_human_mouse_rat/Data/mapped_cds_plus"
# dir3 = "/home/zzhu49/Dropbox (ASU)/Indel_project/test_human_mouse_rat/Raw_data/mapped_cds_mafft_Tmp/"


main=function(dir1,dir2,dir3,ouFile){
  
  File1.total = list.files(dir1, full.names = TRUE)
  File2.total = list.files(dir2, full.names = TRUE)
  
  #Grab outgroup homologs. 
  nameList    = basename(File1.total)
  File3.total = c()
  for (i in 1:length(nameList)) {
    file.path       = paste0(dir3,nameList[i])
    File3.total[i]  = file.path
  }
  
  PST = c()
  for(i in 1:length(File1.total)){
    print(nameList[i])
    phase_score = phasing(File1.total[i],File2.total[i],File3.total[i])
    PST = rbind(PST,phase_score)
    
  }
  
  write.table(PST,file = ouFile,sep = "\t",append = FALSE,quote = FALSE,row.names = FALSE,col.names = TRUE)
  
}
  
args = commandArgs(trailingOnly = TRUE)
main(args[1],args[2],args[3],args[4])
  

