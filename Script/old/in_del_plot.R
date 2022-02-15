library(readr)

#setwd("~/Dropbox (ASU)/Indel_project/test_human_mouse_rat/Script")


Record = function(PST){
  if(is.null(PST)){return(PST)}
  rem   = PST%%3
  phase.0 = length(which(rem==1))
  phase.1 = length(which(rem==2)) 
  phase.2 = length(which(rem==0)) 
  
  df.phase = data.frame(Phase_0 = phase.0, Phase_1 = phase.1, Phase_2 = phase.2)
  return(df.phase)
}  


#filename = "../Data/indelCount.mafft.txt"

main = function(filename,ouFile,ouFig){
  
  data = read_delim(filename,"\t", col_names = TRUE)
  
  m.pos = data[data$insert_1==1 | data$delete_1==1,]
  r.pos = data[data$insert_2==1 | data$delete_2==1,]
  
  #Profiling the phases
  m.pos.3  = m.pos[m.pos$wid==3,]$pos
  m.pos.6  = m.pos[m.pos$wid==6,]$pos
  m.pos.9  = m.pos[m.pos$wid==9,]$pos
  m.pos.12 = m.pos[m.pos$wid==12,]$pos
  
  r.pos.3  = r.pos[r.pos$wid==3,]$pos
  r.pos.6  = r.pos[r.pos$wid==6,]$pos
  r.pos.9  = r.pos[r.pos$wid==9,]$pos
  r.pos.12 = r.pos[r.pos$wid==12,]$pos
  
  m.phase.3  = Record(m.pos.3)
  m.phase.6  = Record(m.pos.6)
  m.phase.9  = Record(m.pos.9)
  m.phase.12 = Record(m.pos.12)
  
  r.phase.3  = Record(r.pos.3)
  r.phase.6  = Record(r.pos.6)
  r.phase.9  = Record(r.pos.9)
  r.phase.12 = Record(r.pos.12)
  
  #Put extra insertion/deleteion info.
  m.ins.3  = length(data[data$insert_1==1 & data$wid==3,]$insert_1)
  m.del.3  = length(data[data$delete_1==1 & data$wid==3,]$delete_1)
  m.ins.6  = length(data[data$insert_1==1 & data$wid==6,]$insert_1)
  m.del.6  = length(data[data$delete_1==1 & data$wid==6,]$delete_1)
  m.ins.9  = length(data[data$insert_1==1 & data$wid==9,]$insert_1)
  m.del.9  = length(data[data$delete_1==1 & data$wid==9,]$delete_1)
  m.ins.12 = length(data[data$insert_1==1 & data$wid==12,]$insert_1)
  m.del.12 = length(data[data$delete_1==1 & data$wid==12,]$delete_1)
  
  r.ins.3  = length(data[data$insert_2==1 & data$wid==3,]$insert_2)
  r.del.3  = length(data[data$delete_2==1 & data$wid==3,]$delete_2)
  r.ins.6  = length(data[data$insert_2==1 & data$wid==6,]$insert_2)
  r.del.6  = length(data[data$delete_2==1 & data$wid==6,]$delete_2)
  r.ins.9  = length(data[data$insert_2==1 & data$wid==9,]$insert_2)
  r.del.9  = length(data[data$delete_2==1 & data$wid==9,]$delete_2)
  r.ins.12 = length(data[data$insert_2==1 & data$wid==12,]$insert_2)
  r.del.12 = length(data[data$delete_2==1 & data$wid==12,]$delete_2)
  
  
  #Generate the table
  m.Phase.df = rbind(m.phase.3,m.phase.6,m.phase.9,m.phase.12)
  
  m.insert = c(m.ins.3,m.ins.6,m.ins.9,m.ins.12)
  m.delete = c(m.del.3,m.del.6,m.del.9,m.del.12)
  
  m.Phase.df["insertion"] = m.insert
  m.Phase.df["deletion"]  = m.delete
  
  m.Phase.df = rbind(m.Phase.df,colSums(m.Phase.df))
  
  m.gap      = c('3','6','9','12',"T")
  m.Phase.df["gap.length"] = m.gap
  
  #Second species
  r.Phase.df = rbind(r.phase.3,r.phase.6,r.phase.9,r.phase.12)
  
  r.insert = c(r.ins.3,r.ins.6,r.ins.9,r.ins.12)
  r.delete = c(r.del.3,r.del.6,r.del.9,r.del.12)
  
  r.Phase.df["insertion"] = r.insert
  r.Phase.df["deletion"]  = r.delete
  
  r.Phase.df = rbind(r.Phase.df,colSums(r.Phase.df))
  
  r.gap      = c('3','6','9','12',"T")
  r.Phase.df["gap.length"] = r.gap
  
  #ouFile = "../Data/indelTable.mafft.txt"
  
  write.table(m.Phase.df,file = ouFile,sep = "\t",append = FALSE,quote = FALSE,row.names = FALSE,col.names = TRUE)
  cat(file = ouFile,sep="\n",append = TRUE)
  write.table(r.Phase.df,file = ouFile,sep = "\t",append = TRUE,quote = FALSE,row.names = FALSE,col.names = TRUE)
  
  
  #dev.off(dev.list())
  #ouFig="../Figure/Phase.sep.pdf"
  pdf(ouFig)
  
  name   = colnames(m.Phase.df)[1:3]
  
  m.p.3  = as.numeric(m.Phase.df[1,1:3])
  m.p.6  = as.numeric(m.Phase.df[2,1:3])
  m.p.9  = as.numeric(m.Phase.df[3,1:3])
  m.p.12 = as.numeric(m.Phase.df[4,1:3])
  m.p.T  = as.numeric(m.Phase.df[5,1:3])
  
  r.p.3  = as.numeric(r.Phase.df[1,1:3])
  r.p.6  = as.numeric(r.Phase.df[2,1:3])
  r.p.9  = as.numeric(r.Phase.df[3,1:3])
  r.p.12 = as.numeric(r.Phase.df[4,1:3])
  r.p.T  = as.numeric(r.Phase.df[5,1:3])
  
  par(mfrow=c(3,2))
  
  #Sometimes if all values in P._ is zero, then errors shows.  
  pie(m.p.3,labels=name,main  = "Indel-phases of gap length 3",col=rainbow(length(m.p.3)))
  pie(m.p.6,labels=name,main  = "Indel-phases of gap length 6",col=rainbow(length(m.p.6)))
  pie(m.p.9,labels=name,main  = "Indel-phases of gap length 9",col=rainbow(length(m.p.9)))
  pie(m.p.12,labels=name,main = "Indel-phases of gap length 12",col=rainbow(length(m.p.12)))
  pie(m.p.T,labels=name,main  = "Indel-phases of all gap lengths",col=rainbow(length(m.p.T)))
  
  mtext("Indel phases distribution of 1st focal species", side = 1, outer=TRUE, cex=1, line=-2)
  
  plot.new()
  par(mfrow=c(3,2))

  pie(r.p.3,labels=name,main  = "Indel-phases of gap length 3",col=rainbow(length(r.p.3)))
  pie(r.p.6,labels=name,main  = "Indel-phases of gap length 6",col=rainbow(length(r.p.6)))
  pie(r.p.9,labels=name,main  = "Indel-phases of gap length 9",col=rainbow(length(r.p.9)))
  pie(r.p.12,labels=name,main = "Indel-phases of gap length 12",col=rainbow(length(r.p.12)))
  pie(r.p.T,labels=name,main  = "Indel-phases of all gap lengths",col=rainbow(length(r.p.T)))
  
  mtext("Indel phases distribution of 2nd focal species",side = 1, outer=TRUE, cex=1, line=-2)

}

args = commandArgs(trailingOnly = TRUE)
main(args[1],args[2],args[3])

