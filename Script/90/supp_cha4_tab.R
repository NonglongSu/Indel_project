suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(readr))  
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(ggplot2))  
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(jsonlite))


#setwd("~/Dropbox (ASU)/Indel_project/test_90_species")

#Phase-Proportion-table
#####################################
inD1  = "Results/phases"
ouF1  = "Results/supp4/phase_prop_all.csv"
ouFig1= "Results/supp4/phase_prop.pdf"

Files  = list.files(inD1,full.names=T)
df.max = read.table(Files[1],header=T,sep='\t',row.names=NULL)
spname = df.max[,1]
spnew  = c()
for (i in 1:length(spname)) {
  spnew[i] = str_remove(spname[i],"_aligned_cds")
  spnew[i] = str_remove(spnew[i],"_")
}

df.sam = read.table(Files[2],header=T,sep='\t',row.names=NULL)
df.sw  = read.table(Files[3],header=T,sep='\t',row.names=NULL)

df.max1 = df.max[,-1] %>% round(3)
df.sam1 = df.sam[,-1] %>% round(3)
df.sw1  = df.sw[,-1] %>%  round(3)

df.all = cbind(df.sw1,df.max1,df.sam1)
df.all = sapply(df.all, function(x){sprintf("%.3f",x)})



rownames(df.all) = spnew
write.csv(df.all,ouF1,row.names=T,quote=F)


# qqnorm(scale(df.sw1[,1]))
# abline(0,1,col='#661100',lwd=2)
# shapiro.test(df.sw1[,1])
# 
# qqnorm(scale(df.sw1[,2]))
# abline(0,1,col='#661100',lwd=2)
# shapiro.test(df.sw1[,2])
# 
# qqnorm(scale(df.sw1[,3]))
# abline(0,1,col='#661100',lwd=2)
# shapiro.test(df.sw1[,3])


plt = function(dat,name){
  Species      = rep(1:90,each=3)
  Phases       = rep(c("phase-0","phase-1","phase-2"),90)
  Prop         = c(t(dat))
  df           = data.frame(Species,Phases,Prop)

  mu           = ddply(df, "Phases", summarise, grp.mean=mean(Prop), grp.sd=sd(Prop))
  y            = sapply(1:3, function(x){dnorm(dat[,x],mu$grp.mean[x],mu$grp.sd[x])})
  y            = c(t(y))
  gg1 = ggplot(df,aes(x=Prop,color=Phases,fill=Phases)) + 
        geom_histogram(position="identity",alpha=0.6,bins=40) + geom_density(alpha=0.6) + 
        geom_line(aes(x=Prop,y=y), lwd=1.5) + 
        geom_vline(data=mu, aes(xintercept=grp.mean, color=Phases),linetype="dashed") +
        labs(x="Proportion", y = "Density") + theme_classic()
    
  #shaprio test
  sha = sapply(1:3, function(x) shapiro.test(dat[,x]))
  sha = c(sha[,1]$p.value,sha[,2]$p.value,sha[,3]$p.value)
  gg2 = ggplot(df,aes(sample=Prop,color=Phases)) + stat_qq() + stat_qq_line() + 
        labs(x="Theoretical", y = "Sample") + 
        annotate("text",x=-3,y=0.6,label=sprintf("shapiro.test:%.3f",sha[1]),hjust=0, col='#F8766D') + 
        annotate("text",x=-3,y=0.55,label=sprintf("shapiro.test:%.3f",sha[2]),hjust=0,col='#00BA38') + 
        annotate("text",x=-3,y=0.5,label=sprintf("shapiro.test:%.3f",sha[3]),hjust=0,col='#619CFF') 
    
  gg  = ggarrange(gg1, gg2, labels=c("A","B"), ncol=1,nrow=2)
  gg  = annotate_figure(gg,top=text_grob(name,face="bold",size=15))
  print(gg)
}

pdf(ouFig1,onefile=T,paper='a4')
par(mai=c(1,1.25,1,1.25))

plt(df.sw1,"mafft+sw")
plt(df.max1,"coati-max")
plt(df.sam1,"coati-sampling")
dev.off()


#ZnZs table of 90 species
#########################################
inD2  = "Results/ZD_sum"
ouF2  = "Results/supp4/dnds_ZnZs.csv"

F.sw  = list.files(inD2,full.names=T,pattern='sw')
F.max = list.files(inD2,full.names=T,pattern='max')
F.sam = list.files(inD2,full.names=T,pattern='sample')

#est ZnZs
ZnZs_est = function(x){
  tyN = sum(x$Zn)
  tyS = sum(x$Zs)
  tyN/(tyN+tyS)
}

#est dnds
dnds_est = function(x){
  Nd = sum(x$Nd)
  Sd = sum(x$Sd)
  N  = sum(x$N)
  S  = sum(x$S)
  
  PN = (Nd/N) 
  PS = (Sd/S) 
  
  #jukes cantor formula (1986)
  log(1-4*PN/3)/log(1-4*PS/3)
}

#generate dnds,ZnZs matrix
sum_gen = function(x){
  sum_90 = matrix(0,90,2)
  for(i in 1:length(x)){
    spec.sum     = read_tsv(x[i],col_names=T,show_col_types=F)
    sum_90[i,1] = dnds_est(spec.sum) 
    sum_90[i,2] = ZnZs_est(spec.sum) 
  }
  sum_90 = sum_90 %>% round(3)
  colnames(sum_90) = c('dN/dS','Zn/Zn+Zs')
  sum_90
}

dmat.sw  = sum_gen(F.sw)
dmat.max = sum_gen(F.max)
dmat.sam = sum_gen(F.sam)

dmat.all = cbind(dmat.sw,dmat.max,dmat.sam)
dmat.vec = sapply(dmat.all, function(x){sprintf("%.3f",x)})
dmat.All = matrix(dmat.vec,90,6)
rownames(dmat.All) = spnew
colnames(dmat.All) = colnames(dmat.all)
write.csv(dmat.All,ouF2,row.names=T,quote=F)


#omega/norm-indel-rates table
#########################################
inD3   = "../chapter3/90/Results/PISE"
ouF3   = "Results/supp4/omega_r.csv"
Files3 = list.files(inD3,full.names=T,pattern='est')

###
omega  = rep(0,90)
r012   = matrix(0,90,3)
for(i in 1:90){
  j.tab   = fromJSON(Files3[i])
  omega[i]= j.tab$omega
  g       = j.tab$gap.openning
  g012    = c(g[1]+g[4],g[2]+g[5],g[3]+g[6])
  r012[i,]= -log(1-g012)
}
norm.r    = r012/rowSums(r012)

est.all   = cbind(omega,norm.r) %>% round(3)
est.vec   = sapply(est.all, function(x){sprintf("%.3f",x)})
est.All   = matrix(est.vec,90,4)

colnames(est.All) = c('omega','r0.norm','r1.norm','r2.norm')
rownames(est.All) = spnew

write.csv(est.All,ouF3,row.names=T,quote=F)




