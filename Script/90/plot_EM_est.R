library(ggplot2)

setwd("~/Dropbox (ASU)/Indel_project/Script/90")
getwd()




# hh = read.table(ouF1,header = FALSE,sep ="\t",quote = "")
# hh = data.matrix(hh)
# qqplot(hh[1,],hh[1,])
pdf(ouFig)   
##Draw LL   


##Draw tau

##Draw omega

##Draw simgas
par(mfrow=c(3,2))


#p=plot(x,y)
#q=plot(m,n)
#multiplot(p)
#multiplot(q)
dev.off()