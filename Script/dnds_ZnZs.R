#rolling mean of typeN,typeS vs dnds 
#build model based on raw data and predict on rolling stats

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(zoo))

#setwd("~/Dropbox (ASU)/Indel_project/test_human_mouse_rat")

# inF   = "Results/dnds.IDtype.sum.tsv"
# ouFig = "Figure/znzs/dnds_ZnZs.pdf"
main = function(inF,ouFig){
  dat = read_tsv(inF,show_col_types=F)
  
  #calculate Stats
  dat = dat %>% mutate(Pn = Nd/N, Ps=Sd/S)
  dat = dat %>% mutate(dn = -0.75*log(1-4*Pn/3), 
                       ds = -0.75*log(1-4*Ps/3),
                       w  = dn/ds)
  
  #drop 1 row and sort by w
  dat = dat %>% filter(is.finite(w)) %>% arrange(w)
  
  #Calculating rolling stats(moving average)   >>>>>>>>>>>>>default:k=201
  ma  = function(x) rollmean(x, k=201)
  tab = tibble(w = ma(dat$w),
               Zs= ma(dat$typeS),
               Zn= ma(dat$typeN))
  tablong = tab %>% pivot_longer(cols=Zs:Zn, names_to="Stat", values_to="Count")
  
  #create plots
  gg_a = ggplot(tablong, aes(x=w, y=Count, color=Stat)) + geom_line() + theme(legend.position=c(.05,.95))
  gg_b = ggplot(tab, aes(x=w, y=Zn/(Zn+Zs))) + geom_line(col="#56B4E9")
  
  gg   = ggarrange(gg_a, gg_b,
                   labels=c("A","B"), ncol=1,nrow=2)
  print(gg)
  
  #binomial/logistic regression
  tab2  = dat %>% filter(w<1/3)                                 #raw data
  m     = glm(cbind(typeN,typeS)~w,data=tab2,family=binomial())
  ypred = predict(m,list(w=tab$w), type="response") 
 
  gg_B = gg_b + geom_line(aes(x=w,y=!!ypred),color="#AA4499",inherit=F) + geom_hline(yintercept=0.2, linetype="dashed", 
                                                                                     color = "#009E73", size=0.5)
  print(gg_B)
  
  
  #>>>>>>>>>>>>output design
  gg2   = ggarrange(gg_a, gg_B,
                   labels=c("A","B"), ncol=1,nrow=2)
  print(gg2)
  
  
  #>>>play w, (typeN,typeS plateau as w gets larger). 
  for(i in 1:10){
    tab2 = dat %>% filter(w<1/i)
    m    = glm(cbind(typeN,typeS)~w,data=tab2,family=binomial())
    ypred= predict(m,list(w=tab$w), type="response") 
    print(gg_b + geom_line(aes(x=w,y=!!ypred),color="#AA4499"))
    #print(gg_b + geom_smooth(aes(y),color="#AA4499"))
    Sys.sleep(1)
  }
  dev.off()
  
  
  #poisson distribution 
  lenA = ma(dat$lenA)
  lenB = ma(dat$lenB)
  tab3 = cbind(tab,"lenA"=lenA,"lenB"=lenB)
  
  m1     = glm(typeN~w+offset(log((lenA+lenB)/2)),data=tab2,family=poisson)
  ypred1 = predict(m1,data.frame(w=tab3$w,lenA=tab3$lenA,lenB=tab3$lenB),type="response")
  
  m2     = glm(typeS~w+offset(log((lenA+lenB)/2)),data=tab2,family=poisson)
  ypred2 = predict(m2,list(w=tab3$w,lenA=tab3$lenA,lenB=tab3$lenB),type="response") 
  
  ##compare the slope
  m1$coefficients
  m2$coefficients
  
  gg_a1 = ggplot(tab, aes(x=w, y=Zn)) + geom_line(col="#D55E00") 
  gg_a2 = ggplot(tab, aes(x=w, y=Zs)) + geom_line(col="#0072B2") 
  
  print(gg_a1)
  print(gg_a2)
  
  gg_a1 + geom_line(aes(x=w,y=!!ypred1),color="#AA4499")
  gg_a2 + geom_line(aes(x=w,y=!!ypred2),color="#AA4499")
  
}

gg_A = ggplot(tab,aes(x=w)) + geom_line(aes(y=Zn,color="Zn")) +
                              geom_line(aes(y=Zs,color="Zs")) +
                              geom_line(aes(y=!!ypred1,color="Zn-pred")) +
                              geom_line(aes(y=!!ypred2,color="Zs-pred")) + 
                              scale_color_manual(name="Z series",values=c("Zn" = "#D55E00","Zs" = "#D55E00","Zn-pred" = "#0072B2","Zs-pred" = "#0072B2" ))


#if he doesn't like the w<1 causing the lines rough, swithch to geom_smooth.                        


gg1   = ggarrange(gg_A, gg_B,
                 labels=c("A","B"), ncol=1,nrow=2)
print(gg1)


##############################################
args = commandArgs(trailingOnly=TRUE)
main(args[1], args[2])
