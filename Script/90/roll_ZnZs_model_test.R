#rolling mean of Zn/(Zn+Zs) vs dnds 
#build model based on raw data and predict on rolling stats

suppressWarnings(suppressMessages(library(tidyverse)))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(zoo))
suppressPackageStartupMessages(library(minpack.lm))

#setwd("~/Dropbox (ASU)/Indel_project/test_90_species")

#cal. rolling stats(moving average) default:k=201
ma = function(x) rollmean(x, k=201)

#extract the species name
naming = function(x) sub("\\..*", "", basename(x))


# inD1  = "Results/ZD_sum"
# inD2  = "Results/GC"
# inF   = "Results/ZnZs/neutral_ZnZs.txt"
# pat   = "max"

Files1 = list.files(inD1,full.names=T,pattern=pat)
N      = length(Files1)
namev  = unlist(lapply(Files1, function(x){naming(x)}))
##############################################test y=a*x^b+c : not bad!!

ouFig  = "Figure/ZnZs/roll_ZnZs_test.pdf"
pdf(ouFig,onefile=T)
for(i in 1:N){
  print(i)
  dat = read_tsv(Files1[i],show_col_types=F)
  dat = dat %>% mutate(Pn=Nd/N, Ps=Sd/S)
  dat = dat %>% mutate(dn=-0.75*log(1-4*Pn/3), 
                       ds=-0.75*log(1-4*Ps/3),
                       w =dn/ds)
  
  ##drop 1 row and sort by w
  dat = dat %>% filter(is.finite(w)) %>% arrange(w)
  dat = dat %>% mutate(p=Zn/(Zn+Zs), N=Zn+Zs)
  dat2= dat %>% filter(N>0)

 
  tab   = tibble(w = ma(dat$w),
                 Zs= ma(dat$Zs),
                 Zn= ma(dat$Zn),
                 p = Zn/(Zn+Zs))
  
  nls_fit = nlsLM(p ~ a*w^b+c, data=tab, 
                start=list(a=1,b=1,c=0), control=list(warnOnly=T))
  coe     = coef(nls_fit)
  ypred   = predict(nls_fit,list(tab$w))
  
  gg_b = ggplot(tab, aes(x=w, y=p)) + geom_line(col="#000000") + labs(x='dnds',y='Zn/(Zn+Zs)') + xlim(0,0.6) + ylim(0,0.6)
  gg_B = gg_b + geom_line(aes(x=w,y=!!ypred),color="#AA4499",inherit=F) 
  
  print(gg_B)
}
dev.off()



####################################### a*x^b+c test with combination of k[101,201,301]
m1 = function(x) rollmean(x, k=101)
m2 = function(x) rollmean(x, k=201)
m3 = function(x) rollmean(x, k=301)

ouFig  = "Figure/ZnZs/roll_ZnZs_test_com.pdf"

pdf(ouFig,onefile=T)
for(i in 1:N){
  print(i)
  dat = read_tsv(Files1[i],show_col_types=F)
  dat = dat %>% mutate(Pn=Nd/N, Ps=Sd/S)
  dat = dat %>% mutate(dn=-0.75*log(1-4*Pn/3), 
                       ds=-0.75*log(1-4*Ps/3),
                       w =dn/ds)
  
  ##drop 1 row and sort by w
  dat = dat %>% filter(is.finite(w)) %>% arrange(w)
  dat = dat %>% mutate(p=Zn/(Zn+Zs), N=Zn+Zs)
  dat2= dat %>% filter(N>0)
  

  # tab1 = tibble(w = m1(dat$w),Zs= m1(dat$Zs),Zn= m1(dat$Zn),p = Zn/(Zn+Zs))
  # tab2 = tibble(w = m2(dat$w),Zs= m2(dat$Zs),Zn= m2(dat$Zn),p = Zn/(Zn+Zs))
  # tab3 = tibble(w = m3(dat$w),Zs= m3(dat$Zs),Zn= m3(dat$Zn),p = Zn/(Zn+Zs))
  
  tab = tibble(w = c(m1(dat$w), m2(dat$w), m3(dat$w)),
               Zs= c(m1(dat$Zs),m2(dat$Zs),m3(dat$Zs)),
               Zn= c(m1(dat$Zn),m2(dat$Zn),m3(dat$Zn)),
               p = Zn/(Zn+Zs))
  tab = tab %>% arrange(w)
                       
  nls_fit = nlsLM(p ~ a*w^b+c, data=tab, 
                  start=list(a=1,b=1,c=0), control=list(warnOnly=T))
  coe     = coef(nls_fit)
  ypred   = predict(nls_fit,list(tab$w))
  
  gg_b = ggplot(tab, aes(x=w, y=p)) + geom_line(col="#000000") + labs(x='dnds',y='Zn/(Zn+Zs)') + xlim(0,0.6) + ylim(0,0.6)
  gg_B = gg_b + geom_line(aes(x=w,y=!!ypred),color="#AA4499",inherit=F) 
  
  print(gg_B)
}
dev.off()


####################################################test with nlsLM try, maxiter=100
ouFig = "Figure/ZnZs/roll_test2.pdf"
pdf(ouFig,onefile=T)
for(i in 1:N){
  print(i)
  dat = read_tsv(Files1[i],show_col_types=F)
  dat = dat %>% mutate(Pn=Nd/N, Ps=Sd/S)
  dat = dat %>% mutate(dn=-0.75*log(1-4*Pn/3), 
                       ds=-0.75*log(1-4*Ps/3),
                       w =dn/ds)
  
  ##drop 1 row and sort by w
  dat = dat %>% filter(is.finite(w)) %>% arrange(w)
  dat = dat %>% mutate(p=Zn/(Zn+Zs), N=Zn+Zs)
  dat2= dat %>% filter(N>0)
  
  ##calculate a logistic regression model using glm
  #dat3  = dat2 %>% mutate_at(vars(c('p')), ~ifelse(p>0,1,.))
  m      = glm(p ~ w, family=binomial(), weight=N, data=dat2)
  
  coe    = setNames(coefficients(m),NULL)    
  scale  = 1/coe[2]
  mu     = -coe[1]/coe[2]
  y0     = logis(0,mu,scale)   #y-intercept
  slope0 = d_logis(0,mu,scale) #slope at the y-intercept
  
  ##choose an asymptote that is 20% larger than our y0
  ##this seems help prevent bad fits
  ybig = if(slope0>0) y0*1.2 else y0/1.2
  
  ##estimate starting parameters for our NLS model
  f = 2*abs(ybig-y0)
  s = f/(4*slope0)
  v = ybig - (slope0>0)*f
  
  ##run nls estimation
  nls_fit = try({nlsLM(p ~ f*logis(w,0,s)+v, data=dat2, weights=N,
                  start = list(s=scale, f=f, v=v), control=nls.lm.control(maxiter=100))},silent=T) 
  
  # nls_fit = try({nlsLM(p ~ f*logis(w,0,s)+v, data=dat3, weights=N,
  #                      start = list(s=s, f=f, v=v), control=nls.lm.control(maxiter=100))},silent=T) 
  
  # nls_fit = nls(p ~ f*logis(w,0,s)+v, data=dat3, weights=N,
  #               start = list(s=scale, f=f, v=v), control=list(warnOnly=T))
  
  
  #Check convergence
  if(class(nls_fit)!="try-error"){
    if(nls_fit$convInfo$isConv) {
      #use estimated parameter estimates
      coe2 = coefficients(nls_fit)
      f    = coe2["f"]
      v    = coe2["v"]
      s    = coe2["s"]
    }
  }
 
  
  ##Calculate final parameters
  y0     = f*0.5+v
  slope0 = f/(4*s)
  ymax   =  f+v
  ymin   =  v
  
  
  tab   = tibble(w = ma(dat$w),
                 Zs= ma(dat$Zs),
                 Zn= ma(dat$Zn),
                 p = Zn/(Zn+Zs))
  
  #ypred = predict(m,list(w=tab$w),type='response')
  #ypred = f*logis(tab$w,u,s)+v
  ypred = f*logis(tab$w,0,s)+v
  
  ##rmse
  rmse = sqrt(mean((tab$p-ypred)^2))
  
  tablong = tab %>% pivot_longer(cols=Zs:Zn, names_to="Stat", values_to="Count")
  ##create plots
  gg_a = ggplot(tablong, aes(x=w, y=Count, color=Stat)) + geom_line() + theme(legend.position=c(.05,.95)) + labs(x='dnds') + xlim(0,0.6) 
  
  gg_b = ggplot(tab, aes(x=w, y=p)) + geom_line(col="#000000") + labs(x='dnds',y='Zn/(Zn+Zs)') + xlim(0,0.6) + ylim(0,0.6)
  gg_b = gg_b + geom_hline(yintercept=neuZ[i], linetype="dashed", color="#661100") + 
    annotate("text",x=0,y=0.6,label=sprintf("GC:%.3f",gc.perc[i]),hjust=0) + 
    annotate("text",x=0,y=0.57,label=sprintf("sample size:%i",nrow(dat)),hjust=0) +
    annotate("text",x=0,y=0.54,label=sprintf("neutral:%.3f",neuZ[i]),hjust=0,color="#661100") 
  gg_B  = gg_b + geom_line(aes(x=w,y=!!ypred),color="#AA4499",inherit=F) 
  
  # gg_B1 = gg_B + geom_hline(yintercept=ymax, linetype="dashed", color="#56B4E9") + 
  #         annotate("text",x=0,y=0.51,label=sprintf("asymptote:%.3f",ymax),hjust=0,color="#56B4E9") + 
  #         annotate("text",x=0,y=0.48,label=sprintf("init-slope:%.3f",slope0),hjust=0) 
  
  gg_B1 = gg_B + {if(slope0>0) geom_hline(yintercept=ymax, linetype="dashed", color="#56B4E9")  else 
    geom_hline(yintercept=ymin, linetype="dashed", color="#56B4E9")} + 
    {if(slope0>0) annotate("text",x=0,y=0.51,label=sprintf("asymptote:%.3f",ymax),hjust=0,color="#56B4E9") else
      annotate("text",x=0,y=0.51,label=sprintf("asymptote:%.3f",ymin),hjust=0,color="#56B4E9")}
  
  gg_B1 = gg_B1 + annotate("text",x=0,y=0.48,label=sprintf("init-slope:%.3f",slope0),hjust=0) 
  
  
  gg    = ggarrange(gg_a, gg_B1, labels=c("A","B"), ncol=1,nrow=2)
  gg    = annotate_figure(gg,top=text_grob(namev[i],face="bold",size=10))
  print(gg)
  
}
dev.off()

#####################################################test with combination of k (with previous nlsLM setting)
m1 = function(x) rollmean(x, k=101)
m2 = function(x) rollmean(x, k=201)
m3 = function(x) rollmean(x, k=301)

ouFig = "Figure/ZnZs/roll_test3.pdf"
pdf(ouFig,onefile=T)
for(i in 1:N){
  print(i)
  dat = read_tsv(Files1[i],show_col_types=F)
  dat = dat %>% mutate(Pn=Nd/N, Ps=Sd/S)
  dat = dat %>% mutate(dn=-0.75*log(1-4*Pn/3), 
                       ds=-0.75*log(1-4*Ps/3),
                       w =dn/ds)
  
  ##drop 1 row and sort by w
  dat = dat %>% filter(is.finite(w)) %>% arrange(w)
  dat = dat %>% mutate(p=Zn/(Zn+Zs), N=Zn+Zs)
  dat2= dat %>% filter(N>0)
  
  ##calculate a logistic regression model using glm
  m      = glm(p ~ w, family=binomial(), weight=N, data=dat2)
  
  coe    = setNames(coefficients(m),NULL)    
  scale  = 1/coe[2]
  mu     = -coe[1]/coe[2]
  y0     = logis(0,mu,scale)   #y-intercept
  slope0 = d_logis(0,mu,scale) #slope at the y-intercept
  
  ##choose an asymptote that is 20% larger than our y0
  ##this seems help prevent bad fits
  ybig = if(slope0>0) y0*1.2 else y0/1.2
  
  ##estimate starting parameters for our NLS model
  f = 2*abs(ybig-y0)
  s = f/(4*slope0)
  v = ybig - (slope0>0)*f
  
  ##run nls estimation
  nls_fit = try({nlsLM(p ~ f*logis(w,0,s)+v, data=dat2, weights=N,
                       start = list(s=s, f=f, v=v), control=nls.lm.control(maxiter=100))},silent=T) 
  
  
  
  #Check convergence
  if(class(nls_fit)!="try-error"){
    if(nls_fit$convInfo$isConv) {
      #use estimated parameter estimates
      coe2 = coefficients(nls_fit)
      f    = coe2["f"]
      v    = coe2["v"]
      s    = coe2["s"]
    }
  }
  
  
  ##Calculate final parameters
  y0     = f*0.5+v
  slope0 = f/(4*s)
  ymax   =  f+v
  ymin   =  v
  
  
  tab = tibble(w = c(m1(dat$w), m2(dat$w), m3(dat$w)),
               Zs= c(m1(dat$Zs),m2(dat$Zs),m3(dat$Zs)),
               Zn= c(m1(dat$Zn),m2(dat$Zn),m3(dat$Zn)),
               p = Zn/(Zn+Zs))
  tab = tab %>% arrange(w)
  
  #ypred = predict(m,list(w=tab$w),type='response')
  #ypred = f*logis(tab$w,u,s)+v
  ypred = f*logis(tab$w,0,s)+v
  
  ##rmse
  rmse = sqrt(mean((tab$p-ypred)^2))
  
  tablong = tab %>% pivot_longer(cols=Zs:Zn, names_to="Stat", values_to="Count")
  ##create plots
  gg_a = ggplot(tablong, aes(x=w, y=Count, color=Stat)) + geom_line() + theme(legend.position=c(.05,.95)) + labs(x='dnds') + xlim(0,0.6) 
  
  gg_b = ggplot(tab, aes(x=w, y=p)) + geom_line(col="#000000") + labs(x='dnds',y='Zn/(Zn+Zs)') + xlim(0,0.6) + ylim(0,0.6)
  gg_b = gg_b + geom_hline(yintercept=neuZ[i], linetype="dashed", color="#661100") + 
    annotate("text",x=0,y=0.6,label=sprintf("GC:%.3f",gc.perc[i]),hjust=0) + 
    annotate("text",x=0,y=0.57,label=sprintf("sample size:%i",nrow(dat)),hjust=0) +
    annotate("text",x=0,y=0.54,label=sprintf("neutral:%.3f",neuZ[i]),hjust=0,color="#661100") 
  gg_B  = gg_b + geom_line(aes(x=w,y=!!ypred),color="#AA4499",inherit=F) 
  
  # gg_B1 = gg_B + geom_hline(yintercept=ymax, linetype="dashed", color="#56B4E9") + 
  #         annotate("text",x=0,y=0.51,label=sprintf("asymptote:%.3f",ymax),hjust=0,color="#56B4E9") + 
  #         annotate("text",x=0,y=0.48,label=sprintf("init-slope:%.3f",slope0),hjust=0) 
  
  gg_B1 = gg_B + {if(slope0>0) geom_hline(yintercept=ymax, linetype="dashed", color="#56B4E9")  else 
    geom_hline(yintercept=ymin, linetype="dashed", color="#56B4E9")} + 
    {if(slope0>0) annotate("text",x=0,y=0.51,label=sprintf("asymptote:%.3f",ymax),hjust=0,color="#56B4E9") else
      annotate("text",x=0,y=0.51,label=sprintf("asymptote:%.3f",ymin),hjust=0,color="#56B4E9")}
  
  gg_B1 = gg_B1 + annotate("text",x=0,y=0.48,label=sprintf("init-slope:%.3f",slope0),hjust=0) 
  
  
  gg    = ggarrange(gg_a, gg_B1, labels=c("A","B"), ncol=1,nrow=2)
  gg    = annotate_figure(gg,top=text_grob(namev[i],face="bold",size=10))
  print(gg)
  
  coefv[i,] = c(y0,slope0,ymax,nrow(dat),rmse)
  
}
dev.off()

####################################################test with SSlogis/getInitial [failed]
# logis2 = function(x,location,scale,alpha) {
#     xx = (x-location)/scale
#     alpha/(1+exp(-xx))
# }

ouFig = "Figure/ZnZs/roll_test4.pdf"
pdf(ouFig,onefile=T)
for(i in 1:N){
  print(i)
  dat = read_tsv(Files1[i],show_col_types=F)
  dat = dat %>% mutate(Pn=Nd/N, Ps=Sd/S)
  dat = dat %>% mutate(dn=-0.75*log(1-4*Pn/3), 
                       ds=-0.75*log(1-4*Ps/3),
                       w =dn/ds)
  
  ##drop 1 row and sort by w
  dat = dat %>% filter(is.finite(w)) %>% arrange(w)
  dat = dat %>% mutate(p=Zn/(Zn+Zs), N=Zn+Zs)
  dat2= dat %>% filter(N>0)
  
  ##calculate a logistic regression model using glm
  m      = glm(p ~ w, family=binomial(), weight=N, data=dat2)
  
  coe    = setNames(coefficients(m),NULL)    
  scale  = 1/coe[2]
  mu     = -coe[1]/coe[2]
  y0     = logis(0,mu,scale)   #y-intercept
  slope0 = d_logis(0,mu,scale) #slope at the y-intercept
  
  ##choose an asymptote that is 20% larger than our y0
  ##this seems help prevent bad fits
  ybig = if(slope0>0) y0*1.2 else y0/1.2
  
  ##estimate starting parameters for our NLS model
  # f = 2*abs(ybig-y0)
  # s = f/(4*slope0)
  # v = ybig - (slope0>0)*f
  
  ##test SSlogis 
  x_logis = function(location,scale,y) location + scale*log(y/(1-y))
  xmid    = x_logis(mu,scale,ybig/2)
  
  #SSlogis(dat3$w,ybig,xmid,scale)
  #getInitial(p ~ f*1/(1+exp(-(x-0)/s))+v,data=dat2)
  ss = getInitial(p~ SSlogis(w,ybig,xmid,scale),data=data.frame(p=dat2$p,w=dat2$w))  
  
  s  = unname(ss['scale'])
  f  = unname(2*abs(ss['ybig']-y0))
  v  = unname(ss['ybig']-(slope0>0)*f)
  
  ##run nls estimation
  nls_fit = try({nlsLM(p ~ f*logis(w,0,s)+v, data=dat2, weights=N, na.action=na.omit,
                       start = list(s=s, f=f, v=v), control=nls.lm.control(maxiter=100))},silent=T) 

    
  #Check convergence
  if(class(nls_fit)!="try-error"){
    if(nls_fit$convInfo$isConv) {
      #use estimated parameter estimates
      coe2 = coefficients(nls_fit)
      f    = coe2["f"]
      v    = coe2["v"]
      s    = coe2["s"]
    }
  }
  
  
  ##Calculate final parameters
  y0     = f*0.5+v
  slope0 = f*s*0.25
  ymax   =  f+v
  ymin   =  v
  
  
  tab   = tibble(w = ma(dat$w),
                 Zs= ma(dat$Zs),
                 Zn= ma(dat$Zn),
                 p = Zn/(Zn+Zs))
  
  #ypred = predict(m,list(w=tab$w),type='response')
  #ypred = f*logis(tab$w,u,s)+v
  ypred = f*logis(tab$w,0,s)+v
  
  ##rmse
  rmse = sqrt(mean((tab$p-ypred)^2))
  
  tablong = tab %>% pivot_longer(cols=Zs:Zn, names_to="Stat", values_to="Count")
  ##create plots
  gg_a = ggplot(tablong, aes(x=w, y=Count, color=Stat)) + geom_line() + theme(legend.position=c(.05,.95)) + labs(x='dnds') + xlim(0,0.6) 
  
  gg_b = ggplot(tab, aes(x=w, y=p)) + geom_line(col="#000000") + labs(x='dnds',y='Zn/(Zn+Zs)') + xlim(0,0.6) + ylim(0,0.6)
  gg_b = gg_b + geom_hline(yintercept=neuZ[i], linetype="dashed", color="#661100") + 
    annotate("text",x=0,y=0.6,label=sprintf("GC:%.3f",gc.perc[i]),hjust=0) + 
    annotate("text",x=0,y=0.57,label=sprintf("sample size:%i",nrow(dat)),hjust=0) +
    annotate("text",x=0,y=0.54,label=sprintf("neutral:%.3f",neuZ[i]),hjust=0,color="#661100") 
  gg_B  = gg_b + geom_line(aes(x=w,y=!!ypred),color="#AA4499",inherit=F) 
  
  # gg_B1 = gg_B + geom_hline(yintercept=ymax, linetype="dashed", color="#56B4E9") + 
  #         annotate("text",x=0,y=0.51,label=sprintf("asymptote:%.3f",ymax),hjust=0,color="#56B4E9") + 
  #         annotate("text",x=0,y=0.48,label=sprintf("init-slope:%.3f",slope0),hjust=0) 
  
  gg_B1 = gg_B + {if(slope0>0) geom_hline(yintercept=ymax, linetype="dashed", color="#56B4E9")  else 
    geom_hline(yintercept=ymin, linetype="dashed", color="#56B4E9")} + 
    {if(slope0>0) annotate("text",x=0,y=0.51,label=sprintf("asymptote:%.3f",ymax),hjust=0,color="#56B4E9") else
      annotate("text",x=0,y=0.51,label=sprintf("asymptote:%.3f",ymin),hjust=0,color="#56B4E9")}
  
  gg_B1 = gg_B1 + annotate("text",x=0,y=0.48,label=sprintf("init-slope:%.3f",slope0),hjust=0) 
  
  
  gg    = ggarrange(gg_a, gg_B1, labels=c("A","B"), ncol=1,nrow=2)
  gg    = annotate_figure(gg,top=text_grob(namev[i],face="bold",size=10))
  print(gg)
  
  coefv[i,] = c(y0,slope0,ymax,nrow(dat),rmse)
  
}
dev.off()

######################################################test with logit-ypread if nls failed. {looks good}
ouFig = "Figure/ZnZs/roll_test5.pdf"
pdf(ouFig,onefile=T)
for(i in 1:N){
  print(i)
  dat = read_tsv(Files1[i],show_col_types=F)
  dat = dat %>% mutate(Pn =Nd/N, Ps=Sd/S)
  dat = dat %>% mutate(dn =-0.75*log(1-4*Pn/3), 
                       ds =-0.75*log(1-4*Ps/3),
                       w  =dn/ds, 
                       len=lenA+lenB)
  
  ##drop 1 row and sort by w
  dat = dat %>% filter(is.finite(w)) %>% arrange(w)
  dat = dat %>% mutate(p=Zn/(Zs+Zn), N=Zn+Zs)           
  dat2= dat %>% filter(N>0)
  
  ##calculate a logistic regression model using glm. p <- [0,1]
  m   = glm(p ~ w, family=binomial(), weights=N, data=dat2)
  
  coe    = setNames(coefficients(m),NULL)    
  scale  = 1/coe[2]
  mu     = -coe[1]/coe[2]
  y0     = logis(0,mu,scale)   #y-intercept
  slope0 = d_logis(0,mu,scale) #slope at the y-intercept
  

  ##choose an asymptote that is 20% larger than our y0
  ##this seems help prevent bad fits
  ybig = if(slope0>0) y0*1.2 else y0/1.2
  
  ##estimate starting parameters for our NLS model
  f = 2*abs(ybig-y0)
  s = f/(4*slope0)
  v = ybig - (slope0>0)*f
  
  ##run nls estimation
  nls_fit = try({nlsLM(p ~ f*logis(w,0,s)+v, data=dat2, weights=N,
                       start = list(s=scale, f=f, v=v), control=nls.lm.control(maxiter=100))},silent=T) 
  
  
  #Check convergence
  if(class(nls_fit)!="try-error"){
    if(nls_fit$convInfo$isConv) {
      #use estimated parameter estimates
      coe2 = coefficients(nls_fit)
      f    = coe2["f"]
      v    = coe2["v"]
      s    = coe2["s"]
    }
  }
  
  
  ##Calculate final parameters
  y0     = f*0.5+v
  slope0 = f/(4*s)
  ymax   =  f+v
  ymin   =  v
  
  
  tab   = tibble(w = ma(dat$w),
                 Zs= ma(dat$Zs),
                 Zn= ma(dat$Zn),
                 p = Zn/(Zn+Zs))
  
  ypred = f*logis(tab$w,0,s)+v
  
  if(class(nls_fit)=="try-error"){#try-logit
    #plot(logis(w,mu,scale)~w,data=tab)
    ypred = logis(tab$w,mu,scale)
  }
  
  
  ##rmse
  rmse = sqrt(mean((tab$p-ypred)^2))
  
  tablong = tab %>% pivot_longer(cols=Zs:Zn, names_to="Stat", values_to="Count")
  ##create plots
  gg_a = ggplot(tablong, aes(x=w, y=Count, color=Stat)) + geom_line() + theme(legend.position=c(.05,.95)) + labs(x='dnds') + xlim(0,0.6) 
  
  gg_b = ggplot(tab, aes(x=w, y=p)) + geom_line(col="#000000") + labs(x='dnds',y='Zn/(Zn+Zs)') + xlim(0,0.6) + ylim(0,0.6)
  gg_b = gg_b + geom_hline(yintercept=neuZ[i], linetype="dashed", color="#661100") + 
    annotate("text",x=0,y=0.6,label=sprintf("GC:%.3f",gc.perc[i]),hjust=0) + 
    annotate("text",x=0,y=0.57,label=sprintf("sample size:%i",nrow(dat)),hjust=0) +
    annotate("text",x=0,y=0.54,label=sprintf("neutral:%.3f",neuZ[i]),hjust=0,color="#661100") 
  gg_B  = gg_b + geom_line(aes(x=w,y=!!ypred),color="#AA4499",inherit=F) 
  
  # gg_B1 = gg_B + geom_hline(yintercept=ymax, linetype="dashed", color="#56B4E9") + 
  #         annotate("text",x=0,y=0.51,label=sprintf("asymptote:%.3f",ymax),hjust=0,color="#56B4E9") + 
  #         annotate("text",x=0,y=0.48,label=sprintf("init-slope:%.3f",slope0),hjust=0) 
  
  gg_B1 = gg_B + {if(slope0>0) geom_hline(yintercept=ymax, linetype="dashed", color="#56B4E9")  else 
    geom_hline(yintercept=ymin, linetype="dashed", color="#56B4E9")} + 
    {if(slope0>0) annotate("text",x=0,y=0.51,label=sprintf("asymptote:%.3f",ymax),hjust=0,color="#56B4E9") else
      annotate("text",x=0,y=0.51,label=sprintf("asymptote:%.3f",ymin),hjust=0,color="#56B4E9")}
  
  gg_B1 = gg_B1 + annotate("text",x=0,y=0.48,label=sprintf("init-slope:%.3f",slope0),hjust=0) 
  
  
  gg    = ggarrange(gg_a, gg_B1, labels=c("A","B"), ncol=1,nrow=2)
  gg    = annotate_figure(gg,top=text_grob(namev[i],face="bold",size=10))
  print(gg)
  
}
dev.off()

#######################################################
