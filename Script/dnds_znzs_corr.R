#rolling mean of znzs vs dnds 
#build model based on raw data and predict on rolling stats

suppressWarnings(suppressMessages(library(tidyverse)))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(zoo))
suppressPackageStartupMessages(library(minpack.lm))

#setwd("~/Dropbox (ASU)/Indel_project/test_human_mouse_rat")

ma = function(x) rollmean(x, k=201)

logis = function(x,location,scale) {
  xx = (x-location)/scale
  1/(1+exp(-xx))
}

d_logis = function(x, location, scale) {
  xx = (x-location)/scale
  (1/scale)*exp(-xx)/((1+exp(-xx))^2)
}

######################################
# inF   = "Results/dnds_znzs_sum.tsv"
# ouFig = "Results/dnds_znzs.pdf"

main = function(inF,ouFig){
  dat = read_tsv(inF,show_col_types=F)
  dat = dat %>% mutate(Pn=Nd/N, Ps=Sd/S)
  dat = dat %>% mutate(dn=-0.75*log(1-4*Pn/3), 
                       ds=-0.75*log(1-4*Ps/3),
                       w =dn/ds)
  
  dat = dat %>% filter(is.finite(w)) %>% arrange(w)
  
  
  tab = tibble(w  = ma(dat$w),
               Zs = ma(dat$Zs),
               Zn = ma(dat$Zn),
               ZS = ma(dat$ZS),
               ZN = ma(dat$ZN),
               zn = Zn/ZN,
               zs = Zs/ZS,
               p  = zn/(zn+zs),
               N  = ZN+ZS)
  
  
  ##fit [weights not working]
  m      = glm(p ~ w, family=binomial(), weights=N, data=tab)
  
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
  
  nls_fit = nlsLM(p ~ f*logis(w,0,s)+v, data=tab, weights=N,
                  start = list(s=scale, f=f, v=v), control=nls.lm.control(maxiter=100))
  
  
  ##final estimates
  coe2 = coefficients(nls_fit)
  f    = coe2["f"]
  v    = coe2["v"]
  s    = coe2["s"]
  
  ypred = predict(nls_fit,list(tab$w))
  
  pdf(ouFig,onefile=T)
  
  tablong = tab %>% pivot_longer(cols=zs:zn, names_to="Stat", values_to="Count")
  gg_a    = ggplot(tablong, aes(x=w, y=Count, color=Stat)) + geom_line() + theme(legend.position=c(.05,.95)) + labs(x='dnds') + xlim(0,0.6) 
  
  gg_b = ggplot(tab, aes(x=w, y=p)) + geom_line(col="#000000") + labs(x='dnds',y='zn/(zn+zs)') + xlim(0,0.6) 
  gg_b = gg_b + geom_line(aes(x=w,y=!!ypred),color="#AA4499",inherit=F) 
  
  gg   = ggarrange(gg_a, gg_b, labels=c("A","B"), ncol=1,nrow=2)
  print(gg)
  dev.off()
}


###########################################
args = commandArgs(trailingOnly=T)
main(args[1], args[2])

