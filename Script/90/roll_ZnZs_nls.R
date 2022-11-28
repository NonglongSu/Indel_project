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

#cal the GC % 
gc_perc = function(x){
  dat    = read_tsv(x,show_col_types=F)
  datsum = unname(colSums(dat))
  gc.perc= (datsum[1]+datsum[2])/(datsum[3]+datsum[4])
  return(gc.perc)
}

#helper functions for logistic curve 
logis = function(x,location,scale) {
  xx = (x-location)/scale
  1/(1+exp(-xx))
}

#helper functions for logistic  derivative
d_logis = function(x, location, scale) {
  xx = (x-location)/scale
  (1/scale)*exp(-xx)/((1+exp(-xx))^2)
}


##################################################################
# ouFig = "../test_human_mouse_rat/Figure/znzs/dnds_ZnZs.pdf"
# 
# inD1   = "Results/ZD_sum"
# inD2   = "Results/GC"
# inF   = "Results/ZnZs/neutral_ZnZs.txt"
# inF1  = "MLE90.tab.csv"
# ouFig = "Figure/ZnZs/roll_ZnZs_max.pdf"
# ouF   = "Results/ZnZs/coef_max.tsv"
# pat   = "max"

main = function(inD1,inD2,inF,inF1,ouFig,ouF,pat){
  Files1 = list.files(inD1,full.names=T,pattern=pat)
  Files2 = list.files(inD2,full.names=T)
  
  N     = length(Files1)
  namev = unlist(lapply(Files1, function(x){naming(x)}))
  spnew  = c()
  for (i in 1:length(namev)) {
    spnew[i] = str_remove(namev[i],"_aligned_cds")
    spnew[i] = str_remove(spnew[i],"_")
  }
  
  #neutral Zn/(Zn+Zs)
  neuZ = read.table(inF,header=T) 
  neuZ = neuZ$neutral_cod_ratio
  
  #cal. the GC%
  gc.perc = unlist(lapply(Files2, function(x){gc_perc(x)}))
  
  #cal. the tau
  dat = read.table(inF1,header=T,sep=",",quote=NULL)
  tau = dat$tau
  
  
  coefv = matrix(0,N,6)
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
    
    
    #############################################################
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
    nls_fit = nls(p ~ f*logis(w,0,s)+v, data=dat2, weights=N,
                  start = list(s=s, f=f, v=v), control=list(warnOnly=T))
    
    
    #Check convergence
    if(nls_fit$convInfo$isConv) {
      #use estimated parameter estimates
      coe2 = coefficients(nls_fit)
      f    = coe2["f"]
      v    = coe2["v"]
      s    = coe2["s"]
    }
    
    ##Calculate final parameters
    y0     = f*0.5+v
    slope0 = f/(4*s)
    ymax   =  f+v
    ymin   =  v
    
    
    #################################
    tab   = tibble(w = ma(dat$w),
                   Zs= ma(dat$Zs),
                   Zn= ma(dat$Zn),
                   p = Zn/(Zn+Zs))
    
    #ypred = predict(m,list(w=tab$w),type='response')
    ypred = f*logis(tab$w,0,s)+v
    
    ##rmse
    rmse = sqrt(mean((tab$p-ypred)^2))
    
    ##correlation
    corr = cor(ypred,tab$p)
    
    ##gg-plot
    tablong = tab %>% pivot_longer(cols=Zs:Zn, names_to="Stat", values_to="Count")
    gg_a    = ggplot(tablong, aes(x=w, y=Count, color=Stat)) + geom_line() + theme(legend.position=c(.05,.95)) + labs(x='dnds') + xlim(0,0.6) 
    
    gg_b    = ggplot(tab, aes(x=w, y=p)) + geom_line(col="#000000") + labs(x='dnds',y='Zn/(Zn+Zs)') + xlim(0,0.6) + ylim(0,0.6)
    gg_b    = gg_b + geom_hline(yintercept=neuZ[i], linetype="dashed", color="#661100") + 
              annotate("text",x=0,y=0.6,label=sprintf("GC:%.3f",gc.perc[i]),hjust=0) + 
              annotate("text",x=0,y=0.57,label=sprintf("sample size:%i",nrow(dat)),hjust=0) +
              annotate("text",x=0,y=0.54,label=sprintf("neutral:%.3f",neuZ[i]),hjust=0,color="#661100") 
    gg_B    = gg_b + geom_line(aes(x=w,y=!!ypred),color="#AA4499",inherit=F) 
  
    gg_B1   = gg_B + {if(slope0>0) geom_hline(yintercept=ymax, linetype="dashed", color="#56B4E9")  else 
              geom_hline(yintercept=ymin, linetype="dashed", color="#56B4E9")} + 
              {if(slope0>0) annotate("text",x=0,y=0.51,label=sprintf("asymptote:%.3f",ymax),hjust=0,color="#56B4E9") else
              annotate("text",x=0,y=0.51,label=sprintf("asymptote:%.3f",ymin),hjust=0,color="#56B4E9")}
            
    gg_B1   = gg_B1 + annotate("text",x=0,y=0.48,label=sprintf("init slope:%.3f",slope0),hjust=0) +
                      annotate("text",x=0,y=0.45,label=sprintf("tau:%.3f",tau[i]),hjust=0) + 
                      annotate("text",x=0,y=0.42,label=sprintf("r:%.3f",corr),hjust=0) 
      
    gg      = ggarrange(gg_a, gg_B1, labels=c("A","B"), ncol=1,nrow=2)
    gg      = annotate_figure(gg,top=text_grob(namev[i],face="bold",size=10))
    print(gg)
    
    coefv[i,] = c(nrow(dat),rmse,corr,y0,slope0,ymax)
  }
  dev.off()
  
  ##output a table
  coe.tab = data.frame(species=spnew,neutral=neuZ,GC=gc.perc, 
            yintercept=coefv[,4],init_slope=coefv[,5],ymax=coefv[,6],sample_size=coefv[,1],rmse=coefv[,2],corr=coefv[,3])
  write.table(coe.tab, ouF, sep="\t",append=F, quote=F, row.names=F, col.names=T)
}


####################################################
args = commandArgs(trailingOnly=T)
main(args[1],args[2],args[3],args[4],args[5],args[6],args[7])
