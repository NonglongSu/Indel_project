#Data analysis of simulation(100) and real data(90)
#setwd("~/Dropbox (ASU)/Indel_project/Script/90")

suppressPackageStartupMessages(library(stats))
suppressPackageStartupMessages(library(jsonlite))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggplot2))



##>>>est90
inD = "../../test_90_species/Results/est90"

main = function(inD){
  Files = list.files(inD, full.names=TRUE)
  
  dmat = matrix(0,length(Files),8)
  for(i in 1:length(Files)){
    Jtmp     = fromJSON(Files[i])
    dmat[i,] = c(Jtmp$par, Jtmp$p.intermed[nrow(Jtmp$p.intermed),10])
  }
  ssv  = dmat[,1:6]/dmat[,8]
  dMat = cbind(ssv,dmat[,7], dmat[,8]) 
  colnames(dMat)=c('s1','s2','s3','s4','s5','s6','omega','tau')
  
  #>boxplot
  boxplot(dMat, ylab='value')
  
  #>sorted boxplot
  dF = as.data.frame(as.table(dMat))
  colnames(dF) = c('haha','paras','value')
  ggplot(dF, aes(x=paras, y=value)) + geom_boxplot()
  dF_sorted = dF %>% mutate(paras = fct_reorder(paras,-value))
  ggplot(dF_sorted, aes(x=paras, y=value)) + geom_boxplot()
  
  #>flip he boxplot
  ggplot(dF_sorted, aes(x=paras,y=value,color=paras)) + geom_boxplot() + coord_flip() 
 
  #>violin plot
  ggplot(dF_sorted, aes(x=paras,y=value,color=paras)) + geom_violin() + coord_flip()
  
  #>lineplot
  ggplot(dF_sorted, aes(x=paras,y=value,color=paras)) + geom_line(size=1) + coord_flip()
   
  #>Dot strip plot
  ggplot(dF_sorted, aes(x=paras,y=value,color=paras)) + geom_point(size=2, alpha=0.6) + coord_flip()
  
  
  #create basic g
  g = ggplot(dF_sorted, aes(x=paras,y=value,color=paras)) + coord_flip()
  
  #>boxplot(median)+dotplot
  ggplot(dF_sorted, aes(x=paras,y=value,color=paras)) + geom_boxplot() + coord_flip()
  g + geom_boxplot(outlier.alpha=0) + geom_point(size=2, alpha=0.6) 
  
  #>hightlight the mean
  set.seed(2021)
  g + geom_jitter(size=2,alpha=0.25,width=0.2) + stat_summary(fun=mean,geom='point',size=5)
  
  #>>>try w vs tau
  plot(dMat[,7]~dMat[,8], xlab='tau',ylab='omega')
  y = log(dMat[,7])
  x = dMat[,8]
  plot(y~x, xlab='tau',ylab='omega')
  abline(lm(y~x),col='red',lwd=2)
  plot(lowess(x,y),col='blue')
  
  cor(x,y)
  
}

dMat[order(dMat[,8]),]





#######################################################
#>>>global
#4*4 GTR matrix
GTR = function(Si, pai){
  r1 = matrix(0, 4, 4)
  r1[lower.tri(r1)] = Si
  r  = r1 + t(r1)
  r  = t(r*pai)
  diag(r) = -rowSums(r)
  return(r)
}

#Construct a 61*61 MG94 matrix
MG94 = function(g, oga, cod){
  R  = matrix(0, 61, 61)
  for (i in 1:61) {
    for (j in 1:61) {
      if(i == j){
        R[i, j] = 0
      }else if(sum(cod[i, ] != cod[j, ]) > 1){#more than 1 nucleotide changes
        R[i, j] = 0
      }else{
        if(codonstrs[j] %in% syn[[codonstrs[i]]]){#syn_subs
          w = 1
        }else{#nonsyn_subs
          w = oga
        }
        pos = which(cod[i, ] != cod[j, ])
        x   = which(DNA_BASES == cod[i, pos])
        y   = which(DNA_BASES == cod[j, pos])
        R[i, j] = w * g[x, y]
      }
    }
  }
  diag(R) = -rowSums(R)
  
  return(R)
}

#Log-likelihood
LL = function(theta){
  rmat  = GTR(theta[1:6], f0)
  Rmat  = MG94(rmat, theta[7], cod)
  llt = sum(log(expm(Rmat)*cf0)*dat1)
  return(llt)
}



#>>>>>>picture LL
inFile = "../../test_90_species/Results/truePar_100.txt"
inD2   = "../../test_90_species/Results/est100"
main = function(inD2, inFile){
  #Construct codons and its degeneracy
  stp = c(49,51,57)
  cod64 = cbind(rep(DNA_BASES, each = 16),
                rep(DNA_BASES, times = 4, each = 4),
                rep(DNA_BASES, 16))
  cod        = cod64[-stp,]
  codonstrs  = apply(cod, 1, stringr::str_c, collapse = "")            
  syn        = syncodons(codonstrs)
  names(syn) = toupper(names(syn))
  syn        = lapply(syn, toupper)
  
  #setup global
  cod <<- cod
  codonstrs <<- codonstrs
  syn <<- syn
  
  #True parameters, unnormalized
  trueP = read.table(inFile, header=T, sep="")
  tP    = unlist(trueP[n,])
  Pi    = tP[1:4]
  Sigma = tP[5:10]
  Tau   = tP[11]
  omega = tP[12]
  
  #Set up GTR matrix
  gtr = GTR(Sigma, Pi)
  #print(-sum(diag(gtr)*Pi))
  mg94 = MG94(gtr, omega, cod)
  
  
  #Set up mg94 matrix and normalize it.
  Pi2 = sapply(seq(61), function(x){prod(Pi[match(cod[x, ], DNA_BASES)])})  
  Pi2 = Pi2/sum(Pi2)
  #print(sum(Pi2))
  
  #Create Symmetric matrix
  o   = outer(sqrt(Pi2), 1/sqrt(Pi2))      
  s94 = mg94 * o
  #print(s94[1, ]) 
  #print(s94[, 1])
  
  p94 = expm(s94)* t(o) 
  P94 = p94* Pi2        
  #print(rowSums(p94))              
  #print(sum(P94))
  
  #Normlize sigma as well
  #Sigma = Sigma/T
  
  
  #Build omega list for M-step.
  ##Locating all non-syn locations in 64*64 R matrix.
  omega.id = c()
  for (i in 1:61) {
    for (j in 1:61) {
      if((i != j) && 
         (sum(cod[i, ] != cod[j, ]) == 1) && 
         (!(codonstrs[j] %in% syn[[codonstrs[i]]])) ){
        omega.id = c(omega.id, i, j)          
      }
    }
  }
  omega.id = split(omega.id, ceiling(seq_along(omega.id) / 2))
  
  #Build sigma list for M-step
  ##Locating all 6-sigma locations in 64*64 R matrix.
  ii = GTR(1:6, rep(1,4))
  diag(ii) = 0
  I = matrix(0, 61, 61)
  for (i in 1:61) {
    for (j in 1:61) {
      if(i == j){
        I[i, j] = 0
      }else if(sum(cod[i, ] != cod[j, ]) > 1){
        I[i, j] = 0
      }else{
        pos = which(cod[i, ] != cod[j, ])
        x   = which(DNA_BASES == cod[i, pos])
        y   = which(DNA_BASES == cod[j, pos])
        I[i, j] = ii[x, y]
      }
    }
  }
  sigma.id = sapply(1:6, function(x){which(I == x)})
  
  ssize = n*10^5
  #Create sample data
  set.seed(8088)
  
  dat  = sample(61*61, ssize, replace=TRUE, prob = P94)
  dat  = table(dat)
  dat  = as.data.frame(dat)
  id1  = as.numeric(as.vector(dat[[1]]))
  id2  = as.numeric(as.vector(dat[[2]]))
  dat1 = matrix(0, 61, 61)
  for (i in 1:length(id1)) {
    dat1[id1[i]] = id2[i]
  }
  #dat1 = 10000 * P94 >>>>>>>>>>>>>>>>>>
  
  
  ##empirical f
  f1     = rowSums(dat1) + colSums(dat1)  
  f1.id  = sapply(seq(61), function(x){match(cod[x, ], DNA_BASES)})
  base.f = c()
  for (i in seq(4)) {
    base.id = sapply(seq(61), function(x){length(which(f1.id[, x] == i))})
    base.f  = c(base.f, sum(base.id*f1))
  }
  f  = base.f/sum(base.f)
  #print(sum(f)) 
  
  #em to obtain the est of f0
  f0=f
  sum_61 = sum(f1)
  for (i in 1:20) {
    Pi_stp = c(f0[1]^2*f0[4],f0[1]*f0[3]*f0[4],f0[3]*f0[1]*f0[4]) 
    Pi_61  = 1 - sum(Pi_stp)
    
    sum_stp  = sum_61/Pi_61 - sum_61      
    stp_norm = Pi_stp/sum(Pi_stp) *sum_stp
    
    #print the ll
    cf0 = sapply(seq(61), function(x){prod(f0[match(cod[x, ], DNA_BASES)])})
    ll  = sum(f1*log(cf0))- sum_61*log(Pi_61)                         
    #cat(sprintf("%i: %.6f\n", i, ll))
    
    denom = 3*(sum_61+sum_stp)
    
    f0 = c(base.f[1]+2*stp_norm[1]+stp_norm[2]+stp_norm[3], 
           base.f[2],
           base.f[3]+stp_norm[2]+stp_norm[3],
           base.f[4]+sum_stp)/denom
    #print(sum(f0))
  }
  cf0 = cf0/sum(cf0)
  
  #>>global
  f0  <<- f0
  cf0 <<- cf0
  
}

Files2 = list.files(inD2, pattern='.5e', full.names=TRUE)
dlist2 = list()
for (i in 1:length(Files2)) {
  Jtmp      = fromJSON(Files2[i])
  dlist2[[i]] = Jtmp$p.intermed
}
for (j in 1:length(dlist2)) {
  emat = dlist2[[j]]
  Emat = cbind(emat[,c(10,7,9)])
  xgrid = seq(min(Emat[,1]),max(Emat[,1]),length.out=100)
  ygrid = seq(min(Emat[,2]),max(Emat[,2]),length.out=100)
  contour(x=xgrid,y=ygrid,z=emat[,9])
}




