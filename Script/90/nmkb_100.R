suppressPackageStartupMessages(library(Biostrings)) 
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(seqinr))
suppressPackageStartupMessages(library(R.utils))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(expm))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(matlib))
suppressPackageStartupMessages(library(jsonlite))
library(dfoptim)

#setwd("~/Dropbox (ASU)/Indel_project/Script")
#getwd()

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

#-Log-likelihood
LL = function(theta){
  rmat  = GTR(theta[1:6], f0)
  Rmat  = MG94(rmat, theta[7], cod)
  -sum(log(expm(Rmat)*cf0)*dat1)
}


#>>>>>>>>>>>>>>>>>>>
main = function(Name,inFile){
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
  
  #indicators
  #Name   = "Results/nmkb100/1.5e.json"
  #inFile = "../test_90_species/Results/truePar_100.txt"
  n     = as.numeric(str_extract(basename(Name), "[^.]+"))
  trueP = read.table(inFile, header=T, sep="")
  tP    = unlist(trueP[n,])
  
  #True parameters, unnormalized
  Pi    = tP[1:4]
  Sigma = tP[5:10]
  Tau   = tP[11]
  omega = tP[12]
  
  gtr  = GTR(Sigma, Pi)
  mg94 = MG94(gtr, omega, cod)
  
  #Set up mg94 matrix and normalize it.
  Pi2 = sapply(seq(61), function(x){prod(Pi[match(cod[x, ], DNA_BASES)])})  
  Pi2 = Pi2/sum(Pi2)
  print(sum(Pi2))
  
  #Create Symmetric matrix
  o   = outer(sqrt(Pi2), 1/sqrt(Pi2))      
  s94 = mg94 * o
  
  p94 = expm(s94)* t(o) 
  P94 = p94* Pi2        
  print(sum(P94))
  
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
  
  ##########################################################Part III simulate data and run EM.
  ssize = 10^5
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
  
  #empirical f
  f1     = rowSums(dat1) + colSums(dat1)  
  f1.id  = sapply(seq(61), function(x){match(cod[x, ], DNA_BASES)})
  base.f = c()
  for (i in seq(4)) {
    base.id = sapply(seq(61), function(x){length(which(f1.id[, x] == i))})
    base.f  = c(base.f, sum(base.id*f1))
  }
  f  = base.f/sum(base.f)
  print(sum(f)) 
  
  #init f0
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
    cat(sprintf("%i: %.6f\n", i, ll))
    
    denom = 3*(sum_61+sum_stp)
    
    f0 = c(base.f[1]+2*stp_norm[1]+stp_norm[2]+stp_norm[3], 
           base.f[2],
           base.f[3]+stp_norm[2]+stp_norm[3],
           base.f[4]+sum_stp)/denom
  }
  
  #sim LL
  ll.sim  = sum(log(expm(mg94)*Pi2)*dat1)
  #emp LL
  r1      = GTR(Sigma, f0)
  cf0     = cf0/sum(cf0)
  R1      = MG94(r1, omega, cod)
  ll.emp  = sum(log(expm(R1)*cf0)*dat1)
  
  if(ll.sim<ll.emp){
    print("Yes!")
  }else{
    print("come on man!")
  }
  
  #init sigma
  #create nuc transition matrix [4-fold-degeneracy-codons only]
  fourD_id = which(sapply(syn, function(x){length(x) == 4}))
  k=1
  ndat = matrix(0,4,4)
  while (k<length(fourD_id)) {
    i = j = fourD_id[k:(k+3)]
    ndat = ndat+dat1[i,j]
    k=k+4
  }
  fn = rowSums(ndat) + colSums(ndat)        #neutral freq
  fn = fn/sum(fn)
  print(sum(fn))
  
  #symmetric average of dat
  sdat  = matrix(0,4,4)
  diag(sdat) = diag(ndat)
  sdat  = (ndat + t(ndat))/2 
  sf    = colSums(sdat)/sum(sdat)
  Dhat  = diag(sf)
  phat  = t(sdat/colSums(sdat))
  
  eigP  = eigen(phat)
  Athat = eigP$vectors %*% diag(log(eigP$values)) %*%  inv(eigP$vectors)
  that  = -sum(diag(Athat*Dhat))
  Ahat  = t(Athat)/that
  s     = Ahat[lower.tri(t(Ahat))]
  #print(that)    
  #print(Ahat)
  if(any(s<0)){
    cat(sprintf("neg-s0:  %i\n",n))
    #s = runif(6)
    mean(diag(Ahat))
    i = c(1,2,3,4)
    j = c(3,4,1,2)
    ts = mean(Ahat[cbind(i,j)])   #transitions
    diag(Ahat)=0 
    tv = (sum(Ahat)-4*ts)/8       #transversions
    
    hky = matrix(0,4,4)
    hky[lower.tri(hky)] = c(tv,ts,tv,tv,ts,tv)
    hky = hky + t(hky)
    hky = t(hky*sf)
    diag(hky) = -rowSums(hky)
    s = hky[lower.tri(t(hky))]
  }
  
 
  
  #init omega
  ##obs non-syn change. 
  obs.non = sum(sapply(omega.id, function(x){dat1[x[1],x[2]]}))
  
  ##exp non-syn change
  gtr2  = GTR(s, f0)
  mg942 = MG94(gtr2, 1, cod)
  P942  = expm(mg942)
  #print(sum(P942))
  
  exp.non = 0
  for (i in 1:61) {
    ith.non = omega.id[which(sapply(omega.id,'[[',1)==i)]
    exp.non = exp.non + f1[i]*sum(sapply(ith.non, function(x){P942[x[1],x[2]]}))
  }
  w = obs.non/(exp.non/2)
  
  
  #init p >>>>>>>>>>>>>>>>>>>>
  f0       <<- f0
  cf0      <<- cf0
  dat1     <<- dat1
  omega.id <<- omega.id
  
  
  p0 = c(s,w)
  if(any(p0>1)){
    len = length(which(p0>1))
    p0[which(p0>1)] = runif(len)
  }
  
  tm = system.time({
    pb = nmkb(fn=LL, par=p0, lower=0.0, upper=2.5, control=list(tol=1e-4,trace=TRUE)) 
  })
  
  # print(tm)
  # print(pb$message)
  # print(pb$par)
  
  rmat = GTR(pb$par[1:6],f0)
  tv   = -sum(diag(rmat)*f0) 
  
  pb$par = c(pb$par,tv,f0)
  pj = toJSON(pb)
  
  write(pj,Name)
}


#################################
args = commandArgs(trailingOnly=T)
main(args[1],args[2])


