# Apply the importance sampling Wi = g/f 
# g(x)--ziqi's phase_coati model; f(x)--juan's coatiM 

# A -- ancestor, B -- descendent
# insertion and deletion length distribution seperated.
# inserting nucleotides are ignored like deletions.


suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(seqinr))
suppressPackageStartupMessages(library(stringi))
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(jsonlite))
suppressPackageStartupMessages(library(dfoptim))
suppressPackageStartupMessages(library(tidyverse))

#setwd("~/Dropbox (ASU)/Indel_project/chapter3")

######################################
#Count the nucleotide freq
count_freq = function(input){
  nuc.count = 0
  for (i in input){
    dna       = readDNAStringSet(i)
    nuc.count = nuc.count + oligonucleotideFrequency(dna,width=1)
  }
  nuc.freq  = colSums(nuc.count)/sum(nuc.count)
  if(sum(nuc.freq)==1){
    return(nuc.freq)
  }else{
    print("The nucleotide frequency sums up not 1!")
  }
}


#-Log-likelihood of subs
LL_min  = function(theta){
  nuc_codon_freq = init_f(cod,datw,64)
  fw  = nuc_codon_freq[[1]]
  cfw = nuc_codon_freq[[2]]
  
  rmat = GTR(theta[1:6], fw)
  Rmat = MG94(rmat,theta[7],cod,codonstrs,syn,64)
  -sum(log(expm(Rmat)*cfw)*datw)
}


##-Log-likelihood of indels
LL_min_Z = function(gv){
  g012 = gv[1:6]
  wz   = gv[7]
  #gap length LL
  if(sum(lgw)==0){#no gap
    scoreE = 0
    lk     = rep(0,4)
  }else{
    lk     = lgw/sum(lgw)
    lkreal = lk[which(lk!=0)]
    lgreal = lgw[which(lgw!=0)]
  }
  #gap LL
  scoreG = sum(phasew[,1]*log(g012*wz) + phasew[,2]*log(g012))
  
  #non-gap LL
  k = rowSums(t(t(pw)*lk))
  
  scoreM = sum(iw[1,]*log(1-g012[1] + (1-wz)*g012[1]*k)) + 
    sum(iw[2,]*log(1-g012[2] + (1-wz)*g012[2]*k)) + 
    sum(iw[3,]*log(1-g012[3] + (1-wz)*g012[3]*k)) + 
    sum(iw[4,]*log(1-g012[4] + (1-wz)*g012[4]*k)) + 
    sum(iw[5,]*log(1-g012[5] + (1-wz)*g012[5]*k)) + 
    sum(iw[6,]*log(1-g012[6] + (1-wz)*g012[6]*k)) 
  
  -(scoreG + scoreM)
}


##grab the number of all 16 indicator patterns
grabI = function(Iz0,pw){
  iw0 = rep(0,16)
  for (p in 1:16) {
    lenI = length(which(sapply(1:nrow(Iz0),function(x){all.equal(Iz0[x,], pw[p,])==TRUE})==TRUE))
    if(lenI==0){
      next
    }else{
      iw0[p] = lenI
    }
  }
  return(iw0)
}  



#setwd("~/Dropbox (ASU)/Indel_project/chapter4")
##############################################



main = function(dir, ouFile){
  
  sub1 = "../Script/sources/"
  sub2 = "../Script/chapter3/"
  source(paste0(sub1,"codon_call.R"))
  source(paste0(sub1,"gtr.R"))
  source(paste0(sub1,"mg94.R"))
  source(paste0(sub1,"init_f.R"))
  source(paste0(sub1,"LL.R"))
  source(paste0(sub2,"haha1.R"))
  
  #construct codons and its degeneracy
  co.res    = codon_call()
  codons    = co.res[[1]]
  codonstrs = co.res[[2]]
  syn       = co.res[[3]]
  
  #nmkb method
  cod      <<- codons
  codonstrs<<- codonstrs
  syn      <<- syn
  
  
  ################################################
  #inD  = "Gs/1"
  Files = list.files(inD, full.names=T)
  index = as.numeric(str_extract(basename(Files),'[^.]+'))
  Files = Files[order(index)]
  Files = Files[order(index)]
  K     = length(Files)
  
  #generate an alignment list
  Alist = list()
  for (i in 1:K) {
    Alist[[i]] = readBStringSet(Files[i])
  }
  
  #summary statistic
  # Gm           = matrix(0,K,6)
  # Mm           = matrix(0,K,6)
  # codon_array  = array(0,c(64,64,K))
  # llz          = rep(0,K)
  # avg.gap      = matrix(0,K,2)
  
  
  #Set up initial parameters
  #indel rates + omega_z
  set.seed(8088)
  f0   = count_freq(Files)
  p0   = runif(7)
  g0   = runif(7)
  e0   = rep(0.8,2)
  
  
  #summary statistic
  #>seq likelihood
  llz  = c()                      #>>>>>>>>>>>>>>>>>>>>>>>>>>
  #>gap total length
  gtl  = c()                      #>>>>>>>>>>>>>>>>>>>>>>>>>>
  #>number of gaps 
  gnum = matrix(0,K,4)            #>>>>>>>>>>>>>>>>>>>>>>>>>>
  #>gap array
  gArr = array(0,c(6,2,K))
  #>codon array
  codArr= array(0,c(64,64,K))
  #>Indicator list of Non-gap edge.
  Iz   = list()
  #>phase_wz_matrix
  pw = as.matrix(expand.grid(c(0,1),c(0,1),c(0,1),c(0,1)))  
  pw = unname(pw)
  
  ############################>>>>>
  #``
  # Jpar = list(Pi,p0,mean(g0),ext0)
  # JJ   = toJSON(Jpar)  
  #``
  
  
  
  #Iteration through the black box
  iter  = 0
  ll.old= rep(0,K)
  llv   = matrix(0,5,K)
  
  repeat{
    iter=iter+1
    if(iter>5) {break}
    
    
    avg.gap = 1/(1-e0)
    if(any(avg.gap<3)){
      print("Average gap length less than 3!")
      break
    }
    e3 = 1-1/(avg.gap/3)
    for(i in 1:K){#E step
      A = Alist[[i]]
      #summary-stat
      res1        = ziqi_prob_wz(A,g0,e3,p0,f0)
      llz[i]      = res1[[1]]
      Iz[[i]]     = res1[[2]]
      gArr[,,i]   = res1[[3]][,1:2]
      gnum[i,]    = res1[[4]]
      gtl[i]      = res1[[5]]
      codArr[,,i] = res1[[6]]
      print(i)
    }
    
    ##set up tolerence
    # cat(sprintf("%i: %.6f\n", iter, mean(abs(llz-ll.old))))
    # llv[iter,]=llz
    # ll.old    = llz
    
    
    
    #M step: para estimates
    ## assume juan's weight is llJ
    #llj = rnorm(K,llz,.1)
    #````
    # cmd1 = "bash coatim-run.sh"
    # system(cmd1)
    # cmd2 = "bash coatim-weight.sh"
    # llj  = system(cmd2, intern=TRUE)
    #````
    
    ##Cal. the weight
    Wi = rep(1/K,K)
    
    ##nmkb for p0
    #>weighted codon matrix
    datw = matrix(0,64,64)
    for (j in 1:K) {
      dat.wei = Wi[j]*codArr[,,j]/K
      datw = datw + dat.wei
    }
    datw <<- datw
    
    pb = nmkb(fn=LL_min, par=p0, lower=0.0, upper=2.5, control=list(tol=1e-4,trace=F)) 
    if(pb$convergence != 0){
      cat("Warning: failed convergence!")
    }else{
      p0 = pb$par
    }
    
    ##nmkb for g0
    #>weighted gap length
    glw     = mean(Wi*gtl)
    gnw     = mean(Wi*rowSums(gnum))
    Est.ext = (glw-gnw)/glw 
    
    #>weighted gap phases
    phasew = matrix(0,6,2)
    for (j in 1:K) {
      phasew = phasew + gArr[,,j]*Wi[j]
    }
    
    #>weighted gap numbers
    lgw = colSums(gnum*Wi)/K
    
    
    
    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    iw = matrix(0,6,16) #based on pw patterns.
    for (h in 1:K) {
      ipos1 = seq(1,nrow(Iz[[h]][,,1])) %% 3
      ipos2 = seq(1,nrow(Iz[[h]][,,2])) %% 3
      
      #insertion of phase 0,1,2
      Iz0 = Iz[[h]][which(ipos1==1),,1]
      Iz1 = Iz[[h]][which(ipos1==2),,1]
      Iz2 = Iz[[h]][which(ipos1==0),,1]
      
      #Deletion of phase 0,1,2
      Iz0D = Iz[[h]][which(ipos1==1),,2]
      Iz1D = Iz[[h]][which(ipos1==2),,2]
      Iz2D = Iz[[h]][which(ipos1==0),,2]
      
      iw[1,] = iw[1,] + grabI(Iz0,pw)/K
      iw[2,] = iw[2,] + grabI(Iz1,pw)/K
      iw[3,] = iw[3,] + grabI(Iz2,pw)/K
      iw[4,] = iw[4,] + grabI(Iz0D,pw)/K
      iw[5,] = iw[5,] + grabI(Iz1D,pw)/K
      iw[6,] = iw[6,] + grabI(Iz2D,pw)/K
    }
    
    
    phasew<<-phasew
    lgw   <<-lgw
    pw    <<- pw
    iw    <<- iw
    
    gb  = nmkb(fn=LL_min_Z, par=g0, lower=0.0, upper=1.0, control=list(tol=1e-4,trace=F)) 
    if(gb$convergence != 0){
      cat("Warning: failed convergence!")
    }else{
      g0 = gb$par
    }
    
    
    ##cal. the tau
    nuc_f = init_f(cod,datw,64)
    fnew  = nuc_f[[1]]
    rmat  = GTR(p0[1:6],fnew)
    tnew  = -sum(diag(rmat)*fnew)
  }
  
  #Cal. the extension prob (T-G)/T
  #will the total length be weighted as well? 
  
  ##pipe out the json file.
  ##confirm with gillespie simu: 
  ##p0[1:6]/tnew vs Sigma/T
  ##p0[7] vs W
  ##fnew vs Pi
  ##Est.ext vs ext  
  ##tnew vs brlen    
  ##r1 vs (g0[1]+g0[4])/2, (g0[2]+g0[5])/2, (g0[3]+g0[6])/2
  ##g0[7] vs Wz
  
  
  Jpar = list(fnew,pnew,mean(gnew),ext.new)
  JJ   = toJSON(Jpar)  
  write(JJ,"../../chapter3/JJ.json")
  JJ0  = fromJSON(file="../../chapter3/JJ.json")
  
  
  
}