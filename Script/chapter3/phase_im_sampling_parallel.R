#Run parallel in phase importance sampling. 

suppressPackageStartupMessages(library(doParallel))
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

#Read sample
read_sample = function(D,K){
  Alist = list()
  for (i in 1:K) {
    Alist[[i]] = DNAStringSet(c(D$Seq1[i],D$Seq2[i]))
  }
  return(Alist)
}

# -Log-likelihood of
LL_min  = function(theta){
  nuc_codon_freq = init_f(cod,DD,64)
  fw  = nuc_codon_freq[[1]]
  cfw = nuc_codon_freq[[2]]
  
  rmat = GTR(theta[1:6], fw)
  Rmat = MG94(rmat,theta[7],cod,codonstrs,syn,64)
  -sum(log(expm(Rmat)*cfw)*DD)
}



##############################################
main = function(inD,ouD,ouF1,ouF2,nc){
  
  sub1 = "../Script/sources/"
  sub2 = "../Script/chapter3/"
  source(paste0(sub1,"codon_call.R"))
  source(paste0(sub1,"gtr.R"))
  source(paste0(sub1,"mg94.R"))
  source(paste0(sub1,"init_f.R"))
  source(paste0(sub1,"LL.R"))
  source(paste0(sub2,"phase_indel_prob.R"))
  
  #construct codons and its degeneracy
  co.res    = codon_call()
  codons    = co.res[[1]]
  codonstrs = co.res[[2]]
  syn       = co.res[[3]]
  
  #nmkb method
  cod      <<- codons
  codonstrs<<- codonstrs
  syn      <<- syn
  
  
  ##############################>>>>>>>>>>>>>>>>>>>
  # ouD   = "Json/"
  # inD   = "fasta"
  Files = list.files(inD, full.names=T)
  index = as.numeric(str_extract(basename(Files),'[^.]+')) #rectify the order
  Files = Files[order(index)]
  n     = length(Files)
  
  
  #Initial parameters
  ##gap
  #g0   = runif(6) 
  #ext0 = runif(1)
  #p0   = runif(7)
  # ##sub
  # brlen=0.01
  # s0   = c(0.7135135,2.9447236,0.3246753,1.1608040,2.9123377,0.6428571)*brlen
  # p0   = c(s0,0.1)
  # g0   = rep(0.02,6)
  # ext0 = 0.66
  
  #Initial parameters
  set.seed(8088)
  f0   = count_freq(Files)
  p0   = rep(0.1,7)
  g0   = rep(0.01,6)
  e0   = rep(0.5,2)
  
  rmat = GTR(p0[1:6],f0)
  t0   = -sum(diag(rmat)*f0) 
  
  cat("the number of cores is: ",detectCores(),"\n")
  ncores = as.numeric(nc)
  registerDoParallel(cores=ncores)
  cat("ncores: ",ncores,"\n")
  
  #combine func
  comb = function(x, ...) {
    lapply(seq_along(x),
           function(i) c(x[[i]], lapply(list(...), function(y) y[[i]])))
  }
  
  #Iterate through the black box
  lljv   = matrix(0,n,10)
  llzv   = matrix(0,n,10)
  maxw   = matrix(0,n,10)
  
  N      = matrix(0,n,6)
  M      = matrix(0,n,6)
  E      = matrix(0,n,4)
  colnames(N) = c('i0','i1','i2','d0','d1','d2')
  colnames(M) = c('no_i0','no_i1','no_i2','no_d0','no_d1','no_d2')
  colnames(E) = c('i.len','i.num','d.len','d.num')
  
  ev     = matrix(e0,1,2)
  tv     = c()
  pv     = matrix(p0,1,7)
  gv     = matrix(g0,1,6)
  Wv     = matrix(0,100,100)
  
  
  iter=1
  tm = system.time({
    repeat{
      Results = foreach (j = 1:n, .combine='comb', .multicombine=T, 
                         .init=list(list(),list(),list(),list(),list(),list(),list(),list())) %dopar% {#E-step
                           
                           #>merge the coati-sampler
                           input = Files[j]
                           ouJ   = paste0(ouD,j,'.json')
                           cmd   = paste("bash ../Script/chapter3/coati_sampler.sh", input,ouJ, 
                                         mean(g0),mean(e0),f0[1],f0[2],f0[3],f0[4], p0[1],p0[2], p0[3],p0[4],p0[5],p0[6],p0[7],t0,sep=' ')
                           system(cmd)
                           Input.tmp    = fromJSON(ouJ)
                           Dat.tmp      = Input.tmp %>% unpack(aln)
                           K            = nrow(Dat.tmp)
                           llj          = Dat.tmp$log_weight
                           llj2         = exp(llj)/sum(exp(llj))
                           Alist        = read_sample(Dat.tmp,K)
                           
                           
                           ##summary statistic
                           Gm           = matrix(0,K,6)
                           Mm           = matrix(0,K,6)
                           codon_array  = array(0,c(64,64,K))
                           llz          = rep(0,K)
                           gl           = matrix(0,K,2)
                           
                           ##run sampler in ziqi's model
                           for(i in 1:K){
                             A      = Alist[[i]]
                             res1   = ziqi_prob(A,g0,e0,p0,f0)
                             Gm[i,] = res1[[1]]
                             Mm[i,] = res1[[2]]
                             codon_array[,,i] = res1[[3]]
                             llz[i] = res1[[4]]
                             gl[i,] = res1[[5]]
                           }
                           
                           
                           #cal. the weight
                           llz.avg = mean(llz)
                           llz.diff= llz-llz.avg
                           llz2    = exp(llz.diff)/sum(exp(llz.diff))
                           Wi      = llz2/llj2
                           Wi      = Wi/sum(Wi)
                           
                           ##weighted codon
                           datw = matrix(0,64,64)
                           for (j in 1:K) {
                             dat.wei = Wi[j]*codon_array[,,j]
                             datw    = datw + dat.wei
                           }
                           
                           #ouput a list 
                           x1 = mean(llj)
                           x2 = mean(llz)
                           x3 = max(Wi)
                           x4 = Wi
                           x5 = colSums(Wi*Gm)
                           x6 = colSums(Wi*Mm)
                           x7 = c(sum(Wi*gl[,1]),sum(Wi*rowSums(Gm[,1:3])),sum(Wi*gl[,2]),sum(Wi*rowSums(Gm[,4:6])))
                           x8 = datw
                           list(x1,x2,x3,x4,x5,x6,x7,x8)
                         }
      
      
      ##decompose the list
      lljv[,iter] = unlist(Results[[1]])
      llzv[,iter] = unlist(Results[[2]])
      maxw[,iter] = unlist(Results[[3]])
      
      D = matrix(0,64,64)  #reset to zero 
      for (j in 1:n) {#row1=fasta1
        Wv[j,] = Results[[4]][[j]]
        N[j,]  = Results[[5]][[j]]
        M[j,]  = Results[[6]][[j]]
        E[j,]  = Results[[7]][[j]]
        D      = D + Results[[8]][[j]]/n
      }
      
      
      #M step: para estimates
      ##gap opening 
      nnew = colMeans(N)
      mnew = colMeans(M)
      gnew = nnew/(nnew+mnew)
      
      print((gnew[1]+gnew[4])/2)
      print((gnew[2]+gnew[5])/2)
      print((gnew[3]+gnew[6])/2)
      
      ##gap extension
      g.new = colMeans(E)
      enew  = c((g.new[1]-g.new[2])/g.new[1], (g.new[3]-g.new[4])/g.new[3])
      
      DD <<- D
      pb = nmkb(fn=LL_min, par=p0, lower=0, upper=1, control=list(tol=1e-5,trace=F))  #change the tolerance
      if(pb$convergence != 0){
        cat("Warning: failed convergence!")
      }else{
        pnew = pb$par
      }
      
      ##cal. the tau
      nuc_f = init_f(cod,DD,64)     
      fnew  = nuc_f[[1]]
      rmat  = GTR(pnew[1:6],fnew)
      tnew  = -sum(diag(rmat)*fnew)
      
      print(pnew[1:6]/tnew)
      print(pnew[7])
      print(tnew)
      
      #record the pars
      gv   = rbind(gv,gnew)
      pv   = rbind(pv,pnew)
      ev   = rbind(ev,enew)
      tv   = c(tv,tnew)
      
      #rmse tolerance
      p=c(g0,e0,p0)
      q=c(gnew,enew,pnew)
      delta= (q-p)^2
      rmse = sqrt(mean(delta))
      if(rmse<=1e-3){
        break
      }else{
        cat(sprintf("iter:%i, rmse:%.6f\n",iter, rmse))
        iter=iter+1
        f0  = fnew
        p0  = pnew
        g0  = gnew
        e0  = enew
        t0  = tnew
      }
    }
  })
  
  cat(sprintf("Iter:%i\t time:%.2f mins", iter, tm[3]/60))
  # ouF1 = "pise.json"
  # ouF2 = "weight_matrix.txt" 
  par_lst = list('nuc.freq'=fnew,'gap.openning'=gnew/tnew, 'gap.extension'=enew, 'sigmas'=pnew[1:6]/tnew, 'omega'=pnew[7], 'branch.length'=tnew)
  par_est = toJSON(par_lst)
  write(par_est,ouF1)
  
  #recording the weight matrix
  write.table(Wv,ouF2,quote=F,row.names=F,col.names=F,sep='\t')






}

########################################
args = commandArgs(trailingOnly=T)
main(args[1],args[2],args[3],args[4],args[5])
