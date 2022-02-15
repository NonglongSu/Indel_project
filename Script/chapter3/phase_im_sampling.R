# Take the coatiM (viterbi algo) and weight score as input
# Apply the importance sampling Wi = g/f
# g(x)--ziqi's phase_coati model; f(x)--juan's coatiM
#
# A -- ancestor, B -- descendent
# insertion and deletion length distribution combined.
# inserting bases are ignored since they are treated as deletions.


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
  dna       = readDNAStringSet(input)
  nuc.count = oligonucleotideFrequency(dna,width=1)
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
  nuc_codon_freq = init_f(cod,datw,64)
  fw  = nuc_codon_freq[[1]]
  cfw = nuc_codon_freq[[2]]

  rmat = GTR(theta[1:6], fw)
  Rmat = MG94(rmat,theta[7],cod,codonstrs,syn,64)
  -sum(log(expm(Rmat)*cfw)*datw)
}



##############################################



main = function(input,output,ouF){

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
  #dir     = "hmr/coatiM/Raw"
  # dir   = "haha2"
  # Files = list.files(dir, full.names=T)
  # K     = length(Files)
  # Alist = list()
  # for (i in 1:K) {Alist[[i]] = readBStringSet(Files[i])}

  #gillespie.1 simu
  # input = "sim.fasta"
  # output= "Json/out0.json"
  
  #generate an alignment list
  input = "example.1.fasta"
  #output= "example.1.json"

  # Input = fromJSON(output)
  # Dat   = Input %>% unpack(aln)
  # K     = nrow(Dat)
  # llj   = Dat$log_weight
  # Alist = read_sample(Dat,K)
  # 
  # llj1  = Dat$weight
  # llj2  = llj1/sum(llj1)
  # print(sum(llj2))

  #summary statistic
  # Gm           = matrix(0,K,6)
  # Mm           = matrix(0,K,6)
  # colnames(Gm) = c('ins0','ins1','ins2','del0','del1','del2')
  # colnames(Mm) = c('no_ins0','no_ins1','no_ins2','no_del0','no_del1','no_del2')
  # codon_array  = array(0,c(64,64,K))
  # 
  # ##seq likelihood, gap total length
  # llz = gtl = c()

  #Initial parameters
  ##gap
  set.seed(8088)
  #g0   = runif(6) 
  #ext0 = runif(1)
  #p0   = runif(7)
  # ##sub
  # brlen=0.01
  # s0   = c(0.7135135,2.9447236,0.3246753,1.1608040,2.9123377,0.6428571)*brlen
  # p0   = c(s0,0.1)
  # g0   = rep(0.02,6)
  # ext0 = 0.66
  
  # f0   = count_freq(input)
  # g0   = runif(6)
  # ext0 = runif(1)
  # p0   = runif(7)
  
  f0   = count_freq(input)
  p0   = rep(0.1,7)
  g0   = rep(0.01,6)
  ext0 = runif(1)
  
  rmat = GTR(p0[1:6],f0)
  tnew = -sum(diag(rmat)*f0) 
  
  

  #Iterate through the black box
  niter  = 1
  iter   = 0
  lljv   = llzv  =  c()
  maxw   = c()
  ev     = c()
  tv     = c()
  pv     = matrix(p0,1,7)
  gv     = matrix(g0,1,6)
  
  repeat{
    if(iter>niter){break}
    iter=iter+1
    
    
      #>merge the coati-sampler
      ouJ  = paste0('Json/out',iter,'.json')
      cmd  = paste("bash ../Script/chapter3/coati_sampler.sh", input,ouJ, 
                     mean(g0),ext0,f0[1],f0[2],f0[3],f0[4], p0[1],p0[2], p0[3],p0[4],p0[5],p0[6],p0[7],tnew, sep=' ')
      system(cmd)
      Input.tmp = fromJSON(ouJ)
      Dat.tmp   = Input.tmp %>% unpack(aln)
      K         = nrow(Dat.tmp)
      llj       = Dat.tmp$log_weight
      llj1      = Dat.tmp$weight
      llj2      = llj1/sum(llj1)
      Alist     = read_sample(Dat.tmp,K)
      #record coati weight
      lljv      = c(lljv,mean(llj))

      #>summary statistic
      Gm           = matrix(0,K,6)
      Mm           = matrix(0,K,6)
      colnames(Gm) = c('ins0','ins1','ins2','del0','del1','del2')
      colnames(Mm) = c('no_ins0','no_ins1','no_ins2','no_del0','no_del1','no_del2')
      codon_array  = array(0,c(64,64,K))
      ##seq likelihood, gap total length
      llz = gtl = c()
      

    ###independent repeat twice
    jter = 0
    repeat{
      if(jter>1){break}
      jter=jter+1
      
      for(i in 1:K){#E-step
        A      = Alist[[i]]
        res1   = ziqi_prob(A,g0,p0,f0,ext0)
        Gm[i,] = res1[[1]]
        Mm[i,] = res1[[2]]
        codon_array[,,i] = res1[[3]]
        llz[i] = res1[[4]]
        gtl[i] = res1[[5]]
      }
      llzv = c(llzv,mean(llz))
      
      ##set up tolerence
      # llv[iter,]= llz
      # cat(sprintf("%i: %.6f\n", iter, mean(llz-ll.old)))
      # ll.old    = llz
      
      #M step: para estimates
      ##Cal. the weight
      
      #option 1: 
      llz.avg = mean(llz)
      llz.diff= llz-llz.avg
      llz2    = exp(llz.diff)/sum(exp(llz.diff))
      print(sum(llz2))
      
      Wi = llz2/llj2
      Wi = Wi/sum(Wi)
      hist(Wi)
      #option 2:
      #llj  = rnorm(K,llz,.01)
      #llj2 = exp(llj)/sum(exp(llj))
      #Wi   = exp(llz-llj)  
      #Wi   = Wi/sum(Wi)
      
      ##weighted gap                   
      # nnew = colMeans(Wi*Gm)
      # mnew = colMeans(Wi*Mm)
      # g0   = nnew/(nnew+mnew)
      
      nnew = colSums(Wi*Gm)
      mnew = colSums(Wi*Mm)
      g0   = nnew/(nnew+mnew)
      
      print((g0[1]+g0[4])/2)
      print((g0[2]+g0[5])/2)
      print((g0[3]+g0[6])/2)
      
      ##weighted gap length
      # gl.new = mean(Wi*gtl)
      # gn.new = mean(Wi*rowSums(Gm))
      # ext0   = (gl.new-gn.new)/gl.new
      
      gl.new = sum(Wi*gtl)
      gn.new = sum(Wi*rowSums(Gm))
      ext0   = (gl.new-gn.new)/gl.new
      
      
      ##weighted codon
      datw = matrix(0,64,64)
      # for (j in 1:K) {
      #   dat.wei = Wi[j]*codon_array[,,j]/K
      #   datw    = datw + dat.wei
      # }
      for (j in 1:K) {
        dat.wei = Wi[j]*codon_array[,,j]
        datw    = datw + dat.wei
      }
      
      
      
      datw <<- datw
      pb   = nmkb(fn=LL_min, par=p0, lower=0, upper=1, control=list(tol=1e-5,trace=T))  #change the tolerance
      if(pb$convergence != 0){
        cat("Warning: failed convergence!")
      }else{
        p0 = pb$par
      }

      
      ##cal. the tau
      nuc_f = init_f(cod,datw,64)
      f0    = nuc_f[[1]]
      rmat  = GTR(p0[1:6],f0)
      tnew  = -sum(diag(rmat)*f0)
      
      print(p0[1:6]/tnew)
      print(p0[7])
      print(tnew)
      #record the pars
      maxw = c(maxw, max(Wi))
      gv   = rbind(gv,g0)
      pv   = rbind(pv,p0)
      ev   = c(ev,ext0)
      tv   = c(tv,tnew)
    }
  }

  #Cal. the extension prob (T-G)/T
  #will the total length be weighted as well?

  ##pipe out the json file.
  ##confirm with gillespie simu:
  ##p0[1:6]/tnew vs Sigma/T
  ##p0[7] vs W
  ##Pi
  ##ext0 vs ext
  ##tnew vs brlen
  ##r1*brlen vs (g0[1]+g0[4])/2, (g0[2]+g0[5])/2, (g0[3]+g0[6])/2

 # par_lst = list('nuc.freq'=f0,'indel.rate'=g0, 'subs.rate'=p0, 'branch.length'=tnew)
 # par_est = toJSON(par_lst)
 # write(par_est,ouF)
}


################################
args = commandArgs(trailingOnly=T)
main(args[1],args[2],args[3])
