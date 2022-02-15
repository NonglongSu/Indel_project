suppressPackageStartupMessages(library(Biostrings)) 
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(seqinr))
suppressPackageStartupMessages(library(R.utils))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(expm))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(matlib))
suppressPackageStartupMessages(library(jsonlite))
suppressPackageStartupMessages(library(SQUAREM))

main = function(num_dim){
  
  setwd("/home/tamade/Dropbox (ASU)/Indel_project/Script")
  
  source("sources/codon_call.R")
  source("sources/gtr.R")
  source("sources/mg94.R")
  source("sources/omega_coor.R")
  source("sources/sigma_coor.R")
  source("sources/pseudo_data.R")
  source("sources/init_f.R")
  source("sources/LL.R")
  source("sources/init_sigma.R")
  source("sources/init_omega.R")
  source("sources/phylo_em.R")
  
  #example
  Name   = "Results/rms/1.json"
  inFile = "../test_90_species/Results/truePar_100.txt"
  n     = as.numeric(str_extract(basename(Name), "[^.]+"))
  trueP = read.table(inFile, header=T, sep="")
  tP    = unlist(trueP[n,])
  
  #True parameters, unnormalized
  Pi    = tP[1:4]
  Sigma = tP[5:10]
  Tau   = tP[11]
  omega = tP[12]
  
  #61 or 64?
  #num_dim = '61'
  num = as.numeric(num_dim)
  print(num)
  
  #>
  cdc = codon_call(num)
  cod = cdc[[1]]
  codonstrs = cdc[[2]]
  syn = cdc[[3]]
  
  #>
  gtr = GTR(Sigma,Pi)
  
  #>
  
  mg94 = MG94(gtr,omega,cod,codonstrs,syn,num)

  #>
  omega.id = omega_coor(cod, codonstrs, syn, num)
  sigma.id = sigma_coor(cod, num)
  
  ##generate pmat
  Pi2 = sapply(seq(num), function(x){prod(Pi[match(cod[x, ], DNA_BASES)])})  
  Pi2 = Pi2/sum(Pi2)
  
  o   = outer(sqrt(Pi2),1/sqrt(Pi2))      
  s94 = mg94*o                            #symmetric matrix
  p94 = expm(s94)*t(o) 
  P94 = p94*Pi2        
  print(sum(P94))
  
  
  #>
  set.seed(8088)
  dat = pseudo_data(1e+5,P94,num)
  
  #>
  nuc_codon_freq = init_f(cod,dat,num)
  f0  = nuc_codon_freq[[1]]
  cf0 = nuc_codon_freq[[2]]
  print(sum(f0))
  print(sum(cf0))
  
  #test if sim.LL<emp.LL
  ll.sim = sum(log(expm(mg94)*Pi2)*dat)
  ll.emp = LL(cod,codonstrs,syn,c(Sigma,omega),f0,cf0,dat,num) 
  if(ll.sim<ll.emp){
    print("Yes!")
  }else{
    print("come on man!")
  }
  
  #>
  s = init_sigma(syn,dat)
  w = init_omega(cod,codonstrs,syn,dat,f0,s,omega.id,num)
  
  #running SQUAREM
  num      <<- num
  cod      <<- cod
  codonstrs<<- codonstrs
  syn      <<- syn
  f0       <<- f0
  cf0      <<- cf0
  dat      <<- dat
  omega.id <<- omega.id
  sigma.id <<- sigma.id
  
  
  p0   = c(s,w)
  tm4b = system.time({
    ps4b = fpiter(par=p0, fixptfn=phylo_em, control=list(tol=1e-4, trace=TRUE, intermed=TRUE))
  })
  pd  = ps4b$p.intermed
  
  llv = c()
  tv  = c()
  for (i in 1:nrow(pd)) {
    rmat   = GTR(pd[i,1:6],f0)
    tv[i]  = -sum(diag(rmat)*f0) 
    llv[i] = LL(cod,codonstrs,syn,pd[i,1:7],f0,cf0,dat,num)
  }
  
  nuc.mat         = matrix(f0,nrow(pd),4,byrow=T)
  ps4b$p.intermed = cbind(pd,llv,tv,nuc.mat)
  ps4B            = toJSON(ps4b)
  
  #>>>display
  print(pd)
  
}
  
############################
args=commandArgs(trailingOnly=T)
main(args[1])

#usage: Rscript --vanilla ../Script/haha_sources_test.R 61

  








