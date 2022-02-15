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

#inF  = "~/Dropbox (ASU)/Indel_project/test_90_species/Results/truePar_100.txt"
inF1 = "~/Dropbox (ASU)/Indel_project/test_90_species/Results/est100/1.5e.json"
inF2 = "~/Dropbox (ASU)/Indel_project/test_90_species/Results/1.5e.dumbEM.json"


setwd("~/Dropbox (ASU)/Indel_project/Script")
source("sources/codon_call.R")
source("sources/gtr.R")
source("sources/mg94.R")
source("sources/LL.R")
source("sources/init_f.R")
source("sources/pseudo_data.R")

# #true LL
tLL = -486885.783818
# 
Jtmp1  = fromJSON(inF1)
Jpar1  = Jtmp1$p.intermed
Jtmp2  = fromJSON(inF2)
Jpar2  = Jtmp2$p.intermed

n1 = 1:nrow(Jpar1)
n2 = 1:nrow(Jpar2)
s1 = seq(length(n1)-1)
s2 = seq(length(n2)-1)
plot(n2,Jpar2[,9],col="#000000", ylim=c(-530000,-480000),
     xlab="Iterations",ylab="Log-likelihood")
points(n1,Jpar1[,9],col="#56B4E9",pch=6)
arrows(x0=n1[s1],y0=Jpar1[s1,9],x1=n1[s1+1],y1=Jpar1[(s1+1),9], 
       col="#56B4E9",lwd=2,length=0.05)
arrows(x0=n2[s2],y0=Jpar2[s2,9],x1=n2[s2+1],y1=Jpar2[(s2+1),9], 
       col="#000000",lwd=2,length=0.05)
abline(tLL,0,col="#CC79A7", lty=2)
legend("topleft", legend = c("with initial value setup","random intial value"),
        col = c("#56B4E9","#000000"), lty=1, pch=c(1,6), bty = "n", cex = 1)

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#codon setup
num = 61
cdc = codon_call(num)
cod = cdc[[1]]
codonstrs = cdc[[2]]
syn = cdc[[3]]


trueP = read.table(inF, header=T, sep="")
tMat  = as.matrix(trueP)
tP    = tMat[1,]
#True parameters, unnormalized
Pi    = tP[1:4]
Sigma = tP[5:10]
Tau   = tP[11]
omega = tP[12]

#Generate 100 sequential tau/omega
cv = 1
a  = 1/(cv^2)
b  = 0.4/a
sv = matrix(0,100,6)
tv = c()
for (i in 1:100) {
  si  = rgamma(6, shape=a, scale=b)
  gtr = GTR(si, Pi)
  sv[i,] = si 
  tv[i]  = -sum(diag(gtr)*Pi)
}
range(tv)
plot(density(tv))

wv = seq(0.05,0.50,length.out=100)
pv = cbind(sv,tv)
pv = pv[order(tv),]
pg = expand.grid(pv[,7],wv)

#substituion model building
gtr  = GTR(Sigma,Pi)
mg94 = MG94(gtr,omega,cod,codonstrs,syn,num)

#pmat
Pi2 = sapply(seq(num), function(x){prod(Pi[match(cod[x, ], DNA_BASES)])})  
Pi2 = Pi2/sum(Pi2)

o   = outer(sqrt(Pi2),1/sqrt(Pi2))      
s94 = mg94*o                            #symmetric matrix
p94 = expm(s94)*t(o) 
P94 = p94*Pi2        
print(sum(P94))

#psudeo data
set.seed(8088)
dat = pseudo_data(1e+5,P94,num)

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

llv = c()
j   = 1
for(i in 1:nrow(pg)){
  llv[i] = LL(cod,codonstrs,syn,c(pv[j,1:6],pg[i,2]),f0,cf0,dat,num)
  if(i%%100 == 0){
    cat("i: ",i, "\n")
    j=j%%100  
  }
  j=j+1
}




#2-D contour plot[tau~omega] 
Jtmp1  = fromJSON(inF1)
Jpar1  = Jtmp1$p.intermed
Jtmp2  = fromJSON(inF2)
Jpar2  = Jtmp2$p.intermed

llm = matrix(llv, nrow=100)
contour(x=sort(tv),y=wv,z=llm,
        col="#000000",xlab="tau",ylab="omega",main="log-likelihood ")

points(Jpar1[,c(7,10)],col="#56B4E9")
points(Jpar2[,c(7,10)],col="#D55E00")
points(c(Tau,omega),col="#CC79A7")
