library(expm)
library(stats)
library(Biostrings)

###################function setup
# Set up 4*4 GTR matrix
GTR = function(si, pai){
  r1 = matrix(0, 4, 4)
  r1[lower.tri(r1)] = si
  r  = r1 + t(r1)
  r  = t(r*pai)
  diag(r) = -rowSums(r)
  return(r)
}

# Construct a 64*64 MG94 matrix
MG94 = function(g, oga, cd, m){
  R  = matrix(0, m, m)
  for (i in 1:m) {
    for (j in 1:m) {
      if(i == j){
        R[i, j] = 0
      }else if(sum(cd[i, ] != cd[j, ]) > 1){#more than 1 nucleotide changes
        R[i, j] = 0
      }else{
        if(codonstrs[j] %in% syn[[codonstrs[i]]]){#syn_subs
          w = 1
        }else{#nonsyn_subs
          w = oga
        }
        pos = which(cd[i, ] != cd[j, ])
        x   = which(DNA_BASES == cd[i, pos])
        y   = which(DNA_BASES == cd[j, pos])
        R[i, j] = w * g[x, y]
      }
    }
  }
  diag(R) = -rowSums(R)

  return(R)
}


set.seed(8088)
#######################PART I: parameter space set up
##pi
##2(A+C)=1; highC-G or highA-T case.
nd  = 100
PiA = runif(nd,0.0,0.50)
PiT = PiA
PiC = 0.50 - PiA
PiG = PiC
Pi.all = list(PiA, PiC, PiG, PiT)
Pi.all = sapply(1:nd, function(x){sapply(Pi.all, "[[", x)})
print(dim(Pi.all))
print(colSums(Pi.all)) #1

##omega
##from paper
wv = runif(nd,0.02,0.5) #sim

##six sigmas
##set up the mean rate first; then b =1/a ; 1-sigma rule
##cv = 1; a = 1/(cv^2)

##distance Tau (total mutations load--only determined by neutral gtr matrix)
#Tau = runif(nd,0.05,5.5)  #(put 80% of the simulations in the area of saturated third-positions,seems too unreal!!)
#we setup the a*b=0.4
cv = 1
a  = 1/(cv^2)
b  = 0.4/a
Sigmas = matrix(0,6,100)
tv = c()
for (i in 1:100) {
  pai = Pi.all[,i]
  si  = rgamma(6, shape=a, scale=b)
  gtr = GTR(si, pai)
  Sigmas[,i] = si
  tv[i]  = -sum(diag(gtr)*pai)
}

mean(tv)
hist(tv, prob = TRUE, xlim = c(0,2))
stv = sort(tv)
lines(stv,dgamma(stv,shape=a,scale=b),col='magenta',lwd=3)
plot(density(tv), xlab='tau', ylab='density')


##normalize sigma.
nSigmas = sapply(1:100, function(x){Sigmas[,x]/tv[x]})
print(nSigmas)

##output the 100 parameters
tPar = matrix(0,100,12)
for (i in 1:nd) {
  tPar[i,1:4]  = Pi.all[,i]
  tPar[i,5:10] = Sigmas[,i]
  tPar[i,11]   = tv[i]
  tPar[i,12]   = wv[i]
}
colnames(tPar)=c('A','C','G','T','s1','s2','s3','s4','s5','s6','tau','omega')
write.table(tPar,"Results/truePar_100.txt",quote=F,sep="\t",
            row.names=F)
