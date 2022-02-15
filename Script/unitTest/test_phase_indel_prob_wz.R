library(testthat)
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(seqinr))
suppressPackageStartupMessages(library(stringi))
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(jsonlite))
suppressPackageStartupMessages(library(dfoptim))
suppressPackageStartupMessages(library(tidyverse))

source("../chapter3/phase_indel_prob_wz.R")
source("../sources/codon_call.R")

co.res = codon_call()
codons    = co.res[[1]]
codonstrs = co.res[[2]]
syn       = co.res[[3]]

#pseudo input
A1 = readBStringSet("sample/test_gap_sum.fa")


##expected results
res1 = matrix(0,6,3)
res1[4,2:3] = c(1,1)
colnames(res1) = c('zn','zs','long')
rownames(res1) = c('I0','I1','I2','D0','D1','D2')
res2 = list('ins'=c(),'del'=343,'l.ins'=c(),'l.del'=88)  
res3 = c(0,2)
res4 = c(0,11)




#setwd("~/Dropbox/Indel_project/Script/unitTest")
################################################

#test gap_sum function
test_that("test gap_sum()",{
  expect_equal(gap_sum(A1,syn,12)[[1]], res1)
  expect_equal(gap_sum(A1,syn,12)[[2]], res2)
  expect_equal(gap_sum(A1,syn,12)[[3]], res3)
  expect_equal(gap_sum(A1,syn,12)[[4]], res4)
})


