suppressPackageStartupMessages(library(testthat))

#Test each function in sw_gap.R
testthat::context("testing process start:")

#setup the window and wall size
Window <<- 6
Wall   <<- 12

##Pseduo data 
seq1 = "AAT===AAACAAAGAATGCTTACTGT---ATAAGGCTTACTGTTCTAGCG===ATCACCGCG===TCATGTCTAGTTATGAACGGC------GGTTTAACATTGAATAGCAAGGCACTTCCA---TAATAGGGCCGTC===GTAATTGTCTAATATAG------ATAGTA==="
seq2 = "TAA------AA===AATTTGATGCTACATTGGATGAGTCTACTTCGAGCGCGCCGCATCGATTGCAAGAGCAGTGTTGCCT===AAGAGCCGTTAGATGCGTCGTTG---ATCGCGTCCGATAATTCGGGAGTTGTGC===CCCAATATTTAATATGATGA===TAGCTATAA"

#global vars
lwall = Wall 
rwall = Wall + 2

##Predicted outcomes
#index
res1 = c(27, 123)
res2 = 108

#alignemnt score
res3 = 25
res4 = 25

#left/right slide mode
seq1.l = "AAT===AAACAAAGAATGCTTACTG---TATAAGGCTTACTGTTCTAGCG===ATCACCGCG===TCATGTCTAGTTATGAACGGC------GGTTTAACATTGAATAGCAAGGCACTTCCA---TAATAGGGCCGTC===GTAATTGTCTAATATAG------ATAGTA==="
seq2.l = "TAA------AA===AATTTGATGCTACATTGGATGAGTCTACTTCGAGCGCGCCGCATCGATTGCAAGAGCAGTGTTGCCT===AAGAGCCGTTAGATGCGTCGTT---GATCGCGTCCGATAATTCGGGAGTTGTGC===CCCAATATTTAATATGATGA===TAGCTATAA"
  
seq1.r = "AAT===AAACAAAGAATGCTTACTGTA---TAAGGCTTACTGTTCTAGCG===ATCACCGCG===TCATGTCTAGTTATGAACGGC------GGTTTAACATTGAATAGCAAGGCACTTCCA---TAATAGGGCCGTC===GTAATTGTCTAATATAG------ATAGTA==="
seq2.r = "TAA------AA===AATTTGATGCTACATTGGATGAGTCTACTTCGAGCGCGCCGCATCGATTGCAAGAGCAGTGTTGCCT===AAGAGCCGTTAGATGCGTCGTTGA---TCGCGTCCGATAATTCGGGAGTTGTGC===CCCAATATTTAATATGATGA===TAGCTATAA"

#Best alignment 
#unchanged
best_aligned.m  = "AAT===AAACAAAGAATGCTTACTGT+++ATAAGGCTTACTGTTCTAGCG===ATCACCGCG===TCATGTCTAGTTATGAACGGC------GGTTTAACATTGAATAGCAAGGCACTTCCA---TAATAGGGCCGTC===GTAATTGTCTAATATAG------ATAGTA==="
#changed
best_aligned.m2 = "AAT===AAACAAAGAATGCTTACTGT---ATAAGGCTTACTGTTCTAGCG===ATCACCGCG===TCATGTCTAGTTATGAACGGC------GGTTTAACATTGAATAGCAAGGCACTTCCATAATAG---GGCCGTC===GTAATTGTCTAATATAG------ATAGTA==="
best_aligned.r  = "TAA------AA===AATTTGATGCTACATTGGATGAGTCTACTTCGAGCGCGCCGCATCGATTGCAAGAGCAGTGTTGCCT===AAGAGCCGTTAGATGCGTCGTTGATC---GCGTCCGATAATTCGGGAGTTGTGC===CCCAATATTTAATATGATGA===TAGCTATAA"


#Test Align func. 
wid.1     = substr(seq1, start =res1[1]- lwall, stop = res1[1] + rwall)
wid_ref.1 = substr(seq2, start =res1[1]- lwall, stop = res1[1] + rwall)
wid.2     = substr(seq2, start =res2[1]- lwall, stop = res2[1] + rwall)
wid_ref.2 = substr(seq1, start =res2[1]- lwall, stop = res2[1] + rwall)


##########################################################################
source("Script/sw_gap.R")


test_that("Test Align() func", {
  expect_equal(Align(wid.1, wid_ref.1), res3)
  expect_equal(Align(wid.2, wid_ref.2), res4)
})

#Test left_slide func.
test_that("Test left_slide() func",{
  expect_equal(left_slide(str_convert(seq1),res1[1],3), seq1.l)
  expect_equal(left_slide(str_convert(seq2),res2[1],3), seq2.l)
})

#Test right_slide func.
test_that("Test right_slide() func",{
  expect_equal(right_slide(str_convert(seq1),res1[1],3), seq1.r)
  expect_equal(right_slide(str_convert(seq2),res2[1],3), seq2.r)
})

#Test Merge func.
test_that("Test Merge()",{
   expect_equal(Merge(res1[1], 3,seq1, seq2), best_aligned.m)
   expect_equal(Merge(res1[2], 3,seq1, seq2), best_aligned.m2)
   expect_equal(Merge(res2[1], 3,seq2, seq1), best_aligned.r)
})




