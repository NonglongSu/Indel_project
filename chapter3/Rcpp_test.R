

set.seed(1234)
library(tidyverse)
library(Rcpp)

data.frame(ID = sample(1:50000, 100000, replace = TRUE),
           permit = sample(c("student", "work", "refugee"), 100000, replace = TRUE)) %>%
  dplyr::distinct(ID, permit) %>%
  dplyr::arrange(ID) -> df

head(df)

#####################################
#setwd("~/Dropbox (ASU)/Indel_project/chapter3")

Rcpp::sourceCpp("../Script/chapter3/c-plus-plus-fun.cpp")

system.time(
  {
    df %>%
      dplyr::pull(permit) -> chr_permit
    df %>%
      dplyr::pull(ID) -> num_id
    
    combined_permits(chr_permit, num_id) -> df_cpp
  }
) 
head(df_cpp)

#########################################Test

Rcpp::sourceCpp("../Script/chapter3/c++fun.cpp")
system.time(
  {
    df %>%
      dplyr::pull(permit) -> chr_permit
    df %>%
      dplyr::pull(ID) -> num_id
    
    combined_results(chr_permit, num_id) -> df_cpp
  }
)



