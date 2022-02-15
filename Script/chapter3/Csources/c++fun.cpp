#include <Rcpp.h>
#include <Csources/codon_call.h>
#include <Csources/gtr.h>
#include <Csources/mg94.h>
#include <Csources/init_f.h>
#include <Csources/LL.h>
#include <Csources/phase_indel_prob2.h>

using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::List combined_results(Rcpp::NumericVector ngroups,
                            Rcpp::NumericVector par1,
                            Rcpp::NumericVector par2,
                            Rcpp::NumericVector par3,
                            Rcpp::NumericVector par4,
                            Rcpp::CharacterVector aseq,
                            Rcpp::CharacterVector bseq) {



  for(int i = 0; i < size - 1; i++){
    Rcpp::List res = phase_indel_prob2::ziqi_prob(aseq[i],bseq[i],par1,par2,par3,par4);
    Rcpp::List x1  = res[[0]];
    Rcpp::List x2  = res[[1]];
    Rcpp::List x3  = res[[2]];
    Rcpp::NumericVector x4 = res[[3]];
    Rcpp::List x5 = res[[4]];
  }



  Rcpp::List df;
  df["Gm"]     = x1;
  df["Mm"]     = x2;
  df["codArr"] = x3;
  df["llz"]    = x4;
  df["gl"]     = x5;


  return df;

}
