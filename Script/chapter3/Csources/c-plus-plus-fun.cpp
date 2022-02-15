#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::List combined_permits(Rcpp::CharacterVector permit,
                            Rcpp::NumericVector id) {

  int size = id.length();
  Rcpp::String x("_");

  for(int i = 0; i < size - 1; i++)

  {

    if(id[i] == id[i + 1])

    {
      permit[i + 1] = (permit[i] += x) += permit[i + 1];

    }

  }

  List df = List::create(_["comb"] = permit);
  df.attr("class") = "data.frame";
  df.attr("row.names") = seq(1, size);
  return df;

}
