#include <Rcpp.h>
#include <stdlib.h>

using namespace std;
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
double bpp(int m, double alphat, double betat, double alphac, double betac, double eta, double delta){

  double prob = 0.0;
  double exps = std::numeric_limits<double>::epsilon();

  for (int t = 0; t < m+1; t++)
  {
    for (int c = 0; c < m+1; c++)
    {
      double pxy=0;

      for (double i = 0+0.001; i < 1.0; i+=0.001)
      {
        pxy += R::pbeta(i-delta,alphac+c,betac+m-c,1,0)*R::dbeta(i,alphat+t,betat+m-t,0)*0.001;
      }

      if (pxy > eta || std::abs(pxy - eta) < exps){
        prob += exp(
                    R::lchoose(m,t)+R::lbeta(alphat+t,betat+m-t)-R::lbeta(alphat,betat)
                   +R::lchoose(m,c)+R::lbeta(alphac+c,betac+m-c)-R::lbeta(alphac,betac)
                   );

      }
    }

  }
  return prob;
}
