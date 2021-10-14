// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <Rcpp.h>
using namespace Eigen;
using namespace std;
using namespace Rcpp;

//[[Rcpp::export]]
Eigen::MatrixXd ddloglik_transmdl (int n, Eigen::VectorXd delta, Eigen::MatrixXd Z, Eigen::VectorXd beta, Eigen::VectorXd E )
{
  int p = beta.size();
  MatrixXd L;
  L.setZero(p, p);
  for (int i = 0; i < n; i++)
  {
    double temp = 0;
    VectorXd temp1;
    temp1.setZero(p);
    MatrixXd temp2;
    temp2.setZero(p, p);
    for (int j = i; j < n; j++)
    {
      temp = temp + exp(Z.row(j) * beta) * E(j);
      temp1 = temp1 + Z.row(j).transpose() * exp(Z.row(j) * beta) * E(j);
      temp2 = temp2 + Z.row(j).transpose() * Z.row(j) * exp(Z.row(j) * beta) * E(j);
    }
    L = L + delta(i) * (temp2 / temp - (temp1 * temp1.transpose()) / (temp*temp) );
  }
  return L;
}


