// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <Rcpp.h>
using namespace Eigen;
using namespace std;
using namespace Rcpp;

//[[Rcpp::export]]
Eigen::VectorXd shoot (Eigen::MatrixXd X, Eigen::VectorXd y, Eigen::VectorXd beta_ini_s, double lamb_s )
{

  VectorXd beta_old_s = beta_ini_s;
  VectorXd beta_new_s = beta_old_s;

  for (int ii = 0; ii < 100; ii++)
  {
    VectorXd beta_tmp_s;
    VectorXd beta_tmp2_s;

    beta_tmp2_s.setZero(beta_ini_s.size());

    for (int j = 0; j < beta_ini_s.size(); j++)
    {
      beta_tmp_s = beta_new_s;
      beta_tmp_s(j) = 0;
      //an 1*1 matrix behind 2* is incorrect.
      double S0_sht;
      S0_sht = 2 * ((X.col(j).transpose() * X * beta_tmp_s)(0, 0) - (X.col(j).transpose() * y)(0, 0));

      if (S0_sht > lamb_s / abs(beta_new_s(j)))
      {
        beta_tmp2_s(j) = (lamb_s - S0_sht) / (2 * X.col(j).transpose() * X.col(j));
        //cout << 1;
      }
      else if (S0_sht < -lamb_s / abs(beta_new_s(j)))
      {
        beta_tmp2_s(j) = (-lamb_s - S0_sht) / (2 * X.col(j).transpose() * X.col(j));
        //cout << 2;
      }
      else
      {
        beta_tmp2_s(j) = 0;
        //cout << 3;
      }
    }

    beta_new_s = beta_tmp2_s;
    if ((beta_new_s - beta_old_s).cwiseAbs().maxCoeff() < 0.00001)
    {
      break;
    }
    beta_old_s = beta_new_s;
  }
  //Rcpp::List ll = Rcpp::List::create(Rcpp::Named("X", X));
  return beta_new_s;
}


