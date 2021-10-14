// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <Rcpp.h>
using namespace Eigen;
using namespace std;
using namespace Rcpp;

//[[Rcpp::export]]
Rcpp::List ESSEB2 (int l, int n, int Q, int size, Eigen::VectorXd nodes, Eigen::VectorXd whts, Eigen::VectorXd lambda_Y_new, Eigen::VectorXd Lambda_Y_new,
                   Eigen::VectorXd beta_X, Eigen::VectorXd delta_i, Eigen::VectorXd dgamma_Q, Eigen::VectorXd Lambda_b_nozero,
                   Eigen::MatrixXd Y_eq_Yhat, Eigen::MatrixXd Y_geq_Yhat, Eigen::MatrixXd X)
{
  MatrixXd ESS;
  ESS.setZero((l+size),(l+size));
  MatrixXd EB;
  EB.setZero((l+size),(l+size));
  for(int i = 0; i<n; ++i)
  {
    double ESS_de_ig = 0;
    MatrixXd ESS_nu_ig;
    MatrixXd EB_nu_ig;
    ESS_nu_ig.setZero((l+size),(l+size));
    EB_nu_ig.setZero((l+size),(l+size));
    for(int g = 0; g<Q; ++g)
    {
      double xi_i = nodes(g);
      double de_each;
      de_each = whts(g) *  pow((xi_i* lambda_Y_new(i)* exp(beta_X(i))),delta_i(i))*
        exp((-1)*xi_i* Lambda_Y_new(i)* exp(beta_X(i)))*
        dgamma_Q(g) * exp(xi_i);
      ESS_de_ig = ESS_de_ig + de_each;

      VectorXd ESS_nu_v;
      ESS_nu_v.setZero(l+size);
      ESS_nu_v.segment(0,size) = (delta_i(i)-xi_i*Lambda_Y_new(i)*exp(beta_X(i)))*X.row(i);
      ESS_nu_v.segment(size,l) = delta_i(i)*(Y_eq_Yhat.col(i).array())/Lambda_b_nozero(i)-
        xi_i*exp(beta_X(i))*(Y_geq_Yhat.col(i).array());
      ESS_nu_ig = ESS_nu_ig + de_each * (ESS_nu_v*ESS_nu_v.transpose());

      MatrixXd EB_nu_m;
      EB_nu_m.setZero((l+size),(l+size));
      EB_nu_m.block(0,0,size,size) = (X.row(i).transpose() * X.row(i)) * xi_i * Lambda_Y_new(i) * exp(beta_X(i));
      EB_nu_m.block(0,size,size,l) = xi_i * exp(beta_X(i)) * ( Y_geq_Yhat.col(i) * X.row(i)).transpose();
      EB_nu_m.diagonal().segment(size,l) = delta_i(i) * (Y_eq_Yhat.col(i))/pow(Lambda_b_nozero[i],2);
      EB_nu_ig = EB_nu_ig + de_each * EB_nu_m;
    }
    ESS = ESS + ESS_nu_ig / ESS_de_ig;
    EB = EB + EB_nu_ig / ESS_de_ig;
  }
  EB.block(size,0,l,size) = EB.block(0,size,size,l).transpose();

  Rcpp::List list_ESSEB = Rcpp::List::create(Rcpp::Named("ESS", ESS),
                                             Rcpp::Named("EB", EB));
  return list_ESSEB;
}


