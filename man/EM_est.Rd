% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/EM_est.R
\name{EM_est}
\alias{EM_est}
\title{Estimate parameters and hazard function via EM algorithm.}
\usage{
EM_est(Y, X, delta, alpha, Q = 60, EM_itmax = 250)
}
\arguments{
\item{Y}{observed event times}

\item{X}{design matrix}

\item{delta}{censoring indicator. If \eqn{Y_i} is censored, \code{delta}=0. If
not, \code{delta}=1.}

\item{alpha}{parameter in transformation function}

\item{Q}{number of nodes and weights in Gaussian quadrature. Defaults to 60.}

\item{EM_itmax}{maximum iteration of EM algorithm. Defaults to 250.}
}
\value{
a list containing
\tabular{lccl}{ \code{beta_new} \tab\tab\tab
 estimator of \eqn{\beta} \cr \code{Lamb_Y} \tab\tab\tab estimator of
 \eqn{\Lambda(Y)} \cr \code{lamb_Y} \tab\tab\tab estimator of
 \eqn{\lambda(Y)} \cr \code{lamb_Ydot} \tab\tab\tab estimator of
 \eqn{\lambda(Y')} \cr \code{Y_eq_Yhat} \tab\tab\tab a matrix used in
 \code{trans.m} and \code{Simu} \cr \code{Y_geq_Yhat} \tab\tab\tab a matrix
 used in \code{trans.m} and \code{Simu} \cr }
}
\description{
Estimate the vector of parameters for baseline covariates
 \eqn{\beta} and baseline cumulative hazard function \eqn{\Lambda(\cdot)}
 using the expectation-maximization algorithm. \eqn{\Lambda(t)} is estimated
 as a step function with jumps only at the observed failure times. Typically,
 it would only be used in a call to \code{trans.m} or \code{Simu}.
}
\examples{
  gen_data = generate_data(200, 1, 0.5, c(-0.5, 1))
  delta = gen_data$delta
  Y = gen_data$Y
  X = gen_data$X
  EM_est(Y, X, delta, alpha = 1)$beta_new - c(-0.5, 1)


}
\references{
Abramowitz, M., and Stegun, I.A. (1972). Handbook of Mathematical
 Functions (9th ed.). Dover Publications, New York.

 Evans, M. and Swartz, T. (2000). Approximating Integrals via Monte Carlo and
 Deterministic Methods. Oxford University Press.

 Liu, Q. and Pierce, D.A. (1994). A note on Gauss-Hermite quadrature.
 Biometrika 81: 624-629.
}
