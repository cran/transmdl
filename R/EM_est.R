#'@title  Estimate parameters and hazard function via EM algorithm.
#'
#'@description   Estimate the vector of parameters for baseline covariates
#'  \eqn{\beta} and baseline cumulative hazard function \eqn{\Lambda(\cdot)}
#'  using the expectation-maximization algorithm. \eqn{\Lambda(t)} is estimated
#'  as a step function with jumps only at the observed failure times. Typically,
#'  it would only be used in a call to \code{trans.m} or \code{Simu}.
#'
#'@return a list containing
#'\tabular{lccl}{ \code{beta_new} \tab\tab\tab
#'  estimator of \eqn{\beta} \cr \code{Lamb_Y} \tab\tab\tab estimator of
#'  \eqn{\Lambda(Y)} \cr \code{lamb_Y} \tab\tab\tab estimator of
#'  \eqn{\lambda(Y)} \cr \code{lamb_Ydot} \tab\tab\tab estimator of
#'  \eqn{\lambda(Y')} \cr \code{Y_eq_Yhat} \tab\tab\tab a matrix used in
#'  \code{trans.m} and \code{Simu} \cr \code{Y_geq_Yhat} \tab\tab\tab a matrix
#'  used in \code{trans.m} and \code{Simu} \cr }
#'
#'
#' @examples
#'   gen_data = generate_data(200, 1, 0.5, c(-0.5, 1))
#'   delta = gen_data$delta
#'   Y = gen_data$Y
#'   X = gen_data$X
#'   EM_est(Y, X, delta, alpha = 1)$beta_new - c(-0.5, 1)
#'
#'
#'@references Abramowitz, M., and Stegun, I.A. (1972). Handbook of Mathematical
#'  Functions (9th ed.). Dover Publications, New York.
#'
#'  Evans, M. and Swartz, T. (2000). Approximating Integrals via Monte Carlo and
#'  Deterministic Methods. Oxford University Press.
#'
#'  Liu, Q. and Pierce, D.A. (1994). A note on Gauss-Hermite quadrature.
#'  Biometrika 81: 624-629.
#'
#'
#'
#'@param Y observed event times
#'@param X design matrix
#'@param delta censoring indicator. If \eqn{Y_i} is censored, \code{delta}=0. If
#'  not, \code{delta}=1.
#'@param alpha parameter in transformation function
#'@param Q number of nodes and weights in Gaussian quadrature. Defaults to 60.
#'@param EM_itmax maximum iteration of EM algorithm. Defaults to 250.
#'
#'
#'@export
EM_est = function(Y, X, delta, alpha, Q = 60, EM_itmax = 250) {


  X.size = ncol(X)
  gq = statmod::gauss.quad(Q, kind = "laguerre")
  nodes = gq$nodes
  whts = gq$weights

  Y_hat = sort(unique(Y[delta == 1]))
  l = length(Y_hat)

  Y_dot = sort(unique(Y))
  n = length(Y)
  m = length(Y_dot)

  lambda_Ydot = matrix(c(1/l), ncol = m)

  # Y1 Y2 ... Yn Y1 Y2 ... Yn l*n
  Y_mat = matrix(rep(Y, l), nrow = l, byrow = TRUE)

  # Y1 Y2 ... Yn Y1 Y2 ... Yn m*n
  Y_mat_mn = matrix(rep(Y, m), nrow = m, byrow = TRUE)

  # Yhat 1 ... Yhat 1 Yhat 2 ... Yhat 2 ... Yhat l ... Yhat l l*n
  Y_hat_mat = matrix(rep(Y_hat, n), nrow = l, byrow = FALSE)

  # Ydot 1 ... Ydot 1 Ydot 2 ... Ydot 2 ... Ydot m ... Ydot m m*n
  Y_dot_mat = matrix(rep(Y_dot, n), nrow = m, byrow = FALSE)

  Y_eq_Yhat = (Y_mat == Y_hat_mat)
  Y_geq_Yhat = (Y_mat >= Y_hat_mat)

  Lambda_Y_initial = colSums(t(lambda_Ydot[1, ] * t(Y_dot <= Y_mat_mn)))

  lambda_Y_initial = colSums(t(lambda_Ydot[1, ] * t(Y_dot == Y_mat_mn)))

  Lambda_Y = matrix(Lambda_Y_initial, ncol = n)

  lambda_Y = matrix(lambda_Y_initial, ncol = n)

  beta = matrix(0, ncol = X.size)

  m_Y_k = matrix(Y, ncol = n, nrow = n, byrow = TRUE)
  m_Y_i = matrix(Y, ncol = n, nrow = n, byrow = FALSE)
  m_Y_kgeqi = (m_Y_k >= m_Y_i)

  dgamma_Q = stats::dgamma(nodes, shape = 1/alpha, scale = alpha)

  t = 1
  for (EM_repeat in 1:EM_itmax) {


    beta_X = c(X %*% beta[t, ])
    exp_beX = exp(beta_X)
    E_n = rep(0, each = n)
    E_d = rep(0, each = n)

    lam_expbeX_m = matrix(rep(lambda_Y[t, ] * exp_beX, Q), nrow = Q, byrow = TRUE)
    Lam_expbeX_m = matrix(rep(Lambda_Y[t, ] * exp_beX, Q), nrow = Q, byrow = TRUE)
    E_de_m = whts * t(t(nodes * lam_expbeX_m)^delta) * exp(-nodes * Lam_expbeX_m + matrix(rep(log(dgamma_Q) +
                                                                                                nodes, n), nrow = Q))
    E_nu_m = nodes * E_de_m
    E = colSums(E_nu_m)/colSums(E_de_m)


    lambda_Ydot_new = c()

    delta_k = rowSums(t(delta * t(Y_mat_mn == Y_dot_mat)))
    de_lam = rowSums(t(E * exp(beta_X) * t(Y_mat_mn >= Y_dot_mat)))
    lambda_Ydot_new = delta_k/de_lam

    lambda_Ydot = rbind(lambda_Ydot, lambda_Ydot_new)


    Lambda_Y_new = colSums(lambda_Ydot_new * (Y_dot <= Y_mat_mn))

    lambda_Y_new = colSums(lambda_Ydot_new * (Y_dot == Y_mat_mn))

    Lambda_Y = rbind(Lambda_Y, Lambda_Y_new)
    lambda_Y = rbind(lambda_Y, lambda_Y_new)



    beta_X = c(X %*% beta[t, ])

    lam_expbeX_m = matrix(rep(lambda_Y[t + 1, ] * exp_beX, Q), nrow = Q, byrow = TRUE)
    Lam_expbeX_m = matrix(rep(Lambda_Y[t + 1, ] * exp_beX, Q), nrow = Q, byrow = TRUE)
    E_de_m = whts * t(t(nodes * lam_expbeX_m)^delta) * exp(-nodes * Lam_expbeX_m + matrix(rep(log(dgamma_Q) +
                                                                                                nodes, n), nrow = Q))
    E_nu_m = nodes * E_de_m
    E = colSums(E_nu_m)/colSums(E_de_m)



    E_exp = E * exp(beta_X)
    I_E_exp = m_Y_kgeqi * matrix(E_exp, nrow = n, ncol = n, byrow = TRUE)
    IEexp_rowsum = rowSums(I_E_exp)
    f_brace = X - (I_E_exp %*% X)/rowSums(I_E_exp)
    f = colSums(delta * f_brace)

    X.size = ncol(X)
    XXT.lay = matrix(, nrow = n, ncol = (X.size * (X.size + 1))/2)
    aXaXT.lay = matrix(, nrow = n, ncol = (X.size * (X.size + 1))/2)
    aX.m = I_E_exp %*% X
    for (i in 1:n) {
      XXT.c = c()
      XXT.m = X[i, ] %*% t(X[i, ])
      aXaXT.c = c()
      aXaXT.m = aX.m[i, ] %*% t(aX.m[i, ])
      for (k in X.size:1) {
        XXT.c = c(XXT.c, XXT.m[1 + X.size - k, X.size - k + c(1:k)])
        aXaXT.c = c(aXaXT.c, aXaXT.m[1 + X.size - k, X.size - k + c(1:k)])
      }
      XXT.lay[i, ] = XXT.c
      aXaXT.lay[i, ] = aXaXT.c
    }

    # fdot.left.lay
    fd.lay = colSums(-delta * (((I_E_exp %*% XXT.lay)/IEexp_rowsum) - (aXaXT.lay/IEexp_rowsum^2)))

    fd = matrix(, nrow = X.size, ncol = X.size)

    v2m.ind = 0
    for (k in X.size:1) {
      fd[X.size + 1 - k, X.size - k + c(1:k)] = fd.lay[v2m.ind + c(1:k)]
      v2m.ind = v2m.ind + k
    }
    fd[lower.tri(fd)] = fd[upper.tri(fd)]

    fdot = fd


    beta_new = beta[t, ] - f %*% solve(fdot)

    beta = rbind(beta, beta_new)

    beta[t + 1, ]

    if ((sum(abs(beta[t + 1, ] - beta[t, ])) <= 10^(-5)) & (max(abs(lambda_Ydot[t + 1, ] -
                                                                    lambda_Ydot[t, ])) <= 10^(-5))) {
      break
    }

    if (EM_repeat == (EM_itmax)) {
      cat("EM iteration >", EM_repeat)
    }

    t = t + 1
  }

  return(list(beta_new = beta_new, Lamb_Y = Lambda_Y_new, lamb_Y = lambda_Y_new, lamb_Ydot = lambda_Ydot_new,
              Y_eq_Yhat = Y_eq_Yhat, Y_geq_Yhat = Y_geq_Yhat))
}
