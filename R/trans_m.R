#'@title  Regression Analysis of Right-censored Data using Semiparametric
#'  Transformation Models.
#'
#'@description   This function is used to conduct the regression analysis of
#'  right-censored data using semiparametric transformation models. It
#'  calculates the estimators, standard errors and p values. A plot of estimated
#'  baseline cumulative hazard function and confidence intervals can be
#'  produced.
#'
#'
#'@details If \eqn{\alpha} is unknown, we firse set \eqn{\alpha=}\code{alpha}.
#'  Then, for each \eqn{\alpha}, we estimate the parameters and record the value
#'  of observed log-likelihood function. The \eqn{\alpha} that maximizes the
#'  observed log-likelihood function and the corresponding \eqn{\hat\beta} and
#'  \eqn{\hat\Lambda(\cdot)} are chosen as the best estimators. Nonparametric
#'  maximum likelihood estimators are developed for the regression parameters
#'  and cumulative intensity functions of these models based on censored data.
#'
#'@return a list containing \tabular{lccl}{ \code{beta.est } \tab\tab\tab
#'  estimators of \eqn{\beta} \cr \code{SE.beta } \tab\tab\tab standard errors
#'  of the estimated \eqn{\beta} \cr \code{SE.Ydot} \tab\tab\tab standard errors
#'  of the estimated \eqn{\Lambda(Y')} \cr \code{Ydot} \tab\tab\tab vector of
#'  sorted event times with duplicate values removed   \cr \code{Lamb.est}
#'  \tab\tab\tab estimated baseline cumulative hazard  \cr \code{lamb.est}
#'  \tab\tab\tab estimated jump sizes of baseline cumulative hazard function \cr
#'  \code{choose.alpha } \tab\tab\tab the chosen \eqn{\alpha} \cr
#'  \code{Lamb.upper} \tab\tab\tab upper confidence limits for the estimated
#'  baseline cumulative hazard function \cr \code{Lamb.lower} \tab\tab\tab lower
#'  confidence limits for the estimated baseline cumulative hazard function \cr
#'  \code{p.beta} \tab\tab\tab P values of estimated \eqn{\beta} \cr
#'  \code{p.Lamb} \tab\tab\tab P values of estimated baseline cumulative hazard
#'  \cr p.beta }
#'
#'
#' @examples
#'  gen_data = generate_data(200, 1, 0.5, c(-0.5,1))
#'  delta = gen_data$delta
#'  Y = gen_data$Y
#'  X = gen_data$X
#'  res.trans = trans_m(X, delta,  Y, plot.Lamb = TRUE, show_res = FALSE)
#'
#'
#'@references Cheng, S.C., Wei, L.J., and Ying, Z. (1995). Analysis of
#'  transformation models with censored data. Biometrika 82, 835-845.
#'
#'  Zeng, D. and Lin, D.Y. (2007). Maximum likelihood estimation in
#'  semiparametric regression models with censored data. J. R. Statist. Soc. B
#'  69, 507-564.
#'
#'  Abramowitz, M., and Stegun, I.A. (1972). Handbook of Mathematical Functions
#'  (9th ed.). Dover Publications, New York.
#'
#'  Evans, M. and Swartz, T. (2000). Approximating Integrals via Monte Carlo and
#'  Deterministic Methods. Oxford University Press.
#'
#'  Liu, Q. and Pierce, D.A. (1994). A note on Gauss-Hermite quadrature.
#'  Biometrika 81, 624-629.
#'
#'  Louis, T. (1982). Finding the Observed Information Matrix when Using the EM
#'  Algorithm. Journal of the Royal Statistical Society. Series B
#'  (Methodological), 44(2), 226-233.
#'
#'
#'
#'@param Y observed event times
#'@param X design matrix
#'@param delta censoring indicator. If \eqn{Y_i} is censored, \code{delta}=0. If
#'  not, \code{delta}=1.
#'@param alpha parameter in transformation function. Generally, \eqn{\alpha} can
#'  not be observed in medical applications. In that situation, \code{alpha}
#'  indicates the scale of choosing \eqn{\alpha}. The default is
#'  \eqn{(0.1,0.2,...,1.1)}. If \eqn{\alpha} is known, \code{alpha} indicates
#'  the true value of \eqn{\alpha}.
#'@param plot.Lamb If TRUE, plot the estimated baseline cumulative hazard
#'  function and confidence intervals. The default is TRUE.
#'@param EM_itmax maximum iteration of EM algorithm. Defaults to 250.
#'@param trsmodel logical value indicating whether to implement transformation
#'  models. The default is TRUE.
#'@param show_res show results after \code{trans_m}.
#'
#'@seealso \code{\link{EM_est}}
#'
#'@export
trans_m = function( X, delta, Y, plot.Lamb = TRUE,
                    alpha = seq(0.1,1.1,by=0.1), trsmodel = TRUE, EM_itmax = 250, show_res = TRUE)
{

  X.size = ncol(X)

  #alpha.true = 1

  result_SE = c() #beta,lambda


  alpha.c = alpha

  length_alpha = length(alpha.c)
  result_alpha = c()
  result_lik = matrix(0,nrow = 1, ncol = length_alpha)

  H.fun = function(x,alpha)
  {
    if(alpha == 0)
    {
      x
    }else if(alpha > 0)
    {
      log(1+alpha * x) / alpha
    }else
    {
      cat("H:wrong alpha < 0")
    }
  }

  Hdot.fun = function(x,alpha)
  {
    if(alpha == 0)
    {
      1
    }else if(alpha > 0)
    {
      1 / (1+alpha*x)
    }else
    {
      cat("Hdot:wrong alpha < 0")
    }
  }

  Q = 60
  gq = statmod::gauss.quad(Q, kind="laguerre")
  nodes = gq$nodes
  whts  = gq$weights

  # plot_x = seq(0.1,2.5,by=0.1)
  # n_x = length(plot_x)
  # plot_Lam_Ydot = c()

  u_975 = stats::qnorm(0.975)

  I_1 = which(Y ==  1- min(1-Y[Y<=1]))[1]
  I_2 = which(Y ==  2- min(2-Y[Y<=2]))[1]

  Y_hat = sort(unique(Y[delta==1]))
  l = length(Y_hat)

  Y_dot = sort(unique(Y))
  n = length(Y)
  m = length(Y_dot)

  #Y1 Y2 ... Yn
  #Y1 Y2 ... Yn  l*n
  Y_mat = matrix(rep(Y,l), nrow = l, byrow = TRUE)

  #Y1 Y2 ... Yn
  #Y1 Y2 ... Yn  m*n
  Y_mat_mn = matrix(rep(Y,m), nrow = m, byrow = TRUE)

  #Yhat 1 ... Yhat 1
  #Yhat 2 ... Yhat 2
  #...
  #Yhat l ... Yhat l   l*n
  Y_hat_mat = matrix(rep(Y_hat,n), nrow = l, byrow = FALSE)

  #Ydot 1 ... Ydot 1
  #Ydot 2 ... Ydot 2
  #...
  #Ydot m ... Ydot m   m*n
  Y_dot_mat = matrix(rep(Y_dot,n), nrow = m, byrow = FALSE)

  #SE
  Y_eq_Yhat = (Y_mat == Y_hat_mat)
  Y_geq_Yhat = (Y_mat >= Y_hat_mat)

  # beta^t+1
  m_Y_k = matrix(Y, ncol = n, nrow = n, byrow = TRUE)
  m_Y_i = matrix(Y, ncol = n, nrow = n, byrow = FALSE)
  m_Y_kgeqi = (m_Y_k>=m_Y_i)


  if((min(alpha) > 0) & (trsmodel == TRUE))
  {
    #print(1)
    lik.fun = c()
    beta_best = c()
    LambdaY_best = c()
    lambdaY_best = c()
    lambdaYd_best = c()
    for(al_i in 1:length(alpha.c))
    {
      alpha.i = alpha.c[al_i]
      # for(al_i in 1)
      # {
      #   alpha.i = 0.06
      EM.res = EM_est(Y = Y, X = X, delta = delta, alpha = alpha.i, Q = Q, EM_itmax = EM_itmax )
      beta_new = EM.res$beta_new
      Lambda_Y_new = EM.res$Lamb_Y
      lambda_Y_new = EM.res$lamb_Y
      lambda_Ydot_new = EM.res$lamb_Ydot

      beX = c(X %*% c(beta_new))
      Lam.ExpbeX = Lambda_Y_new * exp(c(X %*% c(beta_new)))
      lam.ExpbeX = lambda_Y_new * exp(c(X %*% c(beta_new)))
      #log(prod((lam.ExpbeX*(Hdot.fun(Lam.ExpbeX,alpha.i)))^delta * exp(-(H.fun(Lam.ExpbeX,alpha.i))))  )
      aa = delta*(log(lambda_Y_new) + beX + log(Hdot.fun(Lam.ExpbeX,alpha.i)))
      aa[is.na(aa)] = 0
      lik.new =  sum(aa - H.fun(Lam.ExpbeX,alpha.i))
      lik.fun = c( lik.fun, lik.new)

      if(al_i == 1)
      {
        beta_best = beta_new
        LambdaY_best = Lambda_Y_new
        lambdaY_best = lambda_Y_new
        lambdaYd_best = lambda_Ydot_new
        lik_best = lik.new
        alpha_best = alpha.i
        #print(beta_best)
      }else if((al_i>1)&(lik.fun[al_i-1]<lik.new))
      {
        beta_best = beta_new
        LambdaY_best = Lambda_Y_new
        lambdaY_best = lambda_Y_new
        lambdaYd_best = lambda_Ydot_new
        lik_best = lik.new
        alpha_best = alpha.i
        #print(beta_best)
        #print(alpha.i)
      }

    }
    result_lik = rbind(lik.fun, result_lik)
    result_alpha = c(alpha_best,result_alpha)
    result_lambda_Ydot = lambdaYd_best

    result = matrix(c(beta_best,LambdaY_best[I_1],LambdaY_best[I_2],1-mean(delta)), nrow = 1)
    colnames(result) = c(paste("beta", c(1:X.size)), "Lambda(1)", "Lambda(2)", "censor")
    Lam_Yd_I = matrix(rep(Y_dot,m), nrow = m ,byrow = FALSE)>=( matrix(rep(Y_dot,m), nrow = m, byrow = TRUE) )
    result_Lambda_Ydot = rowSums(Lam_Yd_I * matrix(rep(lambdaYd_best,m), nrow = m, byrow = TRUE))

    cat("EM finished, now start computing standard error")
    #-------------------------SE-------------------------
    Y_eq_Yhat = EM.res$Y_eq_Yhat
    Y_geq_Yhat =  EM.res$Y_geq_Yhat

    dgamma_Q = stats::dgamma(nodes, shape = 1/alpha_best, scale = alpha_best)

    beta_X = c( X%*%matrix(beta_best, ncol = 1))
    #-------------------------E(S)----------------------
    ES_nu = matrix(0, nrow = l+X.size, ncol = n)


    #L1/beta1         L2/beta1     ...   Ln/beta1
    #L1/beta2         L2/beta2     ...   Ln/beta2
    #L1/Lambda{Yhat 1} L2/Lambda{Yhat 1} ... Ln/Lambda{Yhat 1}
    #...
    #L1/Lambda{Yhat l} L2/Lambda{Yhat l} ... Ln/Lambda{Yhat l}
    #


    Y_mat = matrix(rep(Y,l), nrow = l, byrow = TRUE)
    #Y1 Y2 ... Yn
    #Y1 Y2 ... Yn  l*n

    Y_hat_mat = matrix(rep(Y_hat,n), nrow = l, byrow = FALSE)
    #Yhat 1 ... Yhat 1
    #Yhat 2 ... Yhat 2
    #...
    #Yhat l ... Yhat l   l*n
    ES = c()
    ES_g = matrix(0, nrow = l+X.size, ncol = n)
    ES_de_c = rep(0,n)

    Lambda_b_nozero = lambdaY_best
    Lambda_b_nozero[which(Lambda_b_nozero==0)]=10^(-10)
    for(g in 1:Q)
    {
      #print(g)

      xi_i = nodes[g]
      ES_de_each = whts[g] * (xi_i* lambdaY_best* exp(beta_X))^(delta) *
        exp(-xi_i* LambdaY_best* exp(beta_X)) *
        dgamma_Q[g] *
        exp(xi_i)

      ES_nu[1:X.size,] =  ES_nu[1:X.size,] + t((X * (delta-xi_i*LambdaY_best*exp(beta_X)))*
                                                 ES_de_each)
      ES_nu[(X.size+1):(l+X.size),] = ES_nu[(X.size+1):(l+X.size),] +
        t( ES_de_each *
             t(
               matrix(rep(delta,l), nrow = l, byrow = TRUE)*
                 (Y_eq_Yhat)/
                 matrix(rep(Lambda_b_nozero,l), nrow = l, byrow = TRUE)-
                 matrix(rep(xi_i*exp(beta_X),l), nrow = l, byrow = TRUE)*
                 (Y_geq_Yhat)

             )
        )
      ES_de_c = ES_de_c + ES_de_each
    }
    cat("for loop finished.")

    ES_g = ES_nu / matrix(rep(ES_de_c,(l+X.size)), nrow = (l+X.size), byrow = TRUE)
    ESES = ES_g%*%t(ES_g)

    #---------------------E(SST),E(B)---------------

    list_ESSEB = ESSEB2(l, n, Q, X.size, nodes, whts, lambdaY_best, LambdaY_best,
                        beta_X, delta, dgamma_Q, Lambda_b_nozero,
                        Y_eq_Yhat, Y_geq_Yhat, X)

    I = (list_ESSEB$EB-list_ESSEB$ESS+ESES)
    SE_square = diag(solve(I))
    SE = sqrt(SE_square)

    e1 = c(rep(0,X.size),(Y_hat <= 1))
    e2 = c(rep(0,X.size),(Y_hat <= 2))

    #standard error of estimated Lambda(1) and Lambda(2)
    #SE_Lambda = sqrt( diag( rbind(e1,e2) %*% solve(I) %*% cbind(e1,e2) ) )
    SE_Ydot_m = cbind(matrix(0, nrow = m, ncol = X.size),
                      (matrix(rep(Y_hat,m),nrow = m, byrow = TRUE )<= matrix(rep(Y_dot, l), nrow = m, byrow = FALSE) ) )
    #standard error of estimated Lambda(Y')
    SE_Ydot = sqrt(diag(SE_Ydot_m%*% solve(I) %*%t(SE_Ydot_m)))
    cat("SE finished. Start plotting")

  }else if((length_alpha == 1)&(alpha.c ==0)&(trsmodel == FALSE))
  {
    lik.fun = c()
    dataframe = data.frame(Y=Y, delta = delta, X1 = X[,1])
    if(ncol(X)>1)
    {
      for (i in 2:X.size)
      {
        dataframe = cbind(dataframe,X[,i])
      }
    }
    names(dataframe) = c("Y", "delta", paste("X", 1:X.size, sep = ""))
    for.sur1 = survival::Surv(time = Y, event = dataframe$delta)
    for.sur2 = paste("for.sur1~", paste("X",1:X.size ,sep = "", collapse = "+"), sep = "")
    for.sur = stats::as.formula(for.sur2 )
    res.cox = survival::coxph(for.sur, data = dataframe)

    su.res.cox = summary(res.cox)

    bas.frame = data.frame(t(rep(0,X.size)))
    names(bas.frame) = paste("X",1:X.size,sep = "")
    sur = survival::survfit(res.cox, newdata = bas.frame)
    su.sur = summary(sur)

    beta_best = su.res.cox$coefficients[,1]
    result_Lambda_Ydot = sur$cumhaz
    result_lambda_Ydot = result_Lambda_Ydot-c(0,result_Lambda_Ydot[1:(length(result_Lambda_Ydot)-1)])
    SE = c()
    SE[1:X.size] = su.res.cox$coefficients[,3]
    SE_Ydot = sur$std.err
    result_alpha = c(0)
    #for survival package, 0.1819217704 and 0.1819217601  are treated as the same number.
    #but for unique(), they are not.
    Y_dot = sur$time
  }else
  {
    cat("wrong alpha<0")
  }


  u_975 = stats::qnorm(0.975)

  Lambda_Ydot_upper = result_Lambda_Ydot + u_975*SE_Ydot
  Lambda_Ydot_lower = result_Lambda_Ydot - u_975*SE_Ydot

  if(plot.Lamb)
  {
    plot(x = rep(Y_dot,each=2)[-1], y = rep(result_Lambda_Ydot,each=2)[1:(2*m-1)], type = "l",lty=1, col = "blue", ylim = c(0,0.7),xlim = c(0,3),
         ylab = "Lambda(Y)", xlab = "Y")
    graphics::lines(x = rep(Y_dot,each=2)[-1], y = rep(Lambda_Ydot_upper,each=2)[1:(2*m-1)], type = "l",lty=2, col = "red")
    graphics::lines(x = rep(Y_dot,each=2)[-1], y = rep(Lambda_Ydot_lower,each=2)[1:(2*m-1)], type = "l",lty=2, col = "red")
    graphics::legend("topleft",legend = c("estimate","confidence interval"),lty = c(1,2), col = c("blue", "red"))
  }


  p.beta = 2*stats::pnorm(abs(beta_best/SE[1:X.size]) , lower.tail = FALSE)
  p.Lambda = 2*stats::pnorm(abs(result_Lambda_Ydot/SE_Ydot) , lower.tail = FALSE)

  if(show_res)
  {cat(list( beta.est = beta_best, SE.beta = SE[1:X.size], Ydot = c(Y_dot[1:5], "..."),Lambda.Ydot = c(result_Lambda_Ydot[1:5], "..."),
              SE.Ydot = c(SE_Ydot[1:5], "..."), choose.alpha = result_alpha,
              p.Lambda = c(p.Lambda[1:5], "...")))
  }else
  {
    cat("done.")
  }



  return(list( beta.est = beta_best, SE.beta = SE[1:X.size], Ydot = Y_dot,Lamb.est = result_Lambda_Ydot,
               lamb.est = result_lambda_Ydot, SE.Ydot = SE_Ydot, choose.alpha = result_alpha, all.alpha = alpha.c,
               Lamb.upper = Lambda_Ydot_upper, Lamb.lower = Lambda_Ydot_lower, p.beta = p.beta,
               p.Lambda = p.Lambda, lik = lik.fun))

}
