#'@title  Adaptive LASSO for Semiparametric Transformation Models.
#'
#'@description   Select the important variables in semiparametric transformation
#'  models for right censored data using adaptive lasso.
#'
#' @return a list containing
#' \tabular{lccl}{
#'   \code{beta_res} \tab\tab\tab the estimated \eqn{\beta} with the selected tuning parameter \eqn{\lambda} \cr
#'   \code{GCV_res} \tab\tab\tab the value of GCV  with the selected tuning parameter \eqn{\lambda} \cr
#'   \code{lamb_res} \tab\tab\tab the selected tuning parameter \eqn{\lambda} \cr
#'   \code{beta_all} \tab\tab\tab estimated \eqn{\beta} with all tuning parameters \cr
#'   \code{CSV_all} \tab\tab\tab value of GCV  with all tuning parameters \cr
#'   \code{skip_para} \tab\tab\tab a list containing the \eqn{\lambda} and the number of adaptive lasso iteration when adaptive lasso doesn't work. \cr
#' }
#'
#'@details The initial value of the coefficient \eqn{\beta} used as the adapting
#'  weights is EM estimator, which is computed by the function \code{EM_est}.
#'  The tuning parameter \eqn{\lambda} is data-dependent and we select it using
#'  generalized crossvalidation. There may be some errors for small
#'  \eqn{\lambda}, in which case the \eqn{\lambda} and the number of adaptive
#'  lasso iteration are recorded in the \code{skip_para}.
#'
#'@examples
#'
#'
#'  if(!requireNamespace("MASS", quietly = TRUE))
#'  {stop("package MASS needed for this example. Please install it.")}
#'
#'  gen_lasdat = function(n,r,rho,beta_true,a,b,seed=66,std = FALSE)
#' {
#'
#'   set.seed(seed)
#'   beta_len = length(beta_true)
#'   beta_len = beta_len
#'   sigm = matrix(0, nrow = beta_len, ncol = beta_len)
#'   for(i in 1:(beta_len-1))
#'   {
#'     diag(sigm[1:(beta_len+1-i),i:beta_len]) = rho^(i-1)
#'   }
#'   sigm[1,beta_len] = rho^(beta_len-1)
#'   sigm[lower.tri(sigm)] = t(sigm)[lower.tri(sigm)]
#'
#'   Z = MASS::mvrnorm(n, mu = rep(0, beta_len), Sigma = sigm)
#'   beta_Z.true = c(Z %*% beta_true)
#'   U = runif(n)
#'   if(r>0)
#'   {
#'     t = ((U^(-r)-1)/(a*r*exp(beta_Z.true)))^(1/b)
#'   }else if(r == 0)
#'   {
#'     t = (-log(U)/(a*exp(beta_Z.true)))^(1/b)
#'     #t = (exp(-log(U)/(0.5 * exp(beta_Z.true))) - 1)
#'   }
#'   C = runif(n,0,8)
#'   Y = pmin(C,t)
#'   delta_i = ifelse( C >= t, 1, 0)
#'   if(std)
#'   {
#'     Z = apply(Z,2,normalize)
#'   }
#'   return(list(Z = Z, Y = Y, delta_i = delta_i,censor = mean(1-delta_i)))
#' }
#'
#' now_rep=1
#' dat = gen_lasdat(100,1,0.5,c(0.3,0.5,0.7,0,0,0,0,0,0,0),2,5,seed= 6+60*now_rep,std = FALSE)
#' Z = dat$Z
#' Y = dat$Y
#' delta_i = dat$delta_i
#'
#' tra_ala = trans_alasso(Z,Y,delta_i,lamb_vec = c(5,7),r=1)
#' tra_ala$GCV_res
#' tra_ala$beta_res
#' tra_ala$lamb_res
#'
#'
#'@references Xiaoxi, L. , & Donglin, Z. . (2013). Variable selection in
#'  semiparametric transformation models for right-censored data. Biometrika(4),
#'  859-876.
#'
#'
#'@param Z the baseline covariates
#'@param Y observed event times
#'@param delta_i censoring indicator. If \eqn{Y} is censored, \code{delta_i}=0. If not, \code{delta_i}=1.
#'@param r parameter in transformation function
#'@param lamb_vec the grad of the tuning parameter \eqn{\lambda}
#'@param solu determines whether the solution path will be plotted. The default is TRUE.
#'
#'
#'@export
trans_alasso = function(Z,Y,delta_i,r,lamb_vec,solu = TRUE)
{
  fun_E = function(Y, Z, delta, beta_ini, lamb_ini, r, Q = 60)
  {

    gq = statmod::gauss.quad(Q, kind="laguerre")
    nodes = gq$nodes
    whts  = gq$weights
    dgamma_Q = stats::dgamma(nodes, shape = 1/r, scale = r)

    n = length(Y)
    Y_mat_nn = matrix(rep(Y,n), nrow = n, byrow = TRUE)
    Lamb_ini = colSums( lamb_ini * (Y <= Y_mat_nn) )

    beta_Z_ini = c(Z %*% beta_ini)
    exp_beZ = exp(beta_Z_ini)
    lam_expbeZ_m = matrix( rep(lamb_ini* exp_beZ, Q ), nrow = Q, byrow = TRUE )
    Lam_expbeZ_m = matrix( rep(Lamb_ini* exp_beZ, Q ), nrow = Q, byrow = TRUE )
    E_de_m = whts * t( t(nodes * lam_expbeZ_m)^delta ) *
      exp( -nodes * Lam_expbeZ_m + matrix(rep(log(dgamma_Q)+nodes,n), nrow = Q)   )
    E_nu_m = nodes * E_de_m  #there's no whts, because it's included in E_de_m

    E = colSums(E_nu_m) / colSums(E_de_m)
    return(E)
  }
  G = function(x,r)
  {
    if(r>0)
    {return(log(1+r*x)/r)}
    else if(r==0)
    {return(x)}
    else{cat("r<0")}
  }
  dG = function(x,r)
  {
    if(r>0)
    {return(1/(1+r*x))}
    else if(r==0)
    {return(1)}
    else{cat("r<0")}
  }
  ddG = function(x,r)
  {
    if(r>0)
    {return(-r/(1+r*x)^2)}
    else if(r==0)
    {return(0)}
    else{cat("r<0")}
  }
  fun_c = function(Y, Z, delta, beta_ini, lamb_ini, r)
  {
    n = length(Y)
    Y_j.m = matrix(rep(Y,n), nrow = n, byrow = TRUE)
    Y_i.m = matrix(rep(Y,n), nrow = n, byrow = FALSE)
    Y_jleqi = Y_j.m <= Y_i.m
    #Y_jgeqi = Y_j.m >= Y_i.m
    beta_Z_ini = c(Z %*% beta_ini)
    temp1 = Y_jleqi %*% (exp(beta_Z_ini) * lamb_ini)
    til.c = c( (delta == 0) * dG(temp1,r) +
                 (delta == 1) * (-ddG(temp1,r)/dG(temp1,r) + dG(temp1,r)) )
    return(til.c)
  }
  loglik = function(n,delta,z,beta,Y,E)
  {
    Y_m = matrix(rep(Y,n), nrow = n, byrow = TRUE)
    Ie = (t(Y_m)<=Y_m)*matrix(rep( E*exp(Z%*%beta),n  ),nrow = n, byrow = TRUE)
    return( sum( delta*(Z%*%beta - log(rowSums(Ie))) ) )
  }
  normalize = function(x)
  {
    y = (x-mean(x))/stats::sd(x)
    return(y)
  }
  n = length(Y)
  beta_len = ncol(Z)
  #------record skipped for loop------
  skip_iter = rep(0,beta_len)#each row contains all lamb of each of iteration
  skip_para = matrix(0,nrow = 1, ncol = 2)
  colnames(skip_para) = c("lamb","the number of iteration of adaptive lasso")
  #------initual value-----
  possibleError_ini <- tryCatch(
    {
      res_EM = EM_est(Y, Z, delta_i, alpha = r)
      beta_ini = c(res_EM$beta_new)
      lamb_ini = c(res_EM$lamb_Y)
    },error=function(e)
    {
      cat(e)
    })
  if(inherits(possibleError_ini, "error"))
  {
    cat("\n","wrong EM","\n")
    return(list(beta_res = NA,
                GCV_res = NA,
                lamb_res = NA,
                beta_all = NA,
                CSV_all = NA,
                lamb_all = NA,
                skip_iter = NA,
                skip_para = NA))
  }
  #------beta,GCV,loglik for each lamb------
  beta_all = matrix(0, nrow = length(lamb_vec), ncol = beta_len)
  GCV = c()
  loglik_vec = c()
  #------begin for loop: lamb------
  for (lamb_ind in 1:length(lamb_vec))
  {
    lamb = lamb_vec[lamb_ind]
    beta_latmp = beta_ini
    skip_loop = FALSE
    #------start adalasso
    for (ada_rep in 1:200)
    {
      possibleError <- tryCatch(
        {
          til_c = fun_E(Y = Y, Z = Z, beta_ini = beta_latmp, delta = delta_i, lamb_ini = lamb_ini ,r = r )
          #til_c = fun_c(Y = Y, Z = Z, beta_ini = beta_latmp, delta = delta_i, lamb_ini = lamb_ini ,r = r )
          Exp_betaZ_tmp = c(exp(Z%*%beta_latmp))*til_c
          mat_Y = matrix(rep(Y,n), nrow = n, byrow =FALSE)#each row are the same
          Del_m = (Exp_betaZ_tmp*( mat_Y >= t(mat_Y)))
          Del = (t(Z) - ((t(Z) %*% (Del_m))/matrix(rep(colSums(Del_m),beta_len), nrow = beta_len, byrow = TRUE) )) %*% (-delta_i)
          Del2 = ddloglik_transmdl(n,delta_i[order(Y)],Z[order(Y),],beta_latmp,til_c[order(Y)])
          X = chol(Del2)
          W = solve(t(X))%*%(Del2%*%beta_latmp-Del)
          ada_new = shoot(X = X, y = W, beta_ini_s = beta_ini, lamb_s = lamb)
        },error=function(e)
        {
          cat(e)
        })
      if(inherits(possibleError, "error"))
      {
        cat("\n","lamb:",lamb,"\n")
        skip_loop = TRUE
        cat("set skip to be TRUE")
        cat(skip_loop)
        #record the parameters in the end of iteration
        skip_para = rbind(skip_para,c(lamb,ada_rep))
        skip_iter[lamb_ind] = 1
        #stop loop
        ada_new=beta_latmp=rep(10^(5),beta_len)
      }

      if(max(abs(ada_new-beta_latmp)) < 10^(-5)){break}
      beta_latmp = ada_new
      #print(beta_latmp)
    }
    #------GCV------
    #if the loop is skipped, GCV will be set to 100000, which ensures that the corresponding GGV can not be chosen.
    if(skip_loop)
    {
      beta_all[lamb_ind,] = beta_latmp
      GCV[lamb_ind] = 10^(5)
    }else
    {
      beta_all[lamb_ind,] = beta_latmp
      A = diag(1/abs(beta_latmp*beta_ini))
      diag(A)[which(beta_latmp == 0)] = 0
      diag(A)[which(beta_ini == 0)] = 0
      p =sum(diag(  solve(Del2 + lamb*A) %*% Del2  ) )
      GCV[lamb_ind] = (-loglik(n,delta_i,Z,beta_latmp,Y,til_c)) / (n*(1-p/n)^2)
      cat("\n","lamb_ind:",lamb_ind)
      #print(beta_all[lamb_ind,])
      #loglik_vec[lamb_ind] = -loglik(n,delta_i,Z,beta_latmp,Y,til_c)
    }

  }


  #------solution path------
  #save.image("trans gcv n=100 r=0.5 rep=200.RData")
  if(solu)
  {
    solutionpath = t(beta_all)
    col_c = c(2:beta_len)
    plot(x = lamb_vec, y = solutionpath[1,], type = "l", col = col_c[1],
         xlim = c(0,max(lamb_vec)+20), ylim = c(min(solutionpath), max(solutionpath)+0.2))
    for (i in 2:beta_len)
    {
      graphics::lines(x = lamb_vec, y = solutionpath[i,], col = col_c[i])
    }

    graphics::legend("right",
           legend=paste("beta",c(1:beta_len), sep = ""),
           col=col_c,
           lty=1,lwd=2)

    graphics::abline(v=lamb_vec[which(GCV == min(GCV))],lwd=1,col="black")
  }
  return(list(beta_res = beta_all[which(GCV == min(GCV))[1],],
              GCV_res = GCV[which(GCV == min(GCV))[1]],
              lamb_res = lamb_vec[which(GCV == min(GCV))[1]],
              beta_all = beta_all,
              GCV_all = GCV,
              skip_para = skip_para))
}
