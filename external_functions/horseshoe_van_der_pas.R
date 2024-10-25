horseshoe_van_der_pas = function (y, X, method.tau = c("fixed", "truncatedCauchy", "halfCauchy"), 
                                  tau = 1, method.sigma = c("fixed", "Jeffreys"), Sigma2 = 1, 
                                  burn = 1000, nmc = 5000, thin = 1, alpha = 0.05) 
{
  ######################################################################################
  # Defalt choice for current GLT paper:
  #  
  #  - INPUT:
  # res = horseshoe_van_der_pas(y = y, X = X, method.tau = "halfCauchy",
  #                                           method.sigma = "Jeffreys",
  #                                           burn = burn, 
  #                                           nmc = nmc,
  #                                           thin = thin,
  #                                           alpha = 0.05) 
  #
  #
  #  - OUTPUT:
  #  thined.beta.vec = res$BetaSamples
  #  thined.lambda.vec = res$LambdaSamples
  #  thined.tau = res$TauSamples
  #
  ######################################################################################
  
  method.tau = match.arg(method.tau)
  method.sigma = match.arg(method.sigma)
  ptm = proc.time()
  N = burn + nmc
  effsamp = (N - burn)/thin
  n = nrow(X)
  p = ncol(X)
  
  Beta = rep(0, p)
  lambda = rep(1, p)
  sigma_sq = Sigma2
  
  betaout = matrix(0, p, effsamp)
  lambdaout = matrix(0, p, effsamp) # lambda also has to be shown
  tauout = rep(0, effsamp)
  sigmaSqout = rep(0, effsamp)
  if (p > n) 
    algo = 1
  else algo = 2
  I_n = diag(n)
  l0 = rep(0, p)
  l1 = rep(1, n)
  l2 = rep(1, p)
  if (algo == 2) {
    Q_star = t(X) %*% X
  }
  for (i in 1:N) {
    if (algo == 1) {
      lambda_star = tau * lambda
      U = as.numeric(lambda_star^2) * t(X)
      u = stats::rnorm(l2, l0, lambda_star)
      v = X %*% u + stats::rnorm(n)
      v_star = solve((X %*% U + I_n), ((y/sqrt(sigma_sq)) - v))
      Beta = sqrt(sigma_sq) * (u + U %*% v_star)
    }
    else if (algo == 2) {
      lambda_star = tau * lambda
      L = chol((1/sigma_sq) * (Q_star + diag(1/as.numeric(lambda_star^2), p, p)))
      v = solve(t(L), t(t(y) %*% X)/sigma_sq)
      mu = solve(L, v)
      u = solve(L, stats::rnorm(p))
      Beta = mu + u
    }
    eta = 1/(lambda^2)
    upsi = stats::runif(p, 0, 1/(1 + eta))
    tempps = Beta^2/(2 * sigma_sq * tau^2)
    ub = (1 - upsi)/upsi
    Fub = 1 - exp(-tempps * ub)
    Fub[Fub < (1e-04)] = 1e-04
    up = stats::runif(p, 0, Fub)
    eta = -log(1 - up)/tempps
    lambda = 1/sqrt(eta)
    if (method.tau == "halfCauchy") {
      tempt = sum((Beta/lambda)^2)/(2 * sigma_sq)
      et = 1/tau^2
      utau = stats::runif(1, 0, 1/(1 + et))
      ubt = (1 - utau)/utau
      Fubt = stats::pgamma(ubt, (p + 1)/2, scale = 1/tempt)
      Fubt = max(Fubt, 1e-08)
      ut = stats::runif(1, 0, Fubt)
      et = stats::qgamma(ut, (p + 1)/2, scale = 1/tempt)
      tau = 1/sqrt(et)
    }
    if (method.tau == "truncatedCauchy") {
      tempt = sum((Beta/lambda)^2)/(2 * sigma_sq)
      et = 1/tau^2
      utau = stats::runif(1, 0, 1/(1 + et))
      ubt_1 = 1
      ubt_2 = min((1 - utau)/utau, p^2)
      Fubt_1 = stats::pgamma(ubt_1, (p + 1)/2, scale = 1/tempt)
      Fubt_2 = stats::pgamma(ubt_2, (p + 1)/2, scale = 1/tempt)
      ut = stats::runif(1, Fubt_1, Fubt_2)
      et = stats::qgamma(ut, (p + 1)/2, scale = 1/tempt)
      tau = 1/sqrt(et)
    }
    if (method.sigma == "Jeffreys") {
      if (algo == 1) {
        E_1 = max(t(y - X %*% Beta) %*% (y - X %*% Beta), 
                  (1e-10))
        E_2 = max(sum(Beta^2/((tau * lambda))^2), (1e-10))
      }
      else {
        E_1 = max(t(y - X %*% Beta) %*% (y - X %*% Beta), 
                  1e-08)
        E_2 = max(sum(Beta^2/((tau * lambda))^2), 1e-08)
      }
      sigma_sq = 1/stats::rgamma(1, (n + p)/2, scale = 2/(E_1 + 
                                                            E_2))
    }
    if (i%%1000 == 0) {
      print(i)
    }
    if (i > burn && i%%thin == 0) {
      betaout[, (i - burn)/thin] = Beta
      lambdaout[, (i - burn)/thin] = lambda
      tauout[(i - burn)/thin] = tau
      sigmaSqout[(i - burn)/thin] = sigma_sq
    }
  }
  pMean = apply(betaout, 1, mean)
  pMedian = apply(betaout, 1, stats::median)
  pSigma = mean(sigmaSqout)
  pTau = mean(tauout)
  left <- floor(alpha * effsamp/2)
  right <- ceiling((1 - alpha/2) * effsamp)
  BetaSort <- apply(betaout, 1, sort, decreasing = F)
  left.points <- BetaSort[left, ]
  right.points <- BetaSort[right, ]
  result = list(BetaHat = pMean, LeftCI = left.points, RightCI = right.points, 
                BetaMedian = pMedian, Sigma2Hat = pSigma, TauHat = pTau, 
                BetaSamples = betaout, LambdaSamples  = lambdaout, TauSamples = tauout, Sigma2Samples = sigmaSqout)
  
  return(result)
  
}


