# GLT_prior without intercept

GLT_prior = function (y, X, burn = burn, nmc = nmc, thin = thin,BCM_sampling = TRUE){
  
  #################################################################################
  ############################## Model Setting ####################################
  #################################################################################
  # y = X*beta + epsilon, X: N*p, epsilon ~ N[0, sigma.sq]
  # For each j in 1:p, we have
  # beta.j | sigma.sq, lambda.j ~ N[0, lambda.j^2 * sigma.sq]
  # lambda.j|tau, xi ~ GPD[tau, xi], tau, xi > 0 
  # tau | xi ~ IG[p/xi + 1, 1] : Inverse Gamma
  # xi ~ log-Normal[mu, rho.sq]*I[1/2, infty](xi)
  # mu: estimated via POT method
  # rho.sq: fixed with some small number like 0.001
  #################################################################################
  #################################################################################
  
  #############################################################################
  # Setup for simulation environments
  # burn : number of burned samples, nmc = number of poseterior samples
  Tot.S = burn + nmc # Total number of posterior samples 
  S = (Tot.S - burn)/thin # effective number of sample size
  N = nrow(X) # number of data
  p = ncol(X) # number of covariates
  #############################################################################
  
  # Hyperparameter specification (Used for Elliptical Slice Sampler centered by Hill estimator)
  # These values will be automatically optimized, hence, do NOT change these values
  k.index = ncol(X)
  rho.sq = 0.001
  
  #############################################################################
  # Make a room to store MCMC samples
  # 1. Rooms for storages of effective posterior samples
  beta.vec = matrix(0, p, S) 
  lambda.vec = matrix(0, p, S) 
  tau = rep(0,S)
  sigma.sq = rep(0,S)
  xi = rep(0,S)
  
  # 2. Initial values
  beta.vec.out = rep(1,p)
  lambda.vec.out = rep(1,p)
  tau.out = 0.1 
  sigma.sq.out = 1
  eta.vec.out = rep(1,p)
  u.lambda.vec.out = rep(0.5,p)
  v.lambda.vec.out = rep(0.5,p)
  u.tau.vec.out = rep(0.5,p)
  v.tau.out = 0.5
  xi.out = 2
  eta.xi.out = log(2)
  
  # 3. Predefine sampling algorithm for multivariate normal distribution
  # Method 1: Rue's algorithm
  rmvt.Rue.sample = function(Q, b){
    # Goal:
    # Sample from beta = N[mu, Sigma]
    # such that
    # mu = solve(Q)%*%b # p dim vector
    # Sigma = solve(Q) # p by p matrix
    
    # Q : p by p matrix, Precision matrix
    # b : p dim vector
    
    # Useful 1. when n >>> p (i.e. Number of data is more than covariate)
    #        2. Possibly, We want to utilize precision structure
    #        3. Use cholesky and avoiding inverting precision matrix by solving three linear equations
    p = nrow(Q)
    
    # Step 1
    L = t(chol(Q))
    
    # Step 2
    z = stats::rnorm(n = p)
    
    # Step 3
    mu = solve(t(L), z)
    
    # Step 4
    v = solve(L, b)
    
    # Step 5
    theta = solve(t(L), v)
    
    # Step 6
    beta = mu + theta
    
    return(beta)
    
  }
  # Method 2: Anirban's algorithm
  rmvt.Anirban.sample = function(PHI, alpha, D.entries){
    # Goal:
    # Sample from N[mu, Sigma]
    # such that
    # mu = Sigma%*%t(PHI)%*%alpha # p dim vector
    # Sigma = solve( t(PHI)%*%PHI + solve(D) ) # p by p matrix
    # PHI : n by p matrix 
    # alpha : n dim vector 
    # D : p by p diagonal matrix 
    # Useful 1. when p >>> n (i.e. Number of covariates are more than number of data)
    #        2. We know how to sample from N[0,D] and solve(D) is easy to compute
    #        3. We don't want to use precisio matrix based sampling such as Rue's algorithm
    
    p = length(D.entries)
    n = nrow(PHI)
    D = diag(D.entries)
    # Step1
    
    u = sqrt(D.entries)*(stats::rnorm(n = p))
    
    delta = stats::rnorm(n = n)
    
    # Step2
    v = PHI%*%u + delta
    # Step3
    w = solve(PHI%*%D%*%t(PHI) + diag(n), alpha - v )
    # Step4
    theta = u + D%*%t(PHI)%*%w
    
    return(theta)
    
  }
  
  # 4. Bayesian computation (Gibbs sampler + Elliptical Slice Sampler centered by the Hill estimator)
  # This sample follows the Gibbs sampler given in the Supplementary Material Section C. Posterior computation
  for (s in 1:(Tot.S -1)){
    
    # Iterations
    {
      
      # Gibbs sampler - Section C - Step 1
      # Updating beta.vec (batch update)
      I_N = diag(N)
      l0 = rep(0, p)
      l2 = rep(1, p)
      
      if (BCM_sampling == TRUE){
        # Anirban's algorithm (fast-sampling algorithm developed by Anirban Bhattacharya (2016, Biometrika))
        U = as.numeric(lambda.vec.out^2) * t(X)
        u = rnorm(l2, l0, lambda.vec.out)
        v = X %*% u + rnorm(N) 
        v_star = solve( (X %*% U + I_N), (y/sqrt(sigma.sq.out)  - v ))
        beta.vec.out = sqrt(sigma.sq.out) * (u + U %*% v_star) 
        # beta.vec.out = rmvt.Anirban.sample(PHI = X/sqrt(sigma.sq.out), alpha = y/sqrt(sigma.sq.out), D.entries = sigma.sq.out*((lambda.vec.out)^2) )
      }
      if (BCM_sampling == FALSE){
        # Rue's algorithm
        inv.L = diag(c( 1/(lambda.vec.out^2) ))
        L.rue = t(chol(  (1/sigma.sq.out)*(t(X)%*%X + inv.L)  ))
        mu.rue = solve(t(L.rue), stats::rnorm(n = p))
        v.rue = solve(L.rue, (1/sigma.sq.out)*(t(X)%*%y))
        theta.rue = solve(t(L.rue), v.rue)
        beta.vec.out = mu.rue + theta.rue
        # Rue's algorithm
        #inv.L = diag(c( 1/(lambda.vec.out^2) ))
        #beta.vec.out = rmvt.Rue.sample(Q =(1/sigma.sq.out)*(t(X)%*%X + inv.L), b = (1/sigma.sq.out)*(t(X)%*%y))
      }
      
      # # Gibbs sampler - Section C - Step 2
      # Updating sigma.sq 
      # Numerical stability idea from horseshoe function by Van der Pas
      E_1 = max(t(y - X %*% beta.vec.out) %*% (y - X %*% beta.vec.out), (1e-10))
      E_2 = max(sum(beta.vec.out^2/(lambda.vec.out)^2), (1e-10) ) #Alternatively, we can also use: inv.L = diag(1/( lambda.vec.out )^2 ) ;E_2 = max(t(beta.vec.out)%*%inv.L%*%beta.vec.out, (1e-10) )
      sigma.sq.out = 1/rgamma(n = 1, shape = (N+p)/2, rate = (1/2)*( E_1 + E_2 )  )
      
      
      # Gibbs sampler - Section C - Step 3
      # Updating lambda.vec
      # Sketch of idea: 
      #############################################################
      ##### lambda --> eta --> (eta,u) --> (eta, v, u)      #######
      #####      1. T   2. PX         3. PX                 #######
      # Note: 1. T means transformation: eta = lambda^2     #######
      #       2. PX means parameter expansion               #######
      #############################################################
      
      ########################################################################################
      # Par-computing                                                                  #######
      # Idea:                                                                          #######
      # lambda.vec => eta.vec => u.lambda.vec => v.lambda.vec => eta.vec => lambda.vec #######
      #          Step 3-1  Step 3-2       Step 3-3          Step 3-4    Step 3-5       #######
      ########################################################################################
      
      # Define all functions necessary for Step 3
      if (s == 1){
        ft.g = function(eta, tau, xi){
          res = ( ( sqrt(eta) )^(1/xi + 1) )*( ( tau + xi*sqrt(eta) )^(-(1/xi + 1) ) )
          return(res)
        }    
        ft.inv.g = function(u, tau, xi){
          res = ( ( tau )/( u^(-( xi/(1+ xi) )) - xi ) )^2
          return(res)
        }  
        ft.m = function(beta, sigma.sq){
          res = (beta^2)/(2*sigma.sq)
          return(res)
        }
        
        # functions necessary for updating v
        ft.A = function(xi){
          res = (1/2)*(1/xi + 1)
          return(res)
        }
        ft.B = function(beta, sigma.sq){
          res = ft.m(beta = beta, sigma.sq = sigma.sq)
          return(res)
        }
        ft.C = function(u, tau, xi){
          res = ft.inv.g(u = u, tau = tau, xi = xi)
          return(res)
        }
        indicator.greater.than.C = function(eta, C){
          
          if ( eta > C ){
            res = 1
          } else {
            res = 0
          }
          
          return(res)
        }
        ft.g.1 = function(eta,B) {
          res = exp(-B/eta)
          return(res)
        }
        ft.h.1 = function(eta, B, C){
          res = ft.g.1(eta = eta, B = B)*indicator.greater.than.C(eta = eta, C = C)
          return(res)
        }
        ft.inv.g.1 = function(v, B){
          res = -B/log(v)
          return(res)
        }
        hat.C = function(v, B, C){
          res = max(C, ft.inv.g.1(v = v, B = B))
          return(res)
        }
        Z.v = function(v, A, B, C){
          res = ( hat.C(v = v, B = B, C = C )^(-A))/A
          return(res)
        }
        inv.CDF.v = function(U, v, C, A, B){
          res = ( hat.C(v = v, B = B, C = C)^(-A) -
                    Z.v(v = v, A = A, B = B, C = C)*
                    A*U )^(-1/A)
          return(res)
        }
        
        # Updating functions 
        
        # We need three updating functions
        # u ~ pi(u|eta, -)
        u.update = function(previous.eta, tau = tau, xi = xi){
          res = stats::runif(n = 1, min = 0, max = ft.g(eta = previous.eta, tau = tau, xi = xi))
          return(res)
        }
        # v ~ pi(v|eta, u, -)
        v.update = function(previous.eta, previous.u, beta, sigma.sq, tau, xi){
          B = ft.B(beta = beta, sigma.sq = sigma.sq)
          C = ft.C(u = previous.u, tau = tau, xi = xi)
          res = stats::runif(n = 1, min = 0, max = ft.h.1(eta = previous.eta, B = B, C = C) )
          return(res)
        }
        # eta ~ pi(eta|v, u, -)
        eta.update = function(previous.v, previous.u, beta, sigma.sq, tau, xi){
          A = ft.A(xi = xi)
          B = ft.B(beta = beta, sigma.sq = sigma.sq)
          C = ft.C(u = previous.u, tau = tau, xi = xi)
          # Sampling 
          U = stats::runif(n = 1, min = 0, max = 1)  
          res = inv.CDF.v(U = U, v = previous.v, A = A, B = B, C = C)
          return(res)
        }  
        
        # Bulk update functions 
        u.lambda.vec.bulk = function(j){
          res = u.update(previous.eta = eta.vec.out[j], tau = tau.out, xi = xi.out)
          return(res)
        }  
        
        v.lambda.vec.bulk = function(j){
          res = v.update(previous.eta = eta.vec.out[j], previous.u = u.lambda.vec.out[j], 
                         beta = beta.vec.out[j], sigma.sq = sigma.sq.out, tau = tau.out, xi = xi.out) 
          return(res)
        }  
        eta.vec.bulk = function(j){
          res = eta.update(previous.v = v.lambda.vec.out[j], previous.u = u.lambda.vec.out[j], beta = beta.vec.out[j], 
                           sigma.sq = sigma.sq.out, tau = tau.out, xi = xi.out) 
          return(res)
          
        }  
        
      }
      
      # Step 3-1: lambda.vec => eta.vec
      eta.vec.out = (lambda.vec.out)^2
      # Step 3-2: eta.vec => u.lambda.vec
      u.lambda.vec.out =  unlist( lapply(1:p,  u.lambda.vec.bulk) )  # output is p dim vector
      # Step 3-3: u.lambda.vec => v.lambda.vec
      v.lambda.vec.out = unlist( lapply(1:p,  FUN = v.lambda.vec.bulk) )   # output is p dim vector
      # Step 3-4: v.lambda.vec => eta.vec
      eta.vec.out = unlist( lapply(1:p,  FUN = eta.vec.bulk) ) # output is p dim vector
      # Step 3-5: eta.vec => lambda.vec
      lambda.vec.out  = sqrt(eta.vec.out)
      #######################################################################################################################################################################
      
      # Gibbs sampler - Section C - Step 4
      # Updating tau
      
      # Sketch of idea: 
      ########################################################################################################################
      ###### (tau) ===> (tau, u_1:p) ===> (tau, v, u_1:p)                                                              ########
      ######       PX_1              PX_2                                                                              ########
      # Note: 1. PX means parameter expansion                                                                          ########
      #       2. to distinguish from notation used in lambda we explicitely state "tau" as suffix in used function     #######
      #       3. To distinguish functions used in PX_1 and PX_2, we state "1" or "2" to represent the PX_1, PX_2       #######
      #              For e.g., ft.g.tau.1 is the slicer g(tau) for the PX_1 used in slice sampler in tau               #######
      ########################################################################################################################
      
      #############################################################
      # Par-computing                                        ######
      # Idea:                                                ######
      #    ====>   u.tau.vec =====> v.tau =====> tau        ######
      #   Step 4-1           Step 4-2     Step 4-2           ######
      #############################################################
      
      # Define all functions necessary for Step 4
      if (s == 1){
        {
          
          # PX_1
          
          # functions necessary for updating u
          ft.g.tau.1 = function(tau, xi, lambda.j){
            res = ( tau + xi*lambda.j)^(-1/xi - 1)
            return(res)
          }    
          ft.inv.g.tau.1 = function(u, xi, lambda.j){
            res = u^(-xi/(1+xi)) - xi*lambda.j
            return(res)
          }  
          
          hat.C.tau = function(u.tau.vec, xi, lambda.vec){
            temp = ft.inv.g.tau.1(u = u.tau.vec, xi = xi, lambda.j = lambda.vec)
            res = min(c(temp))
            return(res)
          }
          
          # PX_2
          # functions necessary for updating v
          indicator.smaller.than.C = function(tau, C.tau){
            if ( tau < C.tau ){
              res = 1
            } else {
              res = 0
            }
            
            return(res)
          }
          
          ft.g.tau.2 = function(tau) {
            res = exp(-1/tau)
            return(res)
          }
          
          ft.h.tau.2 = function(tau, C.tau){
            res = ft.g.tau.2(tau = tau)*indicator.smaller.than.C(tau = tau, C.tau = C.tau)
            return(res)
          }
          
          ft.inv.g.tau.2 = function(v){
            res = -1/log(v)
            return(res)
          }
          inv.CDF.v.tau = function(U, v, u.tau.vec, xi, lambda.vec){
            res = 1/(  1/ft.inv.g.tau.2(v = v)  + U*( 1/hat.C.tau(u.tau.vec = u.tau.vec, xi = xi, lambda.vec = lambda.vec)   -  1/ft.inv.g.tau.2(v = v)  ) )
            return(res)
          }
          
          
          # Updating functions
          # u.j ~ pi(u.j|u_(-j), tau, -) = pi(u.j|tau, -), j in 1:p, individual u_j update function
          u.tau.update = function(previous.tau, xi, lambda.j){
            res = stats::runif(n = 1, min = 0, max = ft.g.tau.1(tau = previous.tau, xi = xi, lambda.j = lambda.j))
            return(res)
          }
          # v ~ pi(v|tau, u_(1:p), -)
          v.tau.update = function(previous.tau, previous.u.tau.vec, xi, lambda.vec ){
            res = stats::runif(n = 1, min = 0, max = ft.h.tau.2(tau = previous.tau, 
                                                         C.tau = hat.C.tau(u.tau.vec = previous.u.tau.vec, xi = xi, lambda.vec = lambda.vec)))
            return(res)
          }
          # tau ~ pi(tau|v, u_(1:p), -)
          tau.update = function(previous.v, previous.u.tau.vec, xi, lambda.vec){
            
            # Sampling 
            U = stats::runif(n = 1, min = 0, max = 1)  
            res = inv.CDF.v.tau(U = U, v = previous.v, u.tau.vec = previous.u.tau.vec, xi = xi, lambda.vec = lambda.vec)
            return(res)
            
          }
          
        }  
        u.tau.vec.bulk = function(j){
          res = u.tau.update(previous.tau = tau.out, lambda.j = lambda.vec.out[j], xi = xi.out)
          return(res)
        } 
        
      }
      # Step 4-1: ==> u.tau.vec
      u.tau.vec.out = unlist( lapply(1:p,  FUN = u.tau.vec.bulk) ) # output is p dim vector
      # Step 4-2: u.tau.vec =====> v.tau (p. 21: Needs a numerical Stability)
      SMN = 5e-324 # SMN: SMallest Number in R
      v.tau.out = max(v.tau.update(previous.tau = tau.out, previous.u.tau.vec = u.tau.vec.out, lambda.vec = lambda.vec.out, xi = xi.out), SMN ) # v.tau.out = v.tau.update(previous.tau = tau.out, previous.u.tau.vec = u.tau.vec.out, lambda.vec = lambda.vec.out, xi = xi.out)
      # Step 4-3:. v.tau =====> tau
      tau.out = tau.update(previous.v = v.tau.out, previous.u.tau.vec = u.tau.vec.out, lambda.vec = lambda.vec.out, xi = xi.out) 
      
      
      
      # Gibbs sampler - Section C - Step 5
      # Updating xi 
      
      # Sketch of idea: 
      ########################################################################################################################
      ###### xi ===> eta.xi ===> xi                                                                                     ######
      ######     T           T                                                                                          ######
      # Note: 1. T means transformation (eta.xi = log(xi))                                                              ######
      #       2. to distinguish from notation used in lambda we explicitely state "xi" as suffix in used function       ######
      #       3. Obtatin mu.hat by using Hill estimator (Need to specify some number of orderstatistics used)           ######                                     
      #       4. Use ESS: eta.xi ~ pi(eta.xi|-) ~ V_p(eta.xi)*I[-log(2) < eta.xi]*N_1[eta.xi | mu.hat, rho.sq = 0.001 ] ######
      ########################################################################################################################
      
      #############################################################################
      # Par-computing                                                        ######
      # Idea:                                                                ######
      # xi   ====>     eta.xi ==== [ CALIBRATION of mu ====> ESS  =====> xi  ######
      #                                                                      ######
      #############################################################################
      
      # Define all functions necessary for Step 5
      if (s == 1){
        # Hill estimator
        Hill.estimator = function(k) {
          res = mean( log( ordered.lambda.vec.out[1:(k-1)] ) ) -
            log(ordered.lambda.vec.out[k])  
          return(res)
        }
        
        log.ratio.vol.epllipsoids = function(x, y, tau, lambda.vec){
          #  res = log(V_p(x)/V_p(y)) # x : new, y : past
          # s.t V_p(x) = (pi^(p/2))*prod(radiouses)/gamma(p/x + 1)
          # when x = 2, V_p(x) becomes the proper volumn of p-dim ellipsoid
          # note that lgamma is very stable built-in function in R
          p = length(lambda.vec)
          res = lgamma(p/y + 1) - lgamma(p/x + 1) -(1+1/x)*sum(log(tau + x*lambda.vec))+ (1+1/y)*sum(log(tau + y*lambda.vec)) 
          
          return(res)
        }  
      }
      # Step 5-1: xi   ====>   eta.xi
      eta.xi.out = log(xi.out)
      # [ CALIBRATION of mu ]
      {
        ordered.index = order(lambda.vec.out,decreasing = T)
        ordered.lambda.vec.out = lambda.vec.out[ordered.index]
        mu.hat = log( mean(unlist( lapply( floor(k.index*1/10): floor(k.index*9/10) ,  FUN = Hill.estimator) )) )
      }
      # [ Elliptical Slice Sampler ]
      {
      # a. 
      nu = rnorm(n = 1, mean = mu.hat, sd = sqrt(rho.sq))
      # b. 
      U = stats::runif(n = 1, min = 0, max = 1)
      # c. 
      theta = stats::runif(n = 1, min = -pi, max = pi)
      # d. 
      eta.xi.new = (eta.xi.out - mu.hat)*cos(theta) + (nu - mu.hat)*sin(theta) + mu.hat
      # e. Restriction check (1/2 < xi, -log(2) < eta.xi.new)
      # e - i. (Must be rejected)
      if (-log(2) > eta.xi.new){
        MH.prob = 0 
      }
      # e - ii. (Proceed for ESS)
      if (-log(2) < eta.xi.new){
        temp = log.ratio.vol.epllipsoids(x = exp(eta.xi.new), y = exp(eta.xi.out), tau = tau.out, lambda.vec = lambda.vec.out)  
        MH.prob = min(exp(temp), 1)
      }
      # f.
      # f - i. Accepted
      if (U < MH.prob){
        eta.xi.out = eta.xi.new 
      }
      # f - ii. Rejected => Angular correction
      if (U > MH.prob){
        # Idea: If reject, then sample again from narrowing interval which collapsing to zero and eventually, acceptence occurs
        # NOTE: U which were sampled in b. is continuously used
        #       However, acceptance ratio will be change
        # Start by setting theta.min and theta.max
        # Rule1: redefine [theta.min = -pi, theta.max = pi] which were used in c.
        theta.min = -pi
        theta.max = pi
        while( U > MH.prob ){
          # Rule2: redefine [theta.min, theta.max] according to realized theta ~ U[theta.min, theta.max]
          if (theta > 0 ){
            theta.max = theta
          }
          if (theta < 0 ){
            theta.min = theta
          }
          ################################################
          # c.
          theta = stats::runif(n = 1, min = theta.min, max = theta.max)
          # d.
          eta.xi.new = (eta.xi.out - mu.hat)*cos(theta) + (nu - mu.hat)*sin(theta) + mu.hat
          # e.
          # e-i. (Must be rejected)
          if (-log(2) > eta.xi.new){
            MH.prob = 0
          }
          # e-ii. (Proceed for ESS)
          if (-log(2) < eta.xi.new){
            temp = log.ratio.vol.epllipsoids(x = exp(eta.xi.new), y = exp(eta.xi.out), tau = tau.out, lambda.vec = lambda.vec.out)
            MH.prob = min(exp(temp), 1)
          }
          
        }
        eta.xi.out = eta.xi.new  
      
      }
      
      }
      # Step 5-3: [ CALIBRATION of mu ====> ESS]  =====> xi
      xi.out = exp(eta.xi.out)
      
    }
    # Print out progression
    if (s%%1000 == 0){
      print(s)  
    }
    # Store output
    s = s + 1 # Index modification for storing
    if ( s > (burn) && s%%thin == 0){
      beta.vec[,(s - burn)/thin] = beta.vec.out
      lambda.vec[,(s - burn)/thin] = lambda.vec.out
      tau[(s - burn)/thin] = tau.out
      sigma.sq[(s - burn)/thin] = sigma.sq.out
      xi[(s - burn)/thin] = xi.out
      }
    }
  
  # 5. Print out posterior samples
  res = list(beta.vec = beta.vec, lambda.vec = lambda.vec, tau = tau, sigma.sq = sigma.sq, xi = xi )
  return(res)
  
}