### auxiliary mixture sampling for the Poisson distribution
### Based on Frühwirth Schnatter's bayesf MATLAB toolbox
### Y               ...   dependent variable, N x 1 vector
### X               ...   independent variables, as a N x k matrix
### nburnin         ...   Number of burn in iterations in Gibbs sampler
### nsave           ...   Number of saved iterations in Gibbs sampler
### beta_prior_var  ...   Prior variance of beta, k x k matrix
### beta_prior_mean ...   Prior mean of beta k x 1 vector
### family          ...   Can be NegBin or Poisson
### alpha_prop      ...   Proposal density for neg. bin. alpha parameter; either uniform or truncnorm
require(truncnorm)

mixturemcmc <- function(Y,X,
                        nburnin = 1000,nsave = 1000,
                        beta_prior_var=diag(ncol(X))*4,
                        beta_prior_mean=matrix(0,ncol(X),1),
                        family = "Poisson",
                        alpha_prop = "uniform",
                        silent = FALSE) {
  source("compute_mixture.R")
  if (family!="Poisson" && family != "NegBin") {
    stop("Family needs to be either Poisson or NegBin")
  }
  if (family == "NegBin" & (alpha_prop!="uniform" && alpha_prop!="truncnorm")) {
    stop("Alpha proposal densities has to be uniform or truncnorm")
  }
  
  n<-dim(X)[1]
  k<-dim(X)[2]
  
  curr.beta <- matrix(0, k, 1)
  betas <- array(0, c(nsave, k, 1))
  
  TI = solve(beta_prior_var)
  
  pred.Y <- array(0, c(nsave, n, 1))
  ptY <- array(0, c(nsave, n + sum(Y!=0), 1))
  sigi <- array(0, c(nsave, n + sum(Y!=0), 1))
  
  # add a small constant
  loglambda = log(matrix(0.1,n,1))
  loglambda[Y != 0] = log(Y[Y != 0])
  
  if (family == "NegBin") { 
    ## set up neg. bin random fx parameter
    nb_alpha = 10
    ## parameters for RW-MH
    c_nb = 10
    c_nb_adjust = 1.1
    nb_accept = 0
    alphas = matrix(0,nsave, 1)
    alpha_accept_rate = matrix(0,nburnin + nsave, 1)
  }
  
  # sort Y into ascending unique values
  # we'll have to get a mixture for each step
  ys = sort(Y)
  index = c(2:n)
  indexz = c(1,index[!duplicated(ys)[-1]])
  ny = ys[indexz]
  ny = ny[ny!=0] # ny contains the different unique steps
  
  Ny = sum(Y != 0) # number of mixture we have to store
  
  # storage structure for the first arrival step
  mixaux1=compute_mixture(1);
  #note: p = weight, m = mu, sigma = v, K = nc
  
  # the storage structure for all the normal mixtures correspoinding to
  # positive Ys
  Kmax = compute_mixture(ny[1])$nc
  mixaux2 = list(nc = matrix(0,Ny,1),
                 m = matrix(-1000,Ny,Kmax),
                 v = matrix(.001,Ny,Kmax),
                 p = matrix(0,Ny,Kmax))
  
 
  if (!silent) {
    cat("Pre-computing mixture..\n[")
    pb <- txtProgressBar(min = 0,max=length(ny),style = 3)
  }
  for (nj in 1:length(ny)) {
    # compute the corresponding mixture
    taux=compute_mixture(ny[nj]);
    Kaux = taux$nc
    
    # and store it in the same order as the positive values in Y
    # for each entry: p = nr of normals; m = corresponding weights
    #                 v = corr. sigmas; nc = corr. mus
    nm = sum(Y == ny[nj])
    mixaux2$p[Y[Y!=0]==ny[nj],1:Kaux]=kronecker(t(as.vector(taux$p)),matrix(1,nm,1));
    mixaux2$m[Y[Y!=0]==ny[nj],1:Kaux]=kronecker(t(as.vector(taux$m)),matrix(1,nm,1));
    mixaux2$v[Y[Y!=0]==ny[nj],1:Kaux]=kronecker(t(as.vector(taux$v)),matrix(1,nm,1));
    mixaux2$nc[Y[Y!=0]==ny[nj],1]=kronecker(Kaux,matrix(1,1,nm))
    
    if (!silent) {
      setTxtProgressBar(pb, nj)
    }
  }
  if (!silent) {
    close(pb)
    cat("Done!\n")
  }
  
  tX=rbind(X, as.matrix(X[Y!=0,]));
  data.bi = matrix(1,n,1)
  
  ttemp = rep(0,nburnin + nsave)
  
  #########################################################
  # START WITH MCMC
  if (!silent) {
    cat("Starting MCMC.. \n")
    cat("Burn-in\n") 
    pb <- txtProgressBar(min = 0,max=nburnin,style = 3)
  }
  for (mcmcit in 1:(nburnin + nsave)) {
    if (mcmcit != 1) {
      loglambda = X %*% curr.beta # calculate current loglambda
      
      ## if NB, calculate random fx
      if (family == "NegBin") {
        data.bi = rgamma(n,shape = nb_alpha + Y,rate = nb_alpha + exp(loglambda))
        loglambda = loglambda + log(data.bi)
      }
    }
    
    #### sample the auxiliary mixture
    # the variance parameter
    tSigma = matrix(0,n + sum(Y!=0),1)
    
    # the sampled responses  
    tY = matrix(0,n + sum(Y!=0),1)
    tY[(n+1):nrow(tY)] = rbeta(sum(Y!=0),Y[Y!=0],1) # tau2, last jump before 1
    tY[1:n]=1-log(runif(n)) / exp(loglambda);   # tau1, first jump after 1
    
    ttY = tY[1:n]
    ttY[Y!=0] = ttY[Y!=0] - tY[(n+1):length(tY)]
    tY[1:n] = ttY
    tY = -log(tY);
      
    # sample the INDICATORs for tau1: (arrival times)
    tempY = tY[1:n] - loglambda;
    ind = sample_mixture_ind(tempY,mixaux1,0);
    tY[1:n] =tY[1:n] - mixaux1$m[ind];
    
    tSigma[1:n]=mixaux1$v[ind];
      
    # sample the INDICATORs for tau2: (inter-arrival times)
    ny=sum(Y != 0);
    if (ny>0) {
      tempY = tY[(n+1):length(tY)] - loglambda[Y!=0];
      ind = sample_mixture_ind(as.matrix(tempY), mixaux2, 0);
      rep_ind = kronecker(ind,matrix(1,1,max(mixaux2$nc)))
      rep_K = kronecker(matrix(1:max(mixaux2$nc),1,max(mixaux2$nc)),matrix(1,ny,1))
      tY[(n+1):length(tY)] = tY[(n+1):length(tY)] - 
        rowSums( mixaux2$m * (rep_ind == rep_K  ) );
      tSigma[(n+1):length(tSigma)]= rowSums(mixaux2$v*(rep_ind==rep_K)  );
    }
    
    #### Sample for beta
    sigmainv = 1/tSigma^.5;  #squareroot
    yy=(tY*sigmainv); 
    XX=(tX * kronecker(sigmainv, matrix(1,1,ncol(tX))));
    post.B = solve(TI + t(XX) %*% XX);
    post.b  = post.B %*% (TI %*% beta_prior_mean + t(XX) %*% yy);
    curr.beta = mvrnorm(1,post.b,post.B)
    #### Beta sampling done
    
    #### Sample Negative Binomial parameter              ####
    #### done currently via Metropolis within Gibbs step ####
    #### with a truncated normal proposal distribution   ####
    #### or a uniform proposal distribution              ####
    if (family == "NegBin") {
      if (mcmcit == 1) {
        loglik = sum(llikeli_nb(Y,exp(X %*% curr.beta),1/nb_alpha))
      }
      if (alpha_prop == "uniform") {
        # uniform proposal would be an alternative
        Um = 1.1 * (2 * runif(1) -1)
        alpha_new = nb_alpha * exp(Um)
      } else if (alpha_prop == "truncnorm") {
        alpha_new = rtruncnorm(1,nb_alpha,c_nb,a = 0.001,b = 10000)
      } else {stop("Unknown alpha proposal density specified!")}
      if (alpha_new < .05) {alpha_new = .05} # too avoid numerical issues
      loglik_new =   sum(llikeli_nb(Y,exp(X %*% curr.beta),1/alpha_new))
      acc = loglik_new - loglik
      if (log(runif(1)) < acc) {
        nb_alpha = alpha_new
        loglik = loglik_new
        nb_accept = nb_accept + 1
      }
      alpha_accept_rate[mcmcit] = nb_accept / mcmcit
      if (mcmcit < nburnin / 2) {
       if (nb_accept / mcmcit < .3) {
         c_nb = c_nb / c_nb_adjust
       } else if (nb_accept / mcmcit > .6) {
         c_nb = c_nb * c_nb_adjust
       }
      }
      ttemp[mcmcit] = nb_alpha
    }
    #### Sampling NB par done ####
    
    if (mcmcit > nburnin) {
      ## store draws
      kp <- mcmcit - nburnin
      if (!silent) {
        if (kp == 1) {
          close(pb)
          cat("\nBurn-in done. Storing draws..\n")
          pb <- txtProgressBar(min = 0,max=nsave,style = 3)
        }
      }
      betas[kp, ,1] <- curr.beta
      
      if (family == "Poisson") {
        pred.Y[kp,,1] <- exp( X %*% curr.beta)
      } else if (family == "NegBin") {
        alphas[kp,1] = nb_alpha
        pred.Y[kp,,1] <- exp( X %*% curr.beta + log(data.bi))
      }
      
      ## get predictive draws
      ptY[kp,,1] <- tY
      sigi[kp,,1] = sigmainv
      if (!silent) {
        setTxtProgressBar(pb, mcmcit - nburnin )
      }
    } else {
      if (!silent) {
        setTxtProgressBar(pb, mcmcit)
      }
    }
  }
  if (!silent) {
    close(pb)
    cat("End of MCMC \n", date(), "\n")
  }
  
  output = list()
  output$betas = betas
  output$Y = pred.Y
  output$ptY = ptY
  output$sigi = sigi
  
  output$pbeta = apply(betas,c(2,3),mean)
  output$pY = apply(pred.Y,c(2,3),mean)
  if (family == "NegBin") {
    output$alphas = alphas
    output$palpha = mean(alphas)
  }
  
  if (family == "Poisson") {
    output$loglike = sum(llikeli_poisson(Y,output$pY))
  } else if (family == "NegBin") {
    output$loglike = sum(llikeli_nb(Y,output$pY,output$palpha))
  }
  output$origY = Y
  output$origX = X
  
  err = Y - output$pY
  output$RMSE = as.double(sqrt( ( t(err) %*% err  ) / n))
  output$BIC = -2 * output$loglike + (k+1) * log(n)
  
  # R/R-bar squared
  sige = apply(sigi,c(2,3),mean)
  ESS = t(err)  %*% err
  SStot = t(Y - mean(Y)) %*% (Y - mean(Y))
  output$R2 = 1 - ESS/SStot
  output$R2bar = output$R2 - (1-output$R2) * (k)/(n-k-1) 
  
  
  return(output)
}

sample_mixture_ind = function(y,mix,varargin) {
  n = nrow(y);
  
  nst=mix$nc;
  
  lh = mixture_llh(y,mix$m,mix$v)$lh;
  if (size(nst,1)*size(nst,2) != 1) { 
    nst=max(nst)
  }
  
  if  (size(mix$p,1)==1) {
    mix$p = kronecker(matrix(mix$p,1,nst),matrix(1,n,1))
  }
  
  p = mix$p *lh;
  sump = kronecker(rowSums(p),matrix(1,1,nst))
  p = p / sump;
  
  if (max(mix$nc) > 1) {
    rnd = kronecker( runif(n), matrix(1,1,max(mix$nc)));
    # sampling # sst 1 times n
    S = (rowSums( t(apply(p,c(1),cumsum)) < rnd,2) + 1);      
  } else {
    S = matrix(1,n,1);
  }  
  
  return(S)
}

llikeli_poisson = function(y,lmda) {
  n = nrow(y)
  nst = nrow(lmda)
  
  lmda[lmda<.0001] = .0001
  llh = y * log(lmda) - lmda - lgamma(y+1)
  return(llh)
}

llikeli_nb = function(y,lmda, nu ) {
  n = nrow(y)
  #lmda[lmda<.0001] = .0001
  llh = lgamma(nu + y) - lgamma(nu) - lgamma(y+1)
  llh = llh + nu*log(nu) + y*log(lmda) - (nu + y) * log(nu + lmda)
  return(llh)
}

mixture_llh = function(y,beta,sgma2) {
  n=length(y);
  nst=size(beta,2);
  
  if  (size(beta,1)==1) {
    beta = kronecker(matrix(beta,1,nst),matrix(1,n,1))
    sgma2 = kronecker(matrix(sgma2,1,nst),matrix(1,n,1))
  }    
  
  
  llh = -.5*(log(2*pi) + log(sgma2) + (y[,matrix(1,nst,1)] - beta)^2 / sgma2);
  
  
  maxl = apply(llh,c(1),max);
  temp_maxl = kronecker(matrix(maxl,n,1),matrix(1,1,nst))
  lh = exp(llh - temp_maxl);
  
  return(list(lh = lh,llh = llh,maxl = maxl))
}