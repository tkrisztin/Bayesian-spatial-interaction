### auxiliary mixture sampling for the Poisson distribution
### Based on Frühwirth Schnatter's bayesf MATLAB toolbox
### Implements spatial effects with normal-gamma prior like LeSage, Fischer, Scherngell (2007)
### modell is y ~ Po(lambda)  OR y ~ NB(lambda)
###  lambda = X b + V u_o + D u_d
###  u_o = r W u_o + e_o;  e_o = N(0,tau_o)
###  u_d = r W u_d + e_d;  e_d = N(0,tau_d)
### Y                   ...   dependent variable, N x 1 vector; vectorized n x n matrix
### X                   ...   independent variables, as a N x k matrix
### W                   ...   n x n row-standardized W matrix
### nburnin             ...   Number of burn in iterations in Gibbs sampler
### nsave               ...   Number of saved iterations in Gibbs sampler
### beta_prior_var      ...   Prior variance of beta, k x k matrix
### beta_prior_mean     ...   Prior mean of beta k x 1 vector
### ptau_o_s, ptau_o_v  ...   Prior on tau_o the spatial effects variance of origins
### ptau_d_s, ptau_d_v  ...   Prior on tau_d the spatial effects variance of destinations
### TT                  ...   optional, if data are in time panel format this is used to correctly
###                             build spatial effects; assumption: data are blockwise listed 1,..,TT
### family              ...   Can be NegBin or Poisson
### alpha_prop          ...   Proposal density for neg. bin. alpha parameter; 
###                             either uniform or truncnorm
source("lndetPaceBarry.R")
source("compute_mixture.R")
require(pracma)
require(MASS)
require(pscl)
require(akima)
require(truncnorm)
require(mgcv)

mixturemcmc_fx_fdi <- function(Y,X,W,include_flow,
                        nburnin = 1000,nsave = 1000,
                        beta_prior_var=diag(ncol(X) )*4,
                        beta_prior_mean=matrix(0,ncol(X),1),
                        ptau_o_s = 1/(4 * 10^8), ptau_o_v = 5, 
                        ptau_d_s = 1/(4 * 10^8), ptau_d_v = 5,
                        TT = 1,
                        family = "Poisson",
                        alpha_prop = "uniform") {
  if (family!="Poisson" && family != "NegBin") {
    stop("Family needs to be either Poisson or NegBin")
  }
  if (family == "NegBin" & (alpha_prop!="uniform" && alpha_prop!="truncnorm")) {
    stop("Alpha proposal densities has to be uniform or truncnorm")
  }
  
  smalln = nrow(W)
  n<-sum(include_flow)
  N = smalln^2
  
  Y = as.matrix(Y[include_flow,])
  X = X[include_flow,]
  
  k<-dim(X)[2]
  
  curr.beta <- matrix(0, k, 1)
  curr.rho_o = 0
  curr.rho_d = 0
  curr.tau_o = .1
  curr.tau_d = .1
  
  curr.uu_o = solve(diag(smalln) - curr.rho_o * W) %*% (rnorm(smalln) * (curr.tau_o))
  curr.uu_d = solve(diag(smalln) - curr.rho_d * W) %*% (rnorm(smalln) * (curr.tau_d))
  
  rho_os = matrix(1,nsave,1)
  rho_ds = matrix(1,nsave,1)
  tau_os = matrix(1,nsave,1)
  tau_ds = matrix(1,nsave,1)
  uu_os = matrix(1,smalln,nsave)
  uu_ds = matrix(1,smalln,nsave)
  betas <- array(0, c(nsave, k, 1))
  
  TI = solve(beta_prior_var)
  bTI = TI %*% beta_prior_mean
  
  #pred.Y <- array(0, c(nsave, n, 1))
  #ptY <- array(0, c(nsave, n + sum(Y!=0), 1))
  #sigi <- array(0, c(nsave, n + sum(Y!=0), 1))
  
  # spatial prior structure settings
  acc_rate_o = matrix(0,nburnin + nsave,1)
  acc_rate_d = matrix(0,nburnin + nsave,1)
  cc_o = 0.2
  cc_d = 0.2
  rmin = -1
  rmax = 1
  
  # pre-compute stuff for spatial
  detval = data.frame(sim_lndet(W,rmin,rmax))
  colnames(detval) = c("rho","detval")
  
  VVo = kronecker(matrix(1,TT,1),  kronecker(diag(smalln),matrix(1,smalln,1)) )
  VVd = kronecker(matrix(1,TT,1),  kronecker(matrix(1,smalln,1),diag(smalln)) )
  VVo = VVo[include_flow,]
  VVd = VVd[include_flow,]
  
  Ao = (diag(smalln) - curr.rho_o * W)
  Ad = (diag(smalln) - curr.rho_d * W)
  acc_o = 0
  acc_d = 0
  # end of spatial pre-comp
  
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
  
 
  cat("Pre-computing mixture..\n[")
  pb <- txtProgressBar(min = 0,max=length(ny),style = 3)
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
    
    setTxtProgressBar(pb, nj)
  }
  close(pb)
  cat("Done!\n")
  
  
  
  tX=rbind(X, X[Y!=0,]);
  tVVo = Matrix(rbind(VVo, VVo[Y!=0,]))
  tVVd = Matrix(rbind(VVd, VVd[Y!=0,]))
  
  data.bi = matrix(1,n,1)
  
  ttemp = rep(0,nburnin + nsave)
  
  #########################################################
  # START WITH MCMC
  cat("Starting MCMC.. \n")
  cat("Burn-in\n") 
  pb <- txtProgressBar(min = 0,max=nburnin,style = 3)
  start.time = Sys.time()
  for (mcmcit in 1:(nburnin + nsave)) {
    if (mcmcit != 1) {
      loglambda = as.matrix(
        X %*% curr.beta + VVo %*% curr.uu_o + VVd %*% curr.uu_d) # calculate current loglambda
      
      ## if NB, calculate random fx
      if (family == "NegBin") {
        cat("alpha:",nb_alpha,"\n")
        data.bi = rgamma(n,shape = nb_alpha + Y,rate =1 /( nb_alpha + exp(loglambda)))
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
     yy=as.matrix((tY - tVVo %*% curr.uu_o - tVVd %*% curr.uu_d )*sigmainv)
     XX=tX * c(sigmainv);
     V_post <- solve(TI + crossprod(XX) )
     beta_post <- V_post %*% (bTI + crossprod(XX,yy))
     curr.beta <- beta_post + t(chol(V_post))%*% rnorm(k,0,1)
     #### Beta sampling done
    
    # sample for uu_o and uu_d
    # uu_o
     sVVo = tVVo * c(sigmainv)
     sYY = (tY - tX %*% curr.beta - tVVd %*% curr.uu_d )*sigmainv
     PL = 1 / curr.tau_o * crossprod(Ao) + crossprod(sVVo)
     bL = crossprod(sVVo,sYY)
     curr.uu_o = chol_sample(bL,PL)
    
    # uu_d
     sVVd = tVVd * c(sigmainv)
     sYY = (tY - tX %*% curr.beta - tVVo %*% curr.uu_o )*sigmainv
     PL = 1 / curr.tau_d * crossprod(Ad) + crossprod(sVVd)
     bL = crossprod(sVVd,sYY)
     curr.uu_d = chol_sample(bL,PL)
     
    # sample for tau_o and tau_d
    # tau_o
    ttemp2 = as.double(crossprod(Ao %*% curr.uu_o))
    postv = smalln + ptau_o_v
    posts = (ttemp2 + ptau_o_v * ptau_o_s)/postv
    curr.tau_o = 1/rgamma(1,shape = postv/2, rate=.5*postv*posts)
   
    # # tau_d
    ttemp2 = as.double(crossprod(Ad %*% curr.uu_d))
    postv = smalln + ptau_d_v
    posts = (ttemp2 + ptau_d_v * ptau_d_s)/postv
    curr.tau_d = 1/rgamma(1,shape = postv/2, rate=.5*postv*posts)
    
    # # sample for rho_o
    curr.rho_o = griddy_sim(curr.uu_o,curr.tau_o,W,detval)
    Ao = (diag(smalln) - curr.rho_o * W)
    # rhox = as.double(c_sim(curr.rho_o,curr.uu_o,curr.tau_o,W,detval))
    # accept = 0
    # rho2 = curr.rho_o + cc_o * rnorm(1)
    # while (accept == 0) {
    #   if ((rho2>rmin) && (rho2 < rmax)) {
    #     accept = 1
    #   } else {
    #     rho2 = curr.rho_o + cc_o * rnorm(1)
    #   }
    # }
    # rhoy = as.double(c_sim(rho2,curr.uu_o,curr.tau_o,W,detval))
    # ru = runif(1)
    # if ((rhoy - rhox) > exp(1)) {
    #   p = 1
    # } else {
    #   ratio = exp(rhoy - rhox)
    #   p = min(1, ratio)
    # }
    # if (ru < p) {
    #   curr.rho_o = rho2
    #   Bo = (diag(smalln) - curr.rho_o * W)
    #   acc_o = acc_o + 1
    # }
    # acc_rate_o[mcmcit] = acc_o/mcmcit
    # # update cc based on std of rho draws
    # if (acc_rate_o[mcmcit] < .4) {
    #   cc_o = cc_o / 1.1
    # } else if (acc_rate_o[mcmcit] > .6) {
    #   cc_o = cc_o * 1.1
    # }
    
    # # sample for rho_d
    curr.rho_d = griddy_sim(curr.uu_d,curr.tau_d,W,detval)
    Ad = (diag(smalln) - curr.rho_d * W)
    # rhox = as.double(c_sim(curr.rho_d,curr.uu_d,curr.tau_d,W,detval))
    # accept = 0
    # rho2 = curr.rho_d + cc_d * rnorm(1)
    # while (accept == 0) {
    #   if ((rho2>rmin) && (rho2 < rmax)) {
    #     accept = 1
    #   } else {
    #     rho2 = curr.rho_d + cc_d * rnorm(1)
    #   }
    # }
    # rhoy = as.double(c_sim(rho2,curr.uu_d,curr.tau_d,W,detval))
    # ru = runif(1)
    # if ((rhoy - rhox) > exp(1)) {
    #   p = 1
    # } else {
    #   ratio = exp(rhoy - rhox)
    #   p = min(1, ratio)
    # }
    # if (ru < p) {
    #   curr.rho_d = rho2
    #   Bd = (diag(smalln) - curr.rho_d * W)
    #   acc_d = acc_d + 1
    # }
    # acc_rate_d[mcmcit] = acc_d/mcmcit
    # # update cc based on std of rho draws
    # if (acc_rate_d[mcmcit] < .4) {
    #   cc_d = cc_d / 1.1
    # } else if (acc_rate_d[mcmcit] > .6) {
    #   cc_d = cc_d * 1.1
    # }
    
    #### Sample Negative Binomial parameter              ####
    #### done currently via Metropolis within Gibbs step ####
    #### with a truncated normal proposal distribution   ####
    #### or a uniform proposal distribution              ####
    if (family == "NegBin") {
      loglambda = as.matrix(
        X %*% curr.beta + VVo %*% curr.uu_o + VVd %*% curr.uu_d) # calculate current loglambda
      loglik = sum(llikeli_nb(Y,exp(loglambda),1/nb_alpha))
      if (alpha_prop == "uniform") {
        # uniform proposal would be an alternative
        Um = 1.1 * (2 * runif(1) -1)
        alpha_new = nb_alpha * exp(Um)
      } else if (alpha_prop == "truncnorm") {
        alpha_new = rtruncnorm(1,nb_alpha,c_nb,a = 0.1,b = 10000)
      } else {stop("Unknown alpha proposal density specified!")}
      loglik_new =   sum(llikeli_nb(Y,exp(loglambda),1/alpha_new))
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
      if (kp == 1) {
        close(pb)
        cat("\nBurn-in done. Storing draws..\n")
        pb <- txtProgressBar(min = 0,max=nsave,style = 3)
      }
      betas[kp, ,1] <- curr.beta
      
      rho_os[kp] = curr.rho_o
      rho_ds[kp] = curr.rho_d
      tau_os[kp] = curr.tau_o
      tau_ds[kp] = curr.tau_d
      uu_os[,kp] = as.matrix(curr.uu_o)
      uu_ds[,kp] = as.matrix(curr.uu_d)
      
      # loglambda = as.matrix(X %*% curr.beta + VV %*% curr.uu_o + DD %*% curr.uu_d)
      # if (family == "Poisson") {
      #   pred.Y[kp,,1] <- exp( loglambda )
      # } else if (family == "NegBin") {
      #   alphas[kp,1] = nb_alpha
      #   pred.Y[kp,,1] <- exp( loglambda + log(data.bi))
      # }
      
      ## get predictive draws
      #ptY[kp,,1] <- tY
      #sigi[kp,,1] = sigmainv
      setTxtProgressBar(pb, mcmcit - nburnin )
    } else {
      setTxtProgressBar(pb, mcmcit)
    }
  }
  close(pb)
  cat("Time elapsed (sec):",as.numeric(Sys.time() - start.time,units = "secs"),"\n")
  cat("End of MCMC \n", date(), "\n")
  
  output = list()
  output$betas = betas
  #output$Y = pred.Y
  #output$ptY = ptY
  #output$sigi = sigi
  output$tau_os = tau_os
  output$tau_ds = tau_ds
  output$rho_os = rho_os
  output$rho_ds = rho_ds
  output$uu_os = uu_os
  output$uu_ds = uu_ds
  output$acc_rate_o = acc_rate_o
  output$acc_rate_d = acc_rate_d
  
  output$pbeta = apply(betas,c(2,3),mean)
  #output$pY = apply(pred.Y,c(2,3),mean)
  output$puu_o = apply(uu_os,c(1),mean)
  output$puu_d = apply(uu_ds,c(1),mean)
  output$pTau_o = mean(tau_os)
  output$pTau_d = mean(tau_ds)
  output$prho_o = mean(rho_os)
  output$prho_d = mean(rho_ds)
  
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
  
  # # R/R-bar squared
  # sige = apply(sigi,c(2,3),mean)
  # ESS = t(err)  %*% err
  # SStot = t(Y - mean(Y)) %*% (Y - mean(Y))
  # output$R2 = 1 - ESS/SStot
  # output$R2bar = output$R2 - (1-output$R2) * (k)/(n-k-1) 
  # 
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
  lmda[lmda<.0001] = .0001
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

sim_lndet = function(W,rmin,rmax,order = 50,iter = 30) {
  # PURPOSE: compute the log determinant |I_n - rho*W|
  # using the user-selected (or default) method
  # ---------------------------------------------------
  #  USAGE: detval = far_lndet(lflag,W,rmin,rmax)
  # where eflag,rmin,rmax,W contains input flags 
  # and the outputs are either user-inputs or default values
  # ---------------------------------------------------
  
  # use Pace and Barry, 1999 MC approximation  
  out = lndetPaceBarry(order,iter,W,rmin,rmax); 
  tt=seq(rmin,rmax,.001) # interpolate a finer grid
  outi = interp1(out[,2],out[,1],tt,'spline');
  detval = cbind(tt,outi)
  
  return(detval)
}



chol_sample = function(bL,PL) {
  nn = nrow(bL)
  ## Can speed up using Choleksy.
  V1 = chol2inv(chol(as.matrix(PL)));
  m1 = V1 %*% (bL );
  sqrtV1 = t(chol(V1));
  return(m1 + sqrtV1 %*% rnorm(nn))
}

c_sim = function(rho,bb,sige,W,detval) {
  index_r = which(abs(detval$rho-rho)==min(abs(detval$rho-rho)))
  detm = detval$detval[index_r]
  
  z = bb - rho*W%*%bb;
  epe = (t(z)%*%z)/(2*sige);
  
  return(detm - epe)
}

griddy_sim = function(bb,sige,W,detval) {
  rrhos = detval$rho
  Wb = W %*% bb
  bpb = as.double(crossprod(bb))
  WpW = as.double(crossprod(Wb))
  bpW = as.double(t(bb) %*% Wb)
  z = bpb  - 2 * rrhos * bpW + rrhos^2 * WpW
  z = z/(2*sige)
  den = detval$detval - z 
  y = rrhos
  adj = max(den)
  den = den - adj
  x = exp(den)
  isum = sum((y[-1] + y[-length(y)])*(x[-1]  - x[-length(x)])/2)
  z = abs(x/isum)
  den = cumsum(z)
  rnd = runif(1) * sum(z)
  ind = max(c(1,which(den <= rnd)))
  
  return(rrhos[ind])
}


gen_d = function(n,r,p,c) {
  mub = (r*c)/(1-r)
  bbb = rbeta(n,mub + p,c + p)
  return(bbb)
}