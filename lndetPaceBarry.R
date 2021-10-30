require(matrixcalc)

lndetPaceBarry <- function(order=50, iter=30, W, rmin=-1, rmax = 1,qq = 0.01){
  n=dim(W)[1]
  
  # Exact moments from 1 to oexact
  td=matrix(c(0,sum(W^2)/2),length(c(0,sum(W^2)/2)),1)

  oexact=length(td)
  
  # stochastic moments
  mavmomi=matrix(0,order,iter)
  
  for(j in 1:iter)
  {
    u=matrix(rnorm(n,0,1),n,1)
    v=u
    utu=t(u)%*%u
    for (i in 1:order)
    {
      v=W%*%v
      mavmomi[i,j]=as.double(n*((t(u)%*%v)/(i*utu)))
      
    }  
  }
    mavmomi[1:oexact,]=td[,matrix(1,iter,1)]
    
  # averages across iterations
  avmomi=as.matrix(rowMeans(mavmomi))
  
  # alpha matrix
  alpha=seq(rmin,rmax,qq)
  valpha=vandermonde.matrix(alpha,length(alpha))
  alomat=-valpha[,(2:(order+1))]
  
  # Estimated ln|I-aD| using mixture of exact, stochastic moments
  # exact from 1 to oexact, stochastic from (oexact+1) to order
  
  lndetmat=alomat%*%avmomi
  
  
  return(cbind(lndetmat,alpha))
}