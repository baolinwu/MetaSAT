## internal func.
# ' Adaptive variant-set association test
# '
# ' Adaptive test (AT) is defined based on the minimum p-value of weighted averages of generalized burden test (GBT) and sum of squares test (SST).
# ' The required inputs are: U (vector of test statistics, say Score vector); R (covariance matrix of U); eta (GBT coefficient);
# ' rho (vector of weights assigned to the GBT).
# ' SST: \eqn{U^TU}; GBT: \eqn{(\eta^TU)^2}; minimum p-value over weighted tests: \eqn{(1-\rho)U^TU + \rho(\eta^TU)^2}.
# ' 
# ' @param  U    vector of variant test statistics
# ' @param  R    covariance matrix for test statistics
# ' @param  eta  coefficient vector for the GBT. Default to equal weights.
# ' @param  rho  weights for the GBT
# ' @return 
# ' \describe{
# '   \item{p.value}{ p-value for AT }
# '   \item{pval}{ vector of p-values for each rho }
# '   \item{rho.est}{ optimal rho value leading to the minimum p-value }
# ' }
# ' @keywords adaptive test
# ' @export
# ' @references
# ' Guo,B., Massoti,M., Liu,N. and Wu,B. (2018).  Efficient and powerful meta-analysis of variant-set association using MetaSAT.
# ' @examples
# ' R = cor(matrix(rnorm(500),100,5)*sqrt(0.9)+rnorm(100)*sqrt(0.1))
# ' Rh = chol(R)
# ' Z = colSums(Rh*rnorm(5))
# ' AVAT(Z, R)
AVAT <- function(U,R, eta=NULL, rho=c((0:5/10)^2,0.5,1)){
  m = length(U)
  if(is.null(eta)) eta = rep(1,m)
  ##
  eR = eigen(R, sym=TRUE); D = sqrt(abs(eR$val))
  U1 = colSums(eR$vec*eta)*D
  Rs = colSums(R*eta); R1 = sum(Rs*eta); R2 = sum(Rs^2)
  P1 = diag(D^2); P2 = outer(U1,U1)
  ## min-pval
  Q = sum(U^2); B = sum(U*eta)
  Qs = (1-rho)*Q + rho*B^2
  K = length(rho)
  Lamk = vector('list', K)
  pval = rep(0,K)
  for(k in 1:K){
    ak = (1-rho[k])*P1 + rho[k]*P2
    aae = zapsmall( abs( eigen(ak,sym=TRUE,only.val=TRUE)$val ) )
    Lamk[[k]] = aae[aae>0]
    pval[k] = KATpval(Qs[k], Lamk[[k]])
  }
  minp = min(pval)
  ## min-pval
  K = length(rho); K1 = K-1
  qval = rep(0,K1)
  for(k in 1:K1) qval[k] = Liu0.qval(minp, Lamk[[k]])
  V1 = P1 - P2*R2/R1^2
  Lame = eigen(V1, sym=TRUE, only.val=TRUE)$val
  ##
  rho1 = rho[-K]
  tau1 = (1-rho1)*R2/R1^2+rho1
  fint = function(xu){
    sapply(xu, function(x){
      x1 = min( (qval-tau1*x^2*R1)/(1-rho1) )
      p1 = KATpval(x1,Lame)
      p1*dnorm(x)
    })
  }
  qkh = -qnorm(minp/2); prec = 1e-5
  p.value = try({ minp + 2*integrate(fint, 0,qkh,  subdivisions=1e3,abs.tol=minp*prec)$val }, silent=TRUE)
  while(class(p.value)=='try-error'){
    prec = prec*2
    p.value = try({ minp + 2*integrate(fint, 0,qkh, abs.tol=minp*prec)$val }, silent=TRUE)
  }
  return( list(p.value=p.value, pval=pval, rho.est=rho[which.min(pval)]) )
}

#' Fixed-effects meta-analysis of variant-set association
#'
#' Conduct meta-analysis of variant-set association test of m variants assuming constant effects across K studies.
#' These association statistics are typically Score vector, and direct summation asymptotically amounts to inverse
#' variance weighting. In practice, we can always input weighted test statistics. 
#' SST: \eqn{Q=(\sum_kU_k)^T(\sum_kU_k)}; GBT: \eqn{(\sum_k\eta^TU_k)^2}; AT: adaptively weighting SST and GBT. 
#' 
#' @param  Us    matrix of variant test statistics (m by K)
#' @param  Rs    array of covariance matrix for test statistics (mxm by K)
#' @param  eta   coefficient vector for the GBT. Default to equal weights.
#' @param  rho   weights assigned to the GBT
#' @return 
#' \describe{
#'   \item{p.value}{ p-value for AT }
#'   \item{pval}{ vector of p-values for each rho }
#'   \item{rho.est}{ optimal rho value leading to the minimum p-value }
#' }
#' @keywords AT
#' @export
#' @references
#' Guo,B., Massoti,M., Liu,N. and Wu,B. (2018). Efficient and powerful meta-analysis of variant-set association using MetaSAT.
#' @examples
#' K = 3; m=10
#' Rs = array(0, dim=c(m,m,K)); Us = matrix(0, m,K)
#' for(k in 1:K){
#'   ak = matrix(rnorm(100*m),100,m)*sqrt(0.8)+rnorm(100)*sqrt(0.2)
#'   Rs[,,k] = cor(ak)
#'   Rh = chol(Rs[,,k])
#'   Us[,k] = colSums(Rh*rnorm(m))
#' }
#' FMSAT(Us,Rs)
#' U1 = Us + runif(m*K, 0,2)
#' FMSAT(U1,Rs)
FMSAT <- function(Us,Rs, eta=NULL, rho=c((0:5/10)^2,0.5,1)){
  ## summary
  U = rowSums(Us)
  R = apply(Rs, 1:2, sum)
  ## test
  ans = AVAT(U,R, eta, rho)
  return(ans)
}


#' Heterogeneous-effects meta-analysis of variant-set association
#'
#' Conduct meta-analysis of variant-set association test of m variants assuming heterogeous effects across K studies.
#' SST: \eqn{Q=\sum_k U_k^TU_k}; GBT: \eqn{(\sum_k\eta_k^TU_k)^2}; AT: adaptively weighting SST and GBT. 
#' 
#' @param  Us    matrix of variant test statistics (m, K)
#' @param  Rs    array of covariance matrix for test statistics (m, m, K)
#' @param  eta   coefficient vector for the GBT. Default to equal weights.
#' @param  rho   weights assigned to the GBT
#' @return 
#' \describe{
#'   \item{p.value}{ p-value for the adaptive test (AT) }
#'   \item{pval}{ vector of p-values for each rho }
#'   \item{rho.est}{ optimal rho value leading to the minimum p-value }
#' }
#' @keywords SAT
#' @export
#' @references
#' Guo,B., Massoti,M., Liu,N. and Wu,B. (2018). Efficient and powerful meta-analysis of variant-set association using MetaSAT.
#' @examples
#' K = 3; m=10
#' Rs = array(0, dim=c(m,m,K)); Us = matrix(0, m,K)
#' for(k in 1:K){
#'   ak = matrix(rnorm(100*m),100,m)*sqrt(0.8)+rnorm(100)*sqrt(0.2)
#'   Rs[,,k] = cor(ak)
#'   Rh = chol(Rs[,,k])
#'   Us[,k] = colSums(Rh*rnorm(m))
#' }
#' HMSAT(Us,Rs)
#' U1 = Us + rnorm(m*K)
#' HMSAT(U1,Rs)
HMSAT <- function(Us,Rs, eta=NULL, rho=c((0:5/10)^2,0.5,1)){
  m = dim(Us)[1]; K = dim(Us)[2]; mK=m*K
  if(is.null(eta)){
    etah = rep(1,mK)
  } else{
    etah = rep(eta, K)[1:mK]
  }
  ## summary
  Uh = as.vector(Us)
  Vh = matrix(0, mK, mK)
  for(k in 1:K){
    ik = (k-1)*m + 1:m
    Vh[ik,ik] = Rs[,,k]
  }
  ## test
  ans = AVAT(Uh,Vh, etah, rho)
  return(ans)
}


#' Robust heterogeneous-effects meta-analysis of variant-set association 
#'
#' Conduct meta-analysis of variant-set association test of m variants assuming heterogeous effects across K studies.
#' SST: \eqn{Q=\sum_k U_k^TU_k}; SBT: \eqn{\sum_k(\eta_k^TU_k)^2} 
#' 
#' @param  Us    matrix of variant test statistics (m by K)
#' @param  Rs    array of covariance matrix for test statistics (m,m by K)
#' @param  eta   coefficient vector for the SBT. Default to equal weights.
#' @param  rho   weights assigned to the SBT
#' @return 
#' \describe{
#'   \item{p.value}{ p-value for AT }
#'   \item{pval}{ vector of p-values for each rho }
#'   \item{rho.est}{ optimal rho value leading to the minimum p-value }
#' }
#' @keywords SAT
#' @export
#' @references
#' Guo,B., Massoti,M., Liu,N. and Wu,B. (2018). Efficient and powerful meta-analysis of variant-set association using MetaSAT.
#' @examples
#' K = 3; m=10
#' Rs = array(0, dim=c(m,m,K)); Us = matrix(0, m,K)
#' for(k in 1:K){
#'   ak = matrix(rnorm(100*m),100,m)*sqrt(0.8)+rnorm(100)*sqrt(0.2)
#'   Rs[,,k] = cor(ak)
#'   Rh = chol(Rs[,,k])
#'   Us[,k] = colSums(Rh*rnorm(m))
#' }
#' RMSAT(Us,Rs)
#' U1 = Us + rnorm(m*K,1,1.25)
#' RMSAT(U1,Rs)
RMSAT <- function(Us,Rs, eta=NULL, rho=c((0:5/10)^2,0.5,1)){
  m = dim(Us)[1]; K = dim(Us)[2]; mK=m*K
  if(is.null(eta)){
    etas = matrix(1, m,K)
  } else{
    etas = matrix(rep(eta, K)[1:mK], m,K)
  }
  R1 = R2 = rep(0, K)
  for(k in 1:K){
    a1 = colSums(Rs[,,k]*etas[,k])
    R1[k] = sum(a1*etas[,k])
    R2[k] = sum(a1^2)
  }
  ## 
  Q = sum(Us^2); B = colSums(Us*etas)
  Qs = (1-rho)*Q + rho*sum(R2/R1^2*B^2)
  L = length(rho); rho1 = rho[-L]
  Lamk = vector('list', L)
  Lame = NULL
  for(k in 1:K){
    Rk = Rs[,,k]; eta = etas[,k]
    eR = eigen(Rk, sym=TRUE); D = sqrt(abs(eR$val))
    U1 = colSums(eR$vec*eta)*D
    P1 = diag(D^2); P2 = outer(U1,U1)
    for(j in 1:L){
      ak = (1-rho[j])*P1 + rho[j]*R2[k]/R1[k]^2*P2
      aae = zapsmall( abs( eigen(ak,sym=TRUE,only.val=TRUE)$val ) )
      Lamk[[j]] = c(Lamk[[j]], aae[aae>0])
    }
    V1 = P1 - P2*R2[k]/R1[k]^2
    Lame = c(Lame, eigen(V1, sym=TRUE, only.val=TRUE)$val )
  }
  pval = rep(1,L)
  for(j in 1:L) pval[j] = KATpval(Qs[j], Lamk[[j]])
  minp = min(pval)
  L1 = L-1
  qval = rep(0,L1)
  for(j in 1:L1) qval[j] = Liu0.qval(minp, Lamk[[j]])
  ##
  Lamb = R2/R1;  q2 = KATqval(minp, Lamb)
  B = 1e3
  q2x = seq(0, q2, length=B)
  p2x = KATpval(q2x, Lamb)
  q1x = sapply(q2x, function(x) min( (qval-x)/(1-rho1) ) )
  p1x = KATpval(q1x, Lame)
  p.val = minp - sum( diff(p2x)*(p1x[-1]+p1x[-B])/2 )
  return( list(p.value=p.val, pval=pval, rho.est=rho[which.min(pval)]) )
}

