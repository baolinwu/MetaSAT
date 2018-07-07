## Internal func.
# ' Adaptive variant-set association test
# '
# ' Adaptive test (AT) is defined based on the minimum p-value of weighted averages of generalized burden test (BT) and sum of squares test (VT).
# ' The required inputs are: U (vector of test statistics, say Score vector); R (covariance matrix of U); eta (BT coefficient);
# ' rho (vector of weights assigned to the BT).
# ' VT: \eqn{U^TU}; BT: \eqn{(\eta^TU)^2}; minimum p-value over weighted tests: \eqn{(1-\rho)U^TU + \rho(\eta^TU)^2}.
# '
# ' @param  U    vector of variant test statistics
# ' @param  R    covariance matrix for test statistics
# ' @param  eta  coefficient vector for the BT. Default to equal weights.
# ' @param  rho  weights for the BT
# ' @return
# ' \describe{
# '   \item{p.value}{ p-values for AT,VT,BT }
# '   \item{pval}{ vector of p-values for each rho }
# '   \item{rho.est}{ optimal rho value leading to the minimum p-value }
# ' }
# ' @keywords adaptive test
# ' @export
# ' @references
# ' Wu,B. and Zhao,H. (2018). Efficient and powerful meta-analysis of variant-set association tests using MetaSAT.
# ' @examples
# ' R = cor(matrix(rnorm(500),100,5)*sqrt(0.9)+rnorm(100)*sqrt(0.1))
# ' Rh = chol(R)
# ' Z = colSums(Rh*rnorm(5))
# ' AVAT(Z, R)
AVAT <- function(U,R, eta=NULL, rho=(0:10/10)^2){
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
  pval = rep(1,K)
  for(k in 1:K){
    ak = (1-rho[k])*P1 + rho[k]*P2
    aae = zapsmall( abs( eigen(ak,sym=TRUE,only.val=TRUE)$val ) )
    Lamk[[k]] = aae[aae>0]
    pval[k] = KATpval(Qs[k], Lamk[[k]], acc=1e-3)
  }
  minp = min(pval)
  ## sim
  ## if( (minp<1e-9)|(minp>1e-8) )     return( list(p.value=c(minp*2.5,pval[c(1,K)]), pval=pval, rho.est=rho[which.min(pval)]) )
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
      p1 = KATpval(x1,Lame, acc=1e-3)
      p1*dnorm(x)
    })
  }
  qkh = -qnorm(minp/2); prec = 1e-5
  p.value = try({ minp + 2*integrate(fint, 0,qkh,  subdivisions=1e3,abs.tol=minp*prec)$val }, silent=TRUE)
  while(class(p.value)=='try-error'){
    prec = prec*2
    p.value = try({ minp + 2*integrate(fint, 0,qkh, abs.tol=minp*prec)$val }, silent=TRUE)
  }
  ans = c(p.value, pval[rho==0], pval[rho==1])
  names(ans) = c('AT', 'VT', 'BT')
  return( list(p.value=ans, pval=pval, rho.est=rho[which.min(pval)]) )
}

#' Fixed-effects meta-analysis of variant-set association
#'
#' Conduct meta-analysis of variant-set association test of m variants assuming constant effects across K studies.
#' These association statistics are typically Score vector, and a direct summation asymptotically amounts to inverse
#' variance weighting. In practice, we typically input weighted test statistics.
#' FE VT: \eqn{(\sum_kU_k)^T(\sum_kU_k)}; FE BT: \eqn{(\sum_k\eta^TU_k)^2}; FAT: adaptively combine FE VT and FE BT.
#'
#' @param  Us    matrix of variant test statistics (m by K)
#' @param  Vs    array of covariance matrix for test statistics (mxm by K)
#' @param  eta   coefficient vector for the FE BT. Default to equal weights.
#' @param  rho   weights assigned to the FE BT
#' @return
#' \describe{
#'   \item{p.value}{ p-values for FAT, FE VT, FE BT }
#'   \item{pval}{ vector of p-values for each rho }
#'   \item{rho.est}{ optimal rho value leading to the minimum p-value }
#' }
#' @keywords FAT
#' @export
#' @references
#' Wu,B. and Zhao,H. (2018). Efficient and powerful meta-analysis of variant-set association tests using MetaSAT.
#' @examples
#' K = 3; m=10
#' Vs = array(0, dim=c(m,m,K)); Us = matrix(0, m,K)
#' for(k in 1:K){
#'   ak = matrix(rnorm(100*m),100,m)*sqrt(0.8)+rnorm(100)*sqrt(0.2)
#'   Vs[,,k] = cor(ak)
#'   Rh = chol(Vs[,,k])
#'   Us[,k] = colSums(Rh*rnorm(m))
#' }
#' FMSAT(Us,Vs)
#' U1 = Us + runif(m*K, 0,2)
#' FMSAT(U1,Vs)
FMSAT <- function(Us,Vs, eta=NULL, rho=(0:10/10)^2){
  ## summary
  U = rowSums(Us)
  R = apply(Vs, 1:2, sum)
  ## test
  ans = AVAT(U,R, eta, rho)
  names(ans$p.value) = c('FAT', 'FE VT', 'FE BT')
  return(ans)
}


#' Heterogeneous-effects meta-analysis of variant-set association
#'
#' Conduct meta-analysis of variant-set association test of m variants assuming heterogeous effects across K studies.
#' HE VT: \eqn{\sum_k U_k^TU_k}; FE BT: \eqn{(\sum_k\eta_k^TU_k)^2}; HAT: adaptively combine HE VT and FE BT.
#'
#' @param  Us    matrix of variant test statistics (m, K)
#' @param  Vs    array of covariance matrix for test statistics (m, m, K)
#' @param  eta   coefficient vector for the FE BT. Default to equal weights.
#' @param  rho   weights assigned to the FE BT
#' @return
#' \describe{
#'   \item{p.value}{ p-values for HAT, HE VT, FE BT }
#'   \item{pval}{ vector of p-values for each rho }
#'   \item{rho.est}{ optimal rho value leading to the minimum p-value }
#' }
#' @keywords HAT
#' @export
#' @references
#' Wu,B. and Zhao,H. (2018). Efficient and powerful meta-analysis of variant-set association tests using MetaSAT.
#' @examples
#' K = 3; m=10
#' Vs = array(0, dim=c(m,m,K)); Us = matrix(0, m,K)
#' for(k in 1:K){
#'   ak = matrix(rnorm(100*m),100,m)*sqrt(0.8)+rnorm(100)*sqrt(0.2)
#'   Vs[,,k] = cor(ak)
#'   Rh = chol(Vs[,,k])
#'   Us[,k] = colSums(Rh*rnorm(m))
#' }
#' HMSAT(Us,Vs)
#' U1 = Us + rnorm(m*K)
#' HMSAT(U1,Vs)
HMSAT <- function(Us,Vs, eta=NULL, rho=(0:10/10)^2){
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
    Vh[ik,ik] = Vs[,,k]
  }
  ## test
  ans = AVAT(Uh,Vh, etah, rho)
  names(ans$p.value) = c('HAT', 'HE VT', 'FE BT')
  return(ans)
}


#' Robust heterogeneous-effects meta-analysis of variant-set association
#'
#' Conduct meta-analysis of variant-set association test of m variants assuming heterogeous effects across K studies.
#' HE VT: \eqn{\sum_k U_k^TU_k}; RHE BT: \eqn{\sum_k(\eta_k^TU_k)^2}; RAT: adaptively combine HE VT and RHE BT.
#'
#' @param  Us    matrix of variant test statistics (m by K)
#' @param  Vs    array of covariance matrix for test statistics (m,m by K)
#' @param  eta   coefficient matrix for variants (m by K). Default to equal weights.
#' @param  rho   weights assigned to the RHE BT
#' @return
#' \describe{
#'   \item{p.value}{ p-values for RAT, HE VT, RHE BT }
#'   \item{pval}{ vector of p-values for each rho }
#'   \item{rho.est}{ optimal rho value leading to the minimum p-value }
#' }
#' @keywords RAT
#' @export
#' @references
#' Wu,B. and Zhao,H. (2018). Efficient and powerful meta-analysis of variant-set association tests using MetaSAT.
#' @examples
#' K = 3; m=10
#' Vs = array(0, dim=c(m,m,K)); Us = matrix(0, m,K)
#' for(k in 1:K){
#'   ak = matrix(rnorm(100*m),100,m)*sqrt(0.8)+rnorm(100)*sqrt(0.2)
#'   Vs[,,k] = cor(ak)
#'   Rh = chol(Vs[,,k])
#'   Us[,k] = colSums(Rh*rnorm(m))
#' }
#' RMSAT(Us,Vs)
#' U1 = Us + rnorm(m*K,1,1.25)
#' RMSAT(U1,R=Vs)
RMSAT <- function(Us,Vs, eta=NULL, rho=(0:10/10)^2){
  m = dim(Us)[1]; K = dim(Us)[2]; mK=m*K
  if(is.null(eta)){
    etas = matrix(1, m,K)
  } else{
    etas = matrix(rep(eta, K)[1:mK], m,K)
  }
  R1 = R2 = rep(0, K)
  for(k in 1:K){
    a1 = colSums(Vs[,,k]*etas[,k])
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
    Rk = Vs[,,k]; eta = etas[,k]
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
  for(j in 1:L) pval[j] = KATpval(Qs[j], Lamk[[j]], acc=1e-3)
  minp = min(pval)
  ## sim
  ##  if( (minp<1e-9)|(minp>1e-8) )    return( list(p.value=c(minp*2.5,pval[c(1,L)]), pval=pval, rho.est=rho[which.min(pval)]) )
  ##
  L1 = L-1
  qval = rep(0,L1)
  for(j in 1:L1) qval[j] = Liu0.qval(minp, Lamk[[j]])
  ##
  Lamb = R2/R1;  q2 = KATqval(minp, Lamb)
  B = 1e3
  q2x = seq(0, q2, length=B)
  p2x = KATpval(q2x, Lamb, acc=1e-3)
  q1x = sapply(q2x, function(x) min( (qval-x)/(1-rho1) ) )
  p1x = KATpval(q1x, Lame, acc=1e-3)
  p.val = minp - sum( diff(p2x)*(p1x[-1]+p1x[-B])/2 )
  ans = c(p.val, pval[rho==0], pval[rho==1])
  names(ans) = c('RAT', 'HE VT', 'RHE BT')
  return( list(p.value=ans, pval=pval, rho.est=rho[which.min(pval)]) )
}

#' Adaptive variant-set association test based on FE and RHE meta-analysis models
#'
#' Conduct meta-analysis of variant-set association test of m variants assuming similar effects across variants.
#' RHE BT: \eqn{\sum_k(\eta_k^TU_k)^2}; FE BT: \eqn{(\sum_k\eta_k^TU_k)^2}; BAT: adaptively combine RHE BT and FE BT.
#'
#' @param  Us    matrix of variant test statistics (m by K)
#' @param  Vs    array of covariance matrix for test statistics (m,m by K)
#' @param  eta   coefficient matrix for variants (m by K). Default to equal weights.
#' @param  rho   weights assigned to the FE BT
#' @return
#' \describe{
#'   \item{p.value}{ p-values for BAT, RHE BT, FE BT }
#'   \item{pval}{ vector of p-values for each rho }
#'   \item{rho.est}{ optimal rho value leading to the minimum p-value }
#' }
#' @keywords BAT
#' @export
#' @references
#' Wu,B. and Zhao,H. (2018). Efficient and powerful meta-analysis of variant-set association tests using MetaSAT.
#' @examples
#' K = 3; m=10
#' Vs = array(0, dim=c(m,m,K)); Us = matrix(0, m,K)
#' for(k in 1:K){
#'   ak = matrix(rnorm(100*m),100,m)*sqrt(0.8)+rnorm(100)*sqrt(0.2)
#'   Vs[,,k] = cor(ak)
#'   Rh = chol(Vs[,,k])
#'   Us[,k] = colSums(Rh*rnorm(m))
#' }
#' RBAT(Us,Vs)
#' U1 = Us + rnorm(m*K,1,1.25)
#' RBAT(U1,Vs)
RBAT <- function(Us,Vs, eta=NULL, rho=(0:10/10)^2){
  m = dim(Us)[1]; K = dim(Us)[2]; mK=m*K
  if(is.null(eta)){
    etas = matrix(1, m,K)
  } else{
    etas = matrix(rep(eta, K)[1:mK], m,K)
  }
  ## summary
  R1 = R2 = rep(0, K)
  for(k in 1:K){
    a1 = colSums(Vs[,,k]*etas[,k])
    R1[k] = sum(a1*etas[,k])
    R2[k] = sum(a1^2)
  }
  ## BAT
  U = colSums(Us*etas)*sqrt(R2)/R1
  R = R2/R1
  eta = R1/sqrt(R2)
  ans = AVAT(U,diag(R),eta, rho=rho)
  names(ans$p.value) = c('BAT', 'RHE BT', 'FE BT')
  return(ans)
}

## Reference implementation
RBAT0 <- function(Us,Vs, eta=NULL, rho=(0:10/10)^2){
  m = dim(Us)[1]; K = dim(Us)[2]; mK=m*K
  if(is.null(eta)){
    etas = matrix(1, m,K)
  } else{
    etas = matrix(rep(eta, K)[1:mK], m,K)
  }
  ## summary
  R1 = rep(0, K)
  for(k in 1:K){
    a1 = colSums(Vs[,,k]*etas[,k])
    R1[k] = sum(a1*etas[,k])
  }
  U = colSums(Us*etas)
  ## test
  ans = AVAT(U,diag(R1), rho=rho)
  names(ans$p.value) = c('BAT0', 'RHE BT', 'FE BT')
  return(ans)
}
