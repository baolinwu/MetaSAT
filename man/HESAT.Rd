% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/metasat.R
\name{HESAT}
\alias{HESAT}
\title{Heterogeneous-effects (HE) meta-analysis of variant-set association}
\usage{
HESAT(Us, Vs, eta = NULL, rho = (0:10/10)^2)
}
\arguments{
\item{Us}{matrix of variant test statistics (m, K)}

\item{Vs}{array of covariance matrix for test statistics (m, m, K)}

\item{eta}{coefficient vector for the FE BT. Default to equal weights.}

\item{rho}{weights assigned to the FE BT}
}
\value{
\describe{
  \item{p.value}{ p-values for HAT, HE VT, FE BT }
  \item{pval}{ vector of p-values for each rho }
  \item{rho.est}{ optimal rho value leading to the minimum p-value }
}
}
\description{
Conduct meta-analysis of variant-set association test of m variants assuming heterogeous effects across K studies.
HE VT: \eqn{\sum_k U_k^TU_k}; FE BT: \eqn{(\sum_k\eta_k^TU_k)^2}; HAT: adaptively combine HE VT and FE BT.
}
\examples{
K = 3; m=10
Vs = array(0, dim=c(m,m,K)); Us = matrix(0, m,K)
for(k in 1:K){
  ak = matrix(rnorm(100*m),100,m)*sqrt(0.8)+rnorm(100)*sqrt(0.2)
  Vs[,,k] = cor(ak)
  Rh = chol(Vs[,,k])
  Us[,k] = colSums(Rh*rnorm(m))
}
HESAT(Us,Vs)
U1 = Us + rnorm(m*K)
HESAT(U1,Vs)
}
\references{
Wu,B. and Zhao,H. (2018) Efficient and powerful meta-analysis of variant-set association tests using MetaSAT.
}
\keyword{HAT}
