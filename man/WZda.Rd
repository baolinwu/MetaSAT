% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{WZda}
\alias{WZda}
\title{Summary score statistics and associated asymptotic covariance matrices for nine rare variants across three studies}
\format{A list with two components
\describe{
  \item{Us}{9x3 matrix of summary score statistics (9 variants across 3 studies)}
  \item{Vs}{9x9x3 array of covariance matrix (9 variants across 3 studies) }
}}
\usage{
WZda
}
\description{
A list containing the summary statistics matrix (9x3) and the covariance matrix array (9x9x3)
}
\examples{
data(WZda)
Us = WZda$Us; Vs = WZda$Vs
FESAT(Us,Vs)
HESAT(Us,Vs)
RESAT(Us,Vs)
RBAT(Us,Vs)
}
\references{
Wu,B. and Zhao,H. (2018) Efficient and powerful meta-analysis of variant-set association tests using MetaSAT.

Wu,B., Guan,W. and Pankow,J.S. (2016) On Efficient and Accurate Calculation of Significance P-Values for Sequence Kernel Association Testing of Variant Set. Ann. Hum. Genet., 80 (2), 123–135.

Wessel,J., Chu,A.Y., Willems,S.M. and others. (2015) Low-frequency and rare exome chip variants associate with fasting glucose and type 2 diabetes susceptibility. Nature Communications, 6, 5897.
}
\keyword{datasets}
