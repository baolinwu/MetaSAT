#' Summary score statistics of 9 rare variants across 3 studies and associated asymptotic covariance matrices
#'
#' A list containing the matrix of summary statistics (9x3) and covariance matrices (9x9x3)
#'
#' @format A list with two components
#' \describe{
#'   \item{Us}{9 by 3 matrix of summary score statistics (9 variants across 3 studies)}
#'   \item{Rs}{9x9x3 array of covariance matrix (9 variants across 3 studies) }
#' }
#' @source \url{http://www.github.com/baolinwu/MetaSAT}
#' @examples
#' data(MAData)
#' Us = MAData$Us; Rs = MAData$Rs
#' FMSAT(Us,Rs)
#' HMSAT(Us,Rs)
#' RMSAT(Us,Rs)
#' RBAT(Us,Rs)
"MAData"
