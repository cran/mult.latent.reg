#' International Adult Literacy Survey (IALS) for 13 countries
#'
#' The data is obtained from the International Adult Literacy Survey (IALS),
#' collected in 13 countries on Prose, Document, and Quantitative scales between 1994 and 1995.
#' The data are reported as the percentage of individuals who could not reach a basic level of literacy in each country.
#'
#' @docType data
#'
#' @usage data(IALS_data)
#'
#' @format An object of class \code{"data.frame"}
#' \describe{
#'  \item{Prose}{On prose scale, the percentage of individuals who could not reach a basic level of literacy in each country.}
#'  \item{Document}{On document scale, the percentage of individuals who could not reach a basic level of literacy in each country.}
#'  \item{Quantitative}{On quantitative scale, the percentage of individuals who could not reach a basic level of literacy in each country.}
#'  \item{Country}{Specify the country}
#'  \item{Gender}{Specify the gender}
#' }
#' @references Sofroniou, N., Hoad, D., & Einbeck, J. (2008).  League tables for literacy survey data based on random effect models.
#'             In: Proceedings of the 23rd International Workshop on Statistical Modelling, Utrecht; pp. 402-405.
#' @keywords datasets
#' @examples
#'
#' data(IALS_data)
#' head(IALS_data)
"IALS_data"
