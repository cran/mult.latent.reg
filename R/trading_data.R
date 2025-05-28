#' A set of import and export data in 44 countries.
#'
#' The variables are given as the percentage of imports and exports in relation to the overall GDP.
#' The data set comprises data from 44 countries (for our analysis),
#' we specifically selected the time period between 2018 and 2022.
#'
#' @docType data
#'
#' @usage data(trading_data)
#'
#' @format An object of class \code{"data.frame"}
#' \describe{
#'  \item{import}{The country-wise percentages of imports in relation to the overall GDP in each country.}
#'  \item{export}{The country-wise percentages of exports in relation to the overall GDP in each country.}
#'  \item{country}{The name of the countries.}
#' }
#' @source Trade in Goods and Services. \url{https://www.oecd.org/en/data/indicators/trade-in-goods-and-services.html}. Accessed on 2023-05-29.
#' @keywords datasets
#' @examples
#'
#' data(trading_data)
#' head(trading_data)
"trading_data"
