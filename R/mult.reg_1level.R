#' @title Selecting the best results for multivariate one level model
#' @name mult.reg_1level
#' @description Run multiple times the function \link{mult.em_1level} for fitting Zhang and Einbeck's (2024) multivariate response models with one-level random effect,
#'              and select the best results with the smallest AIC value.
#' @param data A data set object; we denote the dimension of a data set to be \eqn{m}.
#' @param v Covariate(s).
#' @param K Number of mixture components, the default is \code{K = 2}.
#' @param steps Number of iterations within each \code{num_runs}, the default is \code{steps = 20}.
#' @param num_runs Number of function iteration runs, the default is \code{num_runs = 10}.
#' @param start Containing parameters involved in the proposed model (\code{p}, \code{alpha}, \code{z}, \code{beta}, \code{sigma}, \code{gamma}) in a list,
#'              the starting values can be obtained through the use of \link{start_em}. More details can be found in \link{start_em}.
#' @param option Four options for selecting the starting values for the parameters in the model. The default is \code{option = 1}.
#'                More details can be found in \link{start_em}.
#' @param var_fun There are four types of variance specifications;
#'                \code{var_fun = 1}, the same diagonal variance specification to all \code{K} components of the mixture;
#'                \code{var_fun = 2}, different diagonal variance matrices for different components.
#'                \code{var_fun = 3}, the same full (unrestricted) variance for all components.
#'                \code{var_fun = 4}, different full (unrestricted) variance matrices for different components.
#'                The default is \code{var_fun = 2}.
#' @return The best estimated result (with the smallest AIC value) in the model (Zhang and Einbeck, 2024) \eqn{x_{i} = \alpha + \beta z_k + \Gamma v_i + \varepsilon_i} obtained through the EM algorithm.
#'  \item{p}{The estimates for the parameter \eqn{\pi_k}, which is a vector of length \eqn{K}.}
#'  \item{alpha}{The estimates for the parameter \eqn{\alpha}, which is a vector of length \eqn{m}.}
#'  \item{z}{The estimates for the parameter \eqn{z_k}, which is a vector of length \eqn{K}.}
#'  \item{beta}{The estimates for the parameter \eqn{\beta}, which is a vector of length \eqn{m}.}
#'  \item{gamma}{The estimates for the parameter \eqn{\Gamma}, which is a matrix.}
#'  \item{sigma}{The estimates for the parameter \eqn{\Sigma_k}.
#'                      When \code{var_fun = 1}, \eqn{\Sigma_k} is a diagonal matrix and \eqn{\Sigma_k = \Sigma}, and we obtain a vector of the diagonal elements;
#'                      When \code{var_fun = 2}, \eqn{\Sigma_k} is a diagonal matrix, and we obtain \code{K} vectors of the diagonal elements;
#'                      When \code{var_fun = 3}, \eqn{\Sigma_k} is a full variance-covariance matrix, \eqn{\Sigma_k = \Sigma}, and we obtain a matrix \eqn{\Sigma};
#'                      When \code{var_fun = 4}, \eqn{\Sigma_k} is a full variance-covariance matrix, and we obtain \code{K} different matrices \eqn{\Sigma_k}.}
#'  \item{W}{The posterior probability matrix.}
#'  \item{loglikelihood}{The approximated log-likelihood of the fitted model.}
#'  \item{disparity}{The disparity (\code{-2logL}) of the fitted model.}
#'  \item{number_parameters}{The number of parameters estimated in the EM algorithm.}
#'  \item{AIC}{The AIC value (\code{-2logL + 2number_parameters}).}
#'  \item{BIC}{The BIC value (\code{-2logL + number_parameters*log(n)}), where n is the number of observations.}
#'  \item{aic_data}{All AIC values in each run.}
#'
#' @seealso \code{\link{mult.em_1level}}.
#' @references Zhang, Y. and Einbeck J.  (2024). A Versatile Model for Clustered and Highly Correlated Multivariate Data. J Stat Theory Pract 18(5).\doi{10.1007/s42519-023-00357-0}
#' @examples
#' \donttest{
#' ##run the mult.em_1level() multiple times and select the best results with the smallest AIC value
#' results <- mult.reg_1level(fetal_covid_data[,c(1:5)],v=fetal_covid_data$status_bi,
#' K=3, num_runs = 5,
#' steps = 20, var_fun = 2, option = 1)
#' }
#' @import "mvtnorm"
#' @import "stats"
#' @import "utils"
#' @import "matrixStats"
#' @export
library(matrixStats)
library(mvtnorm)
library(stats)
library(utils)
mult.reg_1level <- function(data, v, start, K = 2, steps = 10, num_runs = 10, var_fun = 2, option = 1) {

  results <- list()
  res_aic <- numeric(num_runs)

  for (i in 1: num_runs) {
    results[[i]] <- mult.em_1level(data=data, v=v, K=K, start=start, steps=steps, var_fun = var_fun, option = option)
    aic <- results[[i]]$AIC
    res_aic[i] <- aic
  }

  min_aic_index <- which.min(res_aic)
  best_result <- results[[min_aic_index]]
  aic_data <- data.frame(Index = seq(1,num_runs), AIC = res_aic)

  return(list(best_result = best_result, aic_data = aic_data))
}

