#' @title Selecting the best results for multivariate two level model
#' @name mult.reg_2level
#' @description This wrapper function runs multiple times the function \link{mult.em_2level} for fitting Zhang et al.'s (2023) multivariate response models with two-level random effect,
#'               and select the best results with the smallest AIC value.
#' @param data A data set object; we denote the dimension of a data set to be \eqn{m}.
#' @param v Covariate(s).
#' @param K Number of mixture components, the default is \code{K = 2}.
#' @param steps Number of iterations within each \code{num_runs}, the default is \code{steps = 20}.
#' @param num_runs Number of function iteration runs, the default is \code{num_runs = 20}.
#' @param start Containing parameters involved in the proposed model (\code{p}, \code{alpha}, \code{z}, \code{beta}, \code{sigma}, \code{gamma}) in a list,
#'              the starting values can be obtained through the use of \link{start_em}. More details can be found in \link{start_em}.
#' @param option Four options for selecting the starting values for the parameters in the model. The default is \code{option = 1}.
#'               More details can be found in \link{start_em}.
#' @param var_fun There are two types of variance specifications; \code{var_fun = 1}, the same diagonal variance specification to all \code{K} components of the mixture;
#'                \code{var_fun = 2}, different diagonal variance matrices for different components;
#'                The default is \code{var_fun = 2}.
#' @references Zhang, Y., Einbeck, J. and Drikvandi, R. (2023). A multilevel multivariate response model for data with latent structures.
#'             In: Proceedings of the 37th International Workshop on Statistical Modelling, pages 343-348.
#'              Link on RG: \url{https://www.researchgate.net/publication/375641972_A_multilevel_multivariate_response_model_for_data_with_latent_structures}
#' @return The best estimated result (with the smallest AIC value) in the model \eqn{x_{ij} = \alpha + \beta z_k + \Gamma v_{ij} + \varepsilon_{ij} } obtained through the EM algorithm (Zhang et al., 2023),
#'        where the upper-level unit is indexed by \eqn{i}, and the lower-level unit is indexed by \eqn{j}.
#'  \item{p}{The estimates for the parameter \eqn{\pi_k}, which is a vector of length \eqn{K}.}
#'  \item{alpha}{The estimates for the parameter \eqn{\alpha}, which is a vector of length \eqn{m}.}
#'  \item{z}{The estimates for the parameter \eqn{z_k}, which is a vector of length \eqn{K}.}
#'  \item{beta}{The estimates for the parameter \eqn{\beta}, which is a vector of length \eqn{m}.}
#'  \item{gamma}{The estimates for the parameter \eqn{\Gamma}, which is a matrix.}
#'  \item{sigma}{The estimates for the parameter \eqn{\Sigma_k}.
#'                      When \code{var_fun = 1}, \eqn{\Sigma_k} is a diagonal matrix and \eqn{\Sigma_k = \Sigma}, and we obtain a vector of the diagonal elements;
#'                      When \code{var_fun = 2}, \eqn{\Sigma_k} is a diagonal matrix, and we obtain \code{K} vectors of the diagonal elements.
#'                      }
#'  \item{W}{The posterior probability matrix.}
#'  \item{loglikelihood}{The approximated log-likelihood of the fitted model.}
#'  \item{disparity}{The disparity (\code{-2logL}) of the fitted model.}
#'  \item{number_parameters}{The number of parameters estimated in the EM algorithm.}
#'  \item{AIC}{The AIC value (\code{-2logL + 2number_parameters}).}
#'  \item{aic_data}{All AIC values in each run.}
#'  \item{Starting_values}{Lists of starting values for parameters used in each \code{num_runs}.
#'                      It allows reproduction of the best result (obtained from \link{mult.reg_2level}) in a single run
#'                      using \link{mult.em_2level} by setting \code{start} equal to the list of starting values
#'                      that were used to obtain the best result in \link{mult.reg_2level}.}
#' @seealso \code{\link{mult.em_2level}}.
#' @examples
#' \donttest{
#' ##run the mult.em_2level() multiple times and select the best results with the smallest AIC value
#' set.seed(7)
#' results <- mult.reg_2level(trading_data, K=4, steps = 10, num_runs = 5,
#'                            var_fun = 2, option = 1)
#' ## Reproduce the best result: the best result is the 2nd run in the above example.
#' rep_best_result <- mult.em_2level(trading_data, K=4, steps = 10,
#' var_fun = 2, option = 1,
#' start = results$Starting_values[[2]])
#' }
#' @import "mvtnorm"
#' @import "stats"
#' @import "utils"
#' @import "matrixStats"
#' @import "lme4"
#' @export
library(matrixStats)
library(mvtnorm)
mult.reg_2level <- function(data, v, start, K = 2, steps = 20, num_runs = 10, var_fun = 2, option = 1) {

  results <- list()
  res_aic <- numeric(num_runs)
  start_values <- list()

  for (i in 1: num_runs) {
    results[[i]] <-mult.em_2level(data=data, v=v, K=K, start=start, steps=steps, var_fun=var_fun, option=option)
    aic <- results[[i]]$AIC
    res_aic[i] <- aic
    start_values[[i]] <- results[[i]]$starting_values
  }

  min_aic_index <- which.min(res_aic)
  best_result <- results[[min_aic_index]]
  aic_data <- data.frame(Index = seq(1,num_runs), AIC = res_aic)

  return(list(best_result = best_result, aic_data = aic_data, Starting_values = start_values))
}


