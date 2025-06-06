#' @name mult.latent.reg
#' @aliases package
#' @title Regression and Clustering in Multivariate Response Scenarios
#' @description
#'    This package implements methodology for the estimation of multivariate response models with random effects on one or two levels;
#'    whereby the (one-dimensional) random effect represents a latent variable approximating the multivariate space of outcomes,
#'    after possible adjustment for covariates. The estimation methodology makes use of a nonparametric maximum likelihood-type approach,
#'    where the random effect distribution is approximated by a discrete mixture, hence allowing the use of the EM algorithm for the estimation of all model parameters.
#'    The method is particularly useful for multivariate,
#'    highly correlated outcome variables with unobserved heterogeneities. Applications include regression with multivariate responses,
#'    as well as multivariate clustering or ranking problems.
#'    The details of the models can be found in Zhang and Einbeck (2024) and Zhang et al. (2023).
#'    The main functions are \code{\link{mult.em_1level}} and \code{\link{mult.em_2level}} for the fitting of the raw models, as well as envelope functions
#'    \code{\link{mult.reg_1level}} and \code{\link{mult.reg_2level}} which facilitate iterative runs of the algorithm with a view to
#'    finding optimal starting points, with help by function \code{\link{start_em}}.
#'
#' @details
#'Package: mult.latent.reg
#'
#'Type:    Package
#'
#'License: GPL-3
#'
#' @author Yingjuan Zhang <yingjuan.zhang7@gmail.com>
#' @author Jochen Einbeck <jochen.einbeck@durham.ac.uk>
#'
#' @references
#' Zhang, Y., Einbeck, J., and Drikvandi, R. (2023). A multilevel multivariate response model for data with latent structures.
#' In: Proceedings of the 37th International Workshop on Statistical Modelling, Dortmund; pages 343-348.
#' Link on RG: \url{https://www.researchgate.net/publication/375641972_A_multilevel_multivariate_response_model_for_data_with_latent_structures}.
#'
#' @references Zhang, Y. and Einbeck, J. (2024). A Versatile Model for Clustered and Highly Correlated Multivariate Data. J Stat Theory Pract 18(5).\doi{10.1007/s42519-023-00357-0}
#'
NULL
