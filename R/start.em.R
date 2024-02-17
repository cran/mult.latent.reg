#' @title Starting values for parameters
#' @name start_em
#' @description The starting values for parameters used for the EM algorithm in the functions: \link{mult.em_1level}, \link{mult.em_2level}, \link{mult.reg_1level} and \link{mult.reg_2level}.
#' @param data A data set object; we denote the dimension of a data set to be \eqn{m}.
#' @param v Covariate(s); we denote the dimension of it to be \eqn{r}.
#' @param K Number of mixture components, the default is \code{K = 2}.
#' @param steps Number of iterations. This will only be used when using \code{option = 2} for both the 1-level model and the 2-level model.
#'              It should also be used when using \code{option = 3} and \code{option = 4} for the 1-level model, provided \code{var_fun} is set to either 3 or 4;
#'              the default is \code{steps = 20}.
#' @param option Four options for selecting the starting values for the parameters. The default is \code{option = 1}.
#'               When \code{option = 1}: \eqn{\pi_k} = \eqn{\frac{1}{K}}, \eqn{z_k} ~ rnorm(\eqn{K}, mean = 0, sd=1), \eqn{\alpha} = column means, \eqn{\beta} = a random row minus alpha,
#'               \eqn{\Gamma} = coefficient estimates from separate linear models, \eqn{\Sigma} is diagonal matrix where the diagonals take the value of column standard deviations over \eqn{K};
#'               when \code{option = 2}: use a short run (\code{steps = 5}) of the EM function which uses \code{option = 1} with \code{var_fun = 1} and use the estimates as the starting values for all the parameters;
#'               when \code{option = 3}: the starting value of \eqn{\beta} is the first principal component, and the starting values for the rest of the parameters are the same as described when \code{option = 1};
#'               when \code{option = 4}: first, take the scores of the first principal component of the data and perform \eqn{K}-means, \eqn{\pi_k} is the proportion of the clustering assignments, and \eqn{z_k} take the values of the \eqn{K}-means centers,
#'               and the starting values for the rest of the parameters are the same as described when \code{option = 1}.
#' @param var_fun The four variance specifications. When \code{var_fun = 1}, the same diagonal variance specification to all \eqn{K} components of the mixture;
#'                \code{var_fun = 2}, different diagonal variance matrices for different components.
#'                \code{var_fun = 3}, the same full (unrestricted) variance for all components.
#'                \code{var_fun = 4}, different full (unrestricted) variance matrices for different components.
#'                If unspecified, \code{var_fun = 2}. Note that for application propose, in two-level models, \code{var_fun} can only take values of 1 or 2.
#' @param p optional; specifies starting values for \eqn{\pi_k}, it is input as a \eqn{K}-dimensional vector.
#'
#' @param z optional; specifies starting values for \eqn{z_k}, it is input as a \eqn{K}-dimensional vector.
#'
#' @param beta optional; specifies starting values for \eqn{\beta}, it is input as an \eqn{m}-dimensional vector.
#'
#' @param alpha optional; specifies starting values for \eqn{\alpha}, it is input as an \eqn{m}-dimensional vector.
#'
#' @param sigma optional; specifies starting values for \eqn{\Sigma_k} (\eqn{\Sigma}, when \code{var_fun = 1} or \code{var_fun = 3}), when \code{var_fun = 1}, it is input as an \eqn{m}-dimensional vector,
#'          when \code{var_fun = 2}, it is input as a list (of length \eqn{K}) of \eqn{m}-dimensional vectors, when \code{var_fun = 3}, it is input as an \eqn{m \times m} matrix,
#'          when \code{var_fun = 4}, it is input as a list (of length \eqn{K}) of \eqn{m \times m} matrices.
#' @param gamma optional; the coefficients for the covariates; specifies starting values for \eqn{\Gamma}, it is input as an \eqn{m \times r} matrix.
#'
#' @return The starting values (in a list) for parameters in the models \eqn{x_{i} = \alpha + \beta z_k + \Gamma v_i + \varepsilon_i} (Zhang and Einbeck, 2024) and
#'         \eqn{x_{ij} = \alpha + \beta z_k + \Gamma v_{ij} + \varepsilon_{ij}} (Zhang et al., 2023) used in the four fucntions: \link{mult.em_1level}, \link{mult.em_2level}, \link{mult.reg_1level} and \link{mult.reg_2level}.
#'  \item{p}{The starting value for the parameter \eqn{\pi_k}, which is a vector of length \eqn{K}.}
#'  \item{alpha}{The starting value for the parameter \eqn{\alpha}, which is a vector of length \eqn{m}.}
#'  \item{z}{The starting value for the parameter \eqn{z_k}, which is a vector of length \eqn{K}.}
#'  \item{beta}{The starting value for the parameter \eqn{\beta}, which is a vector of length \eqn{m}.}
#'  \item{gamma}{The starting value for the parameter \eqn{\Gamma}, which is a matrix.}
#'  \item{sigma}{The starting value for the parameter \eqn{\Sigma_k}.
#'                      When \code{var_fun = 1}, \eqn{\Sigma_k} is a diagonal matrix and \eqn{\Sigma_k = \Sigma}, and we obtain a vector of the diagonal elements;
#'                      When \code{var_fun = 2}, \eqn{\Sigma_k} is a diagonal matrix, and we obtain \code{K} vectors of the diagonal elements;
#'                      When \code{var_fun = 3}, \eqn{\Sigma_k} is a full variance-covariance matrix, \eqn{\Sigma_k = \Sigma}, and we obtain a matrix \eqn{\Sigma};
#'                      When \code{var_fun = 4}, \eqn{\Sigma_k} is a full variance-covariance matrix, and we obtain \code{K} different matrices \eqn{\Sigma_k}.}
#' @references
#' Zhang, Y., Einbeck, J. and Drikvandi, R. (2023). A multilevel multivariate response model for data with latent structures.
#' In: Proceedings of the 37th International Workshop on Statistical Modelling, pages 343-348.
#' Link on RG: \url{https://www.researchgate.net/publication/375641972_A_multilevel_multivariate_response_model_for_data_with_latent_structures}.
#' @references
#' Zhang, Y. and Einbeck, J. (2024). A Versatile Model for Clustered and Highly Correlated Multivariate Data. J Stat Theory Pract 18(5).\doi{10.1007/s42519-023-00357-0}
#' @examples
#' ##example for the faithful data.
#' data(faithful)
#' start <- start_em(faithful, option = 1)
#' @import "mvtnorm"
#' @import "stats"
#' @import "utils"
#' @import "matrixStats"
#' @export

library(matrixStats)
library(mvtnorm)
start_em <- function(data, v, p, alpha, beta, z, gamma, sigma, steps = 20, K = 2, var_fun = 2, option = 1 ){

  data <- data.frame(data)
  data_numeric <- data[sapply(data, is.numeric)]
  data_factor <- data[sapply(data, is.factor)]
  #data_factor <- as.factor(data[,names(Filter(is.factor,data))])
  m <- length(data_numeric[1,])
  n <- length(data_numeric[,1])

  if (missing(v)){
    if (option == 1) {

      if (missing(p)){
        p <- rep(1/K, K)

      }
      if (missing(z)){
        z <- rnorm(K, mean = 0, sd=1)
      }
      if (missing(alpha)){
        alpha <- colMeans(data_numeric)
      }
      if (missing(beta)){
        beta <- as.numeric(data_numeric[sample(nrow(data_numeric), 1), ] - alpha)
      }
      if (missing(sigma)){
        if(var_fun == 1){
          s_n <- c()
          for (j in 1:m) {
            s_n[j] <- sqrt(1/(n-1)*sum((data_numeric[,j] - mean(data_numeric[,j]))^2))
          }
          sigma <- s_n/K
        }
        if(var_fun == 2){
          s_n <- c()
          sigma <- list()
          for (k in 1:K) {
            for (j in 1:m) {
              s_n[j] <- sqrt(1/(n-1)*sum((data_numeric[,j] - mean(data_numeric[,j]))^2))/K
            }
            sigma[[k]] <- s_n}
        }
        if(var_fun == 3){
          s_n <- c()
          sigma <- matrix(0,m,m)
          for (k in 1:K) {
            for (j in 1:m) {
              s_n[j] <- (1/(n-1)*sum((data[,j] - mean(data[,j]))^2))/K
            }
            mat <- matrix(0,m,m)
            diag(mat) <- s_n
            sigma <- mat
          }
        }
        if(var_fun == 4){
          s_n <- c()
          sigma <- vector("list",length = K)
          for (k in 1:K) {
            for (j in 1:m) {
              s_n[j] <- (1/(n-1)*sum((data[,j] - mean(data[,j]))^2))/K
            }
            mat <- matrix(0,m,m)
            diag(mat) <- s_n
            sigma[[k]] <- mat
          }
        }
      }
    }

    if (option == 2) {
      if (any(sapply(data, is.factor))) {
        if(var_fun == 1){
          res_1 <- em_2level_fun.start(data = data, K = K, steps = 5, var_fun = 1)
          if(missing(p)){p <- res_1$p}
          if(missing(alpha)){alpha <- res_1$alpha}
          if(missing(beta)){beta <- res_1$beta_hat}
          if(missing(z)){z  <- res_1$z_hat}
          if(missing(sigma)){sigma <- res_1$sigma}
        }

        if(var_fun == 2){
          res_2 <- em_2level_fun.start(data = data, K = K, steps = 5, var_fun = 2)
          if(missing(p)){p <- res_2$p}
          if(missing(alpha)){alpha <- res_2$alpha}
          if(missing(beta)){beta <- res_2$beta_hat}
          if(missing(z)){z  <- res_2$z_hat}
          if(missing(sigma)){sigma <- res_2$sigma}
        }

      } else {

        if(var_fun == 1){
          res_1 <- em_fun.start(data, K, steps = 5, var_fun = 1)
          if(missing(p)){p <- res_1$p}
          if(missing(alpha)){alpha <- res_1$alpha}
          if(missing(beta)){beta <- res_1$beta_hat}
          if(missing(z)){z  <- res_1$z_hat}
          if(missing(sigma)){sigma <- res_1$sigma}
        }

        if(var_fun == 2){
          res_2 <- em_fun.start(data, K, steps = 5, var_fun = 1)
          if(missing(p)){p <- res_2$p}
          if(missing(alpha)){alpha <- res_2$alpha}
          if(missing(beta)){beta <- res_2$beta_hat}
          if(missing(z)){z  <- res_2$z_hat}
          if(missing(sigma))
            {
            s_n <- res_2$sigma
            sigma <- list()
             for (k in 1:K) {
              sigma[[k]] <- s_n
              }
          }
        }
        if(var_fun == 3){
          res_3 <- em_fun.start(data, K, steps = 5, var_fun = 1)
          if(missing(p)){p <- res_3$p}
          if(missing(alpha)){alpha <- res_3$alpha}
          if(missing(beta)){beta <- res_3$beta_hat}
          if(missing(z)){z  <- res_3$z_hat}
          if(missing(sigma))
          {
            s_n <- res_3$sigma
            sigma <- matrix(0,m,m)
            for (k in 1:K) {
              mat <- matrix(0,m,m)
              diag(mat) <- s_n
              sigma <- mat
            }
            }
        }

        if(var_fun == 4){
          res_4 <- em_fun.start(data, K, steps = 5, var_fun = 1)
          if(missing(p)){p <- res_4$p}
          if(missing(alpha)){alpha <- res_4$alpha}
          if(missing(beta)){beta <- res_4$beta_hat}
          if(missing(z)){z  <- res_4$z_hat}
          if(missing(sigma))
          {
            s_n <- res_4$sigma
            sigma <- vector("list",length = K)
            for (k in 1:K) {
              mat <- matrix(0,m,m)
              diag(mat) <- s_n
              sigma[[k]] <- mat
            }
            }
        }

      }
    }


    if (option == 3) {
      if(missing(beta)){
        scaled_data <- scale(data_numeric)
        pca_result <- prcomp(scaled_data)
        beta <- pca_result$rotation[,1]
      }
      if (missing(p)){
        p <- rep(1/K, K)
      }
      if (missing(z)){
        z <- rnorm(K, mean = 0, sd=1)
      }
      if (missing(alpha)){
        alpha <- colMeans(data_numeric)
      }
      if (missing(sigma)){
        if(var_fun == 1){
          s_n <- c()
          for (j in 1:m) {
            s_n[j] <- sqrt(1/(n-1)*sum((data_numeric[,j] - mean(data_numeric[,j]))^2))
          }
          sigma <- s_n/K
        }
        if(var_fun == 2){
          s_n <- c()
          sigma <- list()
          for (k in 1:K) {
            for (j in 1:m) {
              s_n[j] <- sqrt(1/(n-1)*sum((data_numeric[,j] - mean(data_numeric[,j]))^2))/K
            }
            sigma[[k]] <- s_n}
        }
        if(var_fun == 3){
          s_n <- c()
          sigma <- matrix(0,m,m)
          for (k in 1:K) {
            for (j in 1:m) {
              s_n[j] <- (1/(n-1)*sum((data[,j] - mean(data[,j]))^2))/K
            }
            mat <- matrix(0,m,m)
            diag(mat) <- s_n
            sigma <- mat
          }
        }
        if(var_fun == 4){
          s_n <- c()
          sigma <- vector("list",length = K)
          for (k in 1:K) {
            for (j in 1:m) {
              s_n[j] <- (1/(n-1)*sum((data[,j] - mean(data[,j]))^2))/K
            }
            mat <- matrix(0,m,m)
            diag(mat) <- s_n
            sigma[[k]] <- mat
          }
        }
      }
    }

    if(option == 4){
      pca_result <- prcomp(data_numeric)
      score <- pca_result$x[,1]
      kmeans_result <- kmeans(score, centers = K)
      cluster_assignments <- kmeans_result$cluster

      if(missing(p)){
        p <- as.vector(prop.table(table(cluster_assignments)))
      }
      if(missing(z)){
        z <- as.vector(kmeans_result$centers)
      }
      if(missing(alpha)){
        alpha <- colMeans(data_numeric)
      }
      if (missing(beta)){
        beta <- as.numeric(data_numeric[sample(nrow(data_numeric), 1), ] - alpha)
      }
      if (missing(sigma)){
        if(var_fun == 1){
          s_n <- c()
          for (j in 1:m) {
            s_n[j] <- sqrt(1/(n-1)*sum((data_numeric[,j] - mean(data_numeric[,j]))^2))
          }
          sigma <- s_n/K
        }
        if(var_fun == 2){
          s_n <- c()
          sigma <- list()
          for (k in 1:K) {
            for (j in 1:m) {
              s_n[j] <- sqrt(1/(n-1)*sum((data_numeric[,j] - mean(data_numeric[,j]))^2))/K
            }
            sigma[[k]] <- s_n}
        }
        if(var_fun == 3){
          s_n <- c()
          sigma <- matrix(0,m,m)
          for (k in 1:K) {
            for (j in 1:m) {
              s_n[j] <- (1/(n-1)*sum((data[,j] - mean(data[,j]))^2))/K
            }
            mat <- matrix(0,m,m)
            diag(mat) <- s_n
            sigma <- mat
          }
        }
        if(var_fun == 4){
          s_n <- c()
          sigma <- vector("list",length = K)
          for (k in 1:K) {
            for (j in 1:m) {
              s_n[j] <- (1/(n-1)*sum((data[,j] - mean(data[,j]))^2))/K
            }
            mat <- matrix(0,m,m)
            diag(mat) <- s_n
            sigma[[k]] <- mat
          }
        }
      }
    }

    start <- list(p = p, alpha = alpha, beta = beta, z = z, sigma = sigma)

  } else{

    v <- as.data.frame(v)

    if (option == 1) {
      if (missing(p)){
        p <- rep(1/K, K)
      }
      if (missing(z)){
        z <- rnorm(K, mean = 0, sd=1)
      }
      if (missing(alpha)){
        alpha <- colMeans(data_numeric)
      }
      if (missing(beta)){
        beta <- as.numeric(data_numeric[sample(nrow(data_numeric), 1), ] - alpha)
      }
      if(missing(gamma)){
        m <- as.numeric(length(data_numeric[1,]))
        q <- as.numeric(length(v[1,]))
        dep_vars <- colnames(data_numeric)
        ind_vars <- colnames(v)
        var_comb <- expand.grid(dep_vars, ind_vars )
        formula_vec <- sprintf("%s ~ %s", var_comb$Var1, var_comb$Var2)
        lm_res <- lapply( formula_vec, function(f) {
          fit1 <- lm( f, data = data.frame(data_numeric,v))
          return(fit1$coefficients[2])
        })
        gamma_new <- as.numeric(unlist(lm_res))
        gamma <- matrix(gamma_new,m,q)
      }
      if (missing(sigma)){
        if(var_fun == 1){
          s_n <- c()
          for (j in 1:m) {
            s_n[j] <- sqrt(1/(n-1)*sum((data_numeric[,j] - mean(data_numeric[,j]))^2))
          }
          sigma <- s_n/K
        }
        if(var_fun == 2){
          s_n <- c()
          sigma <- list()
          for (k in 1:K) {
            for (j in 1:m) {
              s_n[j] <- sqrt(1/(n-1)*sum((data_numeric[,j] - mean(data_numeric[,j]))^2))/K
            }
            sigma[[k]] <- s_n}
        }
        if(var_fun == 3){
          s_n <- c()
          sigma <- matrix(0,m,m)
          for (k in 1:K) {
            for (j in 1:m) {
              s_n[j] <- (1/(n-1)*sum((data[,j] - mean(data[,j]))^2))/K
            }
            mat <- matrix(0,m,m)
            diag(mat) <- s_n
            sigma <- mat
          }
        }
        if(var_fun == 4){
          s_n <- c()
          sigma <- vector("list",length = K)
          for (k in 1:K) {
            for (j in 1:m) {
              s_n[j] <- (1/(n-1)*sum((data[,j] - mean(data[,j]))^2))/K
            }
            mat <- matrix(0,m,m)
            diag(mat) <- s_n
            sigma[[k]] <- mat
          }
        }
      }
    }
    if (option == 2) {
      if (any(sapply(data, is.factor))) {
      if(var_fun == 1){
        res_1 <- em_2level_covs.start(data, v, K, steps = 5, var_fun = 1)
        if(missing(p)){p <- res_1$p}
        if(missing(alpha)){alpha <- res_1$alpha}
        if(missing(beta)){beta <- res_1$beta_hat}
        if(missing(z)){z  <- res_1$z_hat}
        if(missing(sigma)){sigma <- res_1$sigma}
        if(missing(gamma)){gamma <- res_1$gamma}
      }

      if(var_fun == 2){
        res_2 <- em_2level_covs.start(data, v, K, steps = 5, var_fun = 1)
        if(missing(p)){p <- res_2$p}
        if(missing(alpha)){alpha <- res_2$alpha}
        if(missing(beta)){beta <- res_2$beta_hat}
        if(missing(z)){z  <- res_2$z_hat}
        if(missing(sigma)){
          s_n <- res_2$sigma
          sigma <- list()
          for (k in 1:K) {
            sigma[[k]] <- s_n
          }
        }
        if(missing(gamma)){gamma <- res_2$gamma}
      }

      } else {
          if(var_fun == 1){
            res_1 <- em_covs.start(data, v, K, steps = 5, var_fun = 1)
            if(missing(p)){p <- res_1$p}
            if(missing(alpha)){alpha <- res_1$alpha}
            if(missing(beta)){beta <- res_1$beta_hat}
            if(missing(z)){z  <- res_1$z_hat}
            if(missing(sigma)){
              sigma <- res_1$sigma
              }
            if(missing(gamma)){gamma <- res_1$gamma}
          }
        if(var_fun == 2){
          res_2 <- em_covs.start(data, v, K, steps = 5, var_fun = 1)
          if(missing(p)){p <- res_2$p}
          if(missing(alpha)){alpha <- res_2$alpha}
          if(missing(beta)){beta <- res_2$beta_hat}
          if(missing(z)){z  <- res_2$z_hat}
          if(missing(sigma)){
            s_n <- res_2$sigma
            sigma <- list()
            for (k in 1:K) {
              sigma[[k]] <- s_n
            }
            }
          if(missing(gamma)){gamma <- res_2$gamma}
        }
        if(var_fun == 3){
          res_3 <- em_covs.start(data, v, K, steps = 5, var_fun = 1)
          res_3_p <- res_3$p[1]
          if (is.na(res_3_p)) {
            repeat {
              res_3 <- em_covs.start(data, v, K, steps = 5, var_fun = 1)
              if (!anyNA(res_3_p)) {
                break
              }
            }
          }
          if(missing(p)){p <- res_3$p}
          if(missing(alpha)){alpha <- res_3$alpha}
          if(missing(beta)){beta <- res_3$beta}
          if(missing(z)){z  <- res_3$z}
          if(missing(sigma)){
            s_n <- res_3$sigma
            sigma <- matrix(0,m,m)
            for (k in 1:K) {
              mat <- matrix(0,m,m)
              diag(mat) <- s_n
              sigma <- mat
            }
            }
          if(missing(gamma)){gamma <- res_3$gamma}
        }
        if(var_fun == 4){
          res_4 <- em_covs.start(data, v, K, steps = 5, var_fun = 1)
          res_4_p <- res_4$p[1]
          if (is.na(res_4_p)) {
            repeat {
              res_4 <- em_covs.start(data, v, K, steps = 5, var_fun = 1)
              if (!anyNA(res_4_p)) {
                break
              }
            }
          }
          if(missing(p)){p <- res_4$p}
          if(missing(alpha)){alpha <- res_4$alpha}
          if(missing(beta)){beta <- res_4$beta}
          if(missing(z)){z  <- res_4$z}
          if(missing(sigma)){
            s_n <- res_4$sigma
            sigma <- vector("list",length = K)
            for (k in 1:K) {
              mat <- matrix(0,m,m)
              diag(mat) <- s_n
              sigma[[k]] <- mat
            }
            }
          if(missing(gamma)){gamma <- res_4$gamma}
        }
        }
    }

    if (option == 3) {
      if(missing(beta)){
        scaled_data <- scale(data_numeric)
        pca_result <- prcomp(scaled_data)
        beta <- pca_result$rotation[,1]
      }
      if (missing(p)){
        p <- rep(1/K, K)
      }
      if (missing(z)){
        z <- rnorm(K, mean = 0, sd=1)
      }
      if (missing(alpha)){
        alpha <- colMeans(data_numeric)
      }
      if(missing(gamma)){
        m <- as.numeric(length(data_numeric[1,]))
        q <- as.numeric(length(v[1,]))
        dep_vars <- colnames(data_numeric)
        ind_vars <- colnames(v)
        var_comb <- expand.grid(dep_vars, ind_vars )
        formula_vec <- sprintf("%s ~ %s", var_comb$Var1, var_comb$Var2)
        lm_res <- lapply( formula_vec, function(f) {
          fit1 <- lm( f, data = data.frame(data_numeric,v))
          return(fit1$coefficients[2])
        })
        gamma_new <- as.numeric(unlist(lm_res))
        gamma <- matrix(gamma_new,m,q)
      }
      if (missing(sigma)){
        if(var_fun == 1){
          s_n <- c()
          for (j in 1:m) {
            s_n[j] <- sqrt(1/(n-1)*sum((data_numeric[,j] - mean(data_numeric[,j]))^2))
          }
          sigma <- s_n/K
        }
        if(var_fun == 2){
          s_n <- c()
          sigma <- list()
          for (k in 1:K) {
            for (j in 1:m) {
              s_n[j] <- sqrt(1/(n-1)*sum((data_numeric[,j] - mean(data_numeric[,j]))^2))/K
            }
            sigma[[k]] <- s_n}
        }
        if(var_fun == 3){
          s_n <- c()
          sigma <- matrix(0,m,m)
          for (k in 1:K) {
            for (j in 1:m) {
              s_n[j] <- (1/(n-1)*sum((data[,j] - mean(data[,j]))^2))/K
            }
            mat <- matrix(0,m,m)
            diag(mat) <- s_n
            sigma <- mat
          }
        }
        if(var_fun == 4){
          s_n <- c()
          sigma <- vector("list",length = K)
          for (k in 1:K) {
            for (j in 1:m) {
              s_n[j] <- (1/(n-1)*sum((data[,j] - mean(data[,j]))^2))/K
            }
            mat <- matrix(0,m,m)
            diag(mat) <- s_n
            sigma[[k]] <- mat
          }
        }
      }
    }

    if(option == 4){
      pca_result <- prcomp(data_numeric)
      score <- pca_result$x[,1]
      kmeans_result <- kmeans(score, centers = K)
      cluster_assignments <- kmeans_result$cluster

      if(missing(p)){
        p <- as.vector(prop.table(table(cluster_assignments)))
      }
      if(missing(z)){
        z <- as.vector(kmeans_result$centers)
      }
      if(missing(alpha)){
        alpha <- colMeans(data_numeric)
      }
      if (missing(beta)){
        beta <- as.numeric(data_numeric[sample(nrow(data_numeric), 1), ] - alpha)
      }
      if(missing(gamma)){
        m <- as.numeric(length(data_numeric[1,]))
        q <- as.numeric(length(v[1,]))
        dep_vars <- colnames(data_numeric)
        ind_vars <- colnames(v)
        var_comb <- expand.grid(dep_vars, ind_vars )
        formula_vec <- sprintf("%s ~ %s", var_comb$Var1, var_comb$Var2)
        lm_res <- lapply( formula_vec, function(f) {
          fit1 <- lm( f, data = data.frame(data_numeric,v))
          return(fit1$coefficients[2])
        })
        gamma_new <- as.numeric(unlist(lm_res))
        gamma <- matrix(gamma_new,m,q)
      }
      if (missing(sigma)){
        if(var_fun == 1){
          s_n <- c()
          for (j in 1:m) {
            s_n[j] <- sqrt(1/(n-1)*sum((data_numeric[,j] - mean(data_numeric[,j]))^2))
          }
          sigma <- s_n/K
        }
        if(var_fun == 2){
          s_n <- c()
          sigma <- list()
          for (k in 1:K) {
            for (j in 1:m) {
              s_n[j] <- sqrt(1/(n-1)*sum((data_numeric[,j] - mean(data_numeric[,j]))^2))/K
            }
            sigma[[k]] <- s_n}
        }
        if(var_fun == 3){
          s_n <- c()
          sigma <- matrix(0,m,m)
          for (k in 1:K) {
            for (j in 1:m) {
              s_n[j] <- (1/(n-1)*sum((data[,j] - mean(data[,j]))^2))/K
            }
            mat <- matrix(0,m,m)
            diag(mat) <- s_n
            sigma <- mat
          }
        }
        if(var_fun == 4){
          s_n <- c()
          sigma <- vector("list",length = K)
          for (k in 1:K) {
            for (j in 1:m) {
              s_n[j] <- (1/(n-1)*sum((data[,j] - mean(data[,j]))^2))/K
            }
            mat <- matrix(0,m,m)
            diag(mat) <- s_n
            sigma[[k]] <- mat
          }
        }
      }
    }


    start <- list(p = p, alpha = alpha, beta = beta, z = z, gamma = gamma, sigma = sigma)

  }

  return(start)
}




#######em_fun.start()
estep<- function(data, p, z, alpha, beta, sigma, var_fun){

  data <- data.frame(data)

  m <- length(data[1,])

  if(missing(var_fun)){
    var_fun <- 2
  }

  if(m != 1){

    n <- length(data[,1])
    K <- length(p)
    m <- length(data[1,])

    mu <- matrix(0,K,m)
    for (k in 1:K) {
      mu[k,] <- alpha + beta*z[k]
    }
    if(var_fun == 1) {
      sigma_matrix <- diag(sigma^2,m,m)
      W_1 <- matrix(0,n,K)
      for (k in 1:K) {
        W_1[,k] <- p[k]*dmvnorm(x=data, mean=mu[k,], sigma=sigma_matrix, log=FALSE)
      }
    }
    if(var_fun == 2){
      W_1 <- matrix(0,n,K)
      for (k in 1:K) {
        W_1[,k] <- p[k]*dmvnorm(x=data, mean=mu[k,], sigma=diag(unlist(sigma[k])^2,m,m), log=FALSE)
      }
    }
    if(var_fun == 3){
      W_1 <- matrix(0,n,K)
      for (k in 1:K) {
        W_1[,k] <- p[k]*dmvnorm(x=data, mean=mu[k,], sigma=matrix(sigma,m,m), log=FALSE)
      }
    }
    if(var_fun == 4){
      W_1 <- matrix(0,n,K)
      for (k in 1:K) {
        W_1[,k] <- p[k]*dmvnorm(x=data, mean=mu[k,], sigma= matrix(unlist(sigma[k]),m,m), log=FALSE)
      }
    }
    W <- W_1/apply(W_1,1,sum)
    return(W)
  }

  if(m == 1){
    length(data[,1])
    K <- length(p)
    mu <- matrix(0,K,m)
    for (k in 1:K) {
      mu[k,] <- alpha + beta*z[k]
    }
    W_1 <- matrix(0,n,K)
    for (k in 1:K) {
      W_1[,k] <- p[k]*dnorm(x=data, mean=mu[k,], sd=sqrt(sigma), log=FALSE)
    }
    W <- W_1/apply(W_1,1,sum)
    return(W)
  }

}


mstep<- function(data, W, alpha, beta, sigma, var_fun){
  data <- data.frame(data)

  m <- length(data[1,])

  if(missing(var_fun)){
    var_fun <- 2
  }

  if(m != 1){
    m <- length(data[1,])
  }
  n <- dim(W)[1]
  K <- dim(W)[2]


  if(m != 1){

    if(var_fun == 1){
      p <- apply(W,2,sum)/n
      counter <- 0
      while (counter <= 10) {
        z_1 <- matrix(0,K,m)
        for (k in 1:K) {
          z_1[k,] <- apply(W[,k]*data,2,sum)
        }
        z_2 <- matrix(0,K,m)
        for (j in 1:m) {
          z_2[,j]  <- beta[j]*z_1[,j]
        }
        z_3 <- apply(z_2,1,sum)
        z_4 <- matrix(0,K,m)
        for (j in 1:m) {
          z_4[,j] <- alpha[j]*beta[j]*apply(W,2,sum)
        }
        z_5 <- apply(z_4,1,sum)
        z_6 <- matrix(0,K,m)
        for (j in 1:m) {
          z_6[,j] <- ((beta[j])^2)*apply(W,2,sum)
        }
        z_7 <- apply(z_6,1,sum)

        z_new <- (z_3 - z_5)/z_7
        z_j <- z_new - sum(z_new*p)
        z <- z_j/sqrt(sum((z_j^2)*p))

        wz <- matrix(0,n,K)
        for (k in 1:K) {
          wz[,k] <- W[,k]*z[k]
        }

        wz_2 <- matrix(0,n,K)
        for (k in 1:K) {
          wz_2[,k] <- W[,k]*(z[k])^2
        }

        beta_1 <- apply(data*apply(wz,1,sum),2,sum)
        beta <- (beta_1 - (1/n)*(apply(data,2,sum)*sum(apply(wz,1,sum))))/(sum(apply(wz_2,1,sum)) - (1/n)*(sum(apply(wz,1,sum)))^2)

        alpha <- (1/n)*(apply(data,2,sum) - beta*sum(apply(wz,1,sum)))

        counter <- counter +1
      }

      diff<-matrix(0,n, K)
      sigma <- c()
      for (j in 1:m) {
        for (k in 1:K) {
          diff[,k] <- (data[,j] - alpha[j] - beta[j]*z[k])^2
        }
        sigma[j] <- sqrt(mean(apply(W*diff,1,sum)))
      }
    }

    if(var_fun == 2){
      p <- apply(W,2,sum)/n

      counter <- 0
      while (counter <= 10) {
        z_1 <- matrix(0,K,m)
        for (k in 1:K) {
          z_1[k,] <- apply(W[,k]*data,2,sum)
        }
        z_2 <- matrix(0,K,m)
        for (j in 1:m) {
          z_2[,j]  <- beta[j]*z_1[,j]
        }
        z_3 <- apply(z_2,1,sum)
        z_4 <- matrix(0,K,m)
        for (j in 1:m) {
          z_4[,j] <- alpha[j]*beta[j]*apply(W,2,sum)
        }
        z_5 <- apply(z_4,1,sum)
        z_6 <- matrix(0,K,m)
        for (j in 1:m) {
          z_6[,j] <- ((beta[j])^2)*apply(W,2,sum)
        }
        z_7 <- apply(z_6,1,sum)

        z_new <- (z_3 - z_5)/z_7
        z_j <- z_new - sum(z_new*p)
        z <- z_j/sqrt(sum((z_j^2)*p))

        wz <- matrix(0,n,K)
        for (k in 1:K) {
          wz[,k] <- W[,k]*z[k]
        }

        wz_2 <- matrix(0,n,K)
        for (k in 1:K) {
          wz_2[,k] <- W[,k]*(z[k])^2
        }

        beta_1 <- apply(data*apply(wz,1,sum),2,sum)
        beta <- (beta_1 - (1/n)*(apply(data,2,sum)*sum(apply(wz,1,sum))))/(sum(apply(wz_2,1,sum)) - (1/n)*(sum(apply(wz,1,sum)))^2)

        alpha <- (1/n)*(apply(data,2,sum) - beta*sum(apply(wz,1,sum)))

        counter <- counter +1
      }

      diff<- matrix(0,n,m)
      sigma <- c()
      for (k in 1:K) {
        for (j in 1:m) {
          diff[,j] <- (data[,j] - alpha[j] - beta[j]*z[k])^2
        }
        sigma[[k]] <- sqrt(apply(W[,k]*diff,2,sum)/sum(W[,k]))
      }
    }


    if(var_fun == 3){
      p <- apply(W,2,sum)/n

      counter <- 0
      while (counter <= 10) {
        z_1 <- matrix(0,K,m)
        for (k in 1:K) {
          z_1[k,] <- apply(W[,k]*data,2,sum)
        }
        z_2 <- matrix(0,K,m)
        for (j in 1:m) {
          z_2[,j]  <- beta[j]*z_1[,j]
        }
        z_3 <- apply(z_2,1,sum)
        z_4 <- matrix(0,K,m)
        for (j in 1:m) {
          z_4[,j] <- alpha[j]*beta[j]*apply(W,2,sum)
        }
        z_5 <- apply(z_4,1,sum)
        z_6 <- matrix(0,K,m)
        for (j in 1:m) {
          z_6[,j] <- ((beta[j])^2)*apply(W,2,sum)
        }
        z_7 <- apply(z_6,1,sum)

        z_new <- (z_3 - z_5)/z_7
        z_j <- z_new - sum(z_new*p)
        z <- z_j/sqrt(sum((z_j^2)*p))

        wz <- matrix(0,n,K)
        for (k in 1:K) {
          wz[,k] <- W[,k]*z[k]
        }

        wz_2 <- matrix(0,n,K)
        for (k in 1:K) {
          wz_2[,k] <- W[,k]*(z[k])^2
        }

        beta_1 <- apply(data*apply(wz,1,sum),2,sum)
        beta <- (beta_1 - (1/n)*(apply(data,2,sum)*sum(apply(wz,1,sum))))/(sum(apply(wz_2,1,sum)) - (1/n)*(sum(apply(wz,1,sum)))^2)

        alpha <- (1/n)*(apply(data,2,sum) - beta*sum(apply(wz,1,sum)))

        counter <- counter +1
      }

      current  <- vector("list", length = n)
      Sigma_1 <- vector("list",length = K)
      for (k in 1:K) {
        for (i in 1:n) {
          current[[i]] <- W[i,k]*t(as.matrix(data[i,] - alpha - beta*z[k]))%*%(as.matrix(data[i,] - alpha - beta*z[k]))
          sum_matrix <- Reduce("+", current)
        }
        Sigma_1[[k]] <- sum_matrix/n
        sigma <- Reduce("+", Sigma_1)
      }
    }

    if(var_fun == 4){

      p <- apply(W,2,sum)/n

      counter <- 0
      while (counter <= 10) {
        z_1 <- matrix(0,K,m)
        for (k in 1:K) {
          z_1[k,] <- apply(W[,k]*data,2,sum)
        }
        z_2 <- matrix(0,K,m)
        for (j in 1:m) {
          z_2[,j]  <- beta[j]*z_1[,j]
        }
        z_3 <- apply(z_2,1,sum)
        z_4 <- matrix(0,K,m)
        for (j in 1:m) {
          z_4[,j] <- alpha[j]*beta[j]*apply(W,2,sum)
        }
        z_5 <- apply(z_4,1,sum)
        z_6 <- matrix(0,K,m)
        for (j in 1:m) {
          z_6[,j] <- ((beta[j])^2)*apply(W,2,sum)
        }
        z_7 <- apply(z_6,1,sum)

        z_new <- (z_3 - z_5)/z_7
        z_j <- z_new - sum(z_new*p)
        z <- z_j/sqrt(sum((z_j^2)*p))

        wz <- matrix(0,n,K)
        for (k in 1:K) {
          wz[,k] <- W[,k]*z[k]
        }

        wz_2 <- matrix(0,n,K)
        for (k in 1:K) {
          wz_2[,k] <- W[,k]*(z[k])^2
        }

        beta_1 <- apply(data*apply(wz,1,sum),2,sum)
        beta <- (beta_1 - (1/n)*(apply(data,2,sum)*sum(apply(wz,1,sum))))/(sum(apply(wz_2,1,sum)) - (1/n)*(sum(apply(wz,1,sum)))^2)

        alpha <- (1/n)*(apply(data,2,sum) - beta*sum(apply(wz,1,sum)))

        counter <- counter +1
      }

      current  <- vector("list", length = n)
      sigma <- vector("list",length = K)
      for (k in 1:K) {
        for (i in 1:n) {
          current[[i]] <- W[i,k]*t(as.matrix(data[i,] - alpha - beta*z[k]))%*%(as.matrix(data[i,] - alpha - beta*z[k]))
          add <- function(x) Reduce("+", x)
          sum_matrix <- add(current)
        }
        sigma[[k]] <- sum_matrix/sum(W[,k])
      }


    }
  }

  if(m == 1){
    z_1 <- matrix(0,K,m)
    for (k in 1:K) {
      z_1[k,] <- sum(W[,k]*(data-alpha))
    }
    z_2 <- c()
    for (j in 1:m) {
      z_2 <- (beta)*apply(W,2,sum)
    }
    z_new <- z_1/z_2
    z_j <- z_new - sum(z_new*p)
    z <- z_j/sqrt(sum((z_j^2)*p))

    wz <- matrix(0,n,K)
    for (k in 1:K) {
      wz[,k] <- W[,k]*z[k]
    }

    wz_2 <- matrix(0,n,K)
    for (k in 1:K) {
      wz_2[,k] <- W[,k]*(z[k])^2
    }

    beta_1 <- sum(data*apply(wz,1,sum))
    beta <- (beta_1 - (1/n)*(sum(data)*sum(apply(wz,1,sum))))/(sum(apply(wz_2,1,sum)) - (1/n)*(sum(apply(wz,1,sum)))^2)

    alpha <- (1/n)*(sum(data) - beta*sum(apply(wz,1,sum)))

    diff<-matrix(0,n, K)
    sigma <- c()
    for (j in 1:m) {
      for (k in 1:K) {
        diff[,k] <- (data - alpha - beta*z[k])^2
      }
      sigma[j] <- sqrt(mean(apply(W*diff,1,sum)))
    }
  }

  return(list("p"=p, "z"=z, "alpha"=alpha, "beta"=beta, "sigma"=sigma))
}

em_fun.start <- function(data,K,steps,p,z,beta,alpha,sigma,var_fun){
  pb<-txtProgressBar(min=0, max=steps, style=3)
  data <- data.frame(data)
  m <- length(data[1,])
  n <- length(data[,1])

  if(missing(var_fun)){
    var_fun <- 2
  }

  if(m != 1){
    if (missing(p)){
      p <- rep(1/K, K)
    }
    if (missing(z)){
      z <- rnorm(K,mean = 0, sd=1)
    }
    if (missing(alpha)){
      alpha <- colMeans(data)
    }
    if (missing(beta)){
      beta <- as.numeric(data[sample(nrow(data), 1), ] - alpha)
    }
    if (missing(sigma)){
      if(var_fun == 1){
        s_n <- c()
        for (j in 1:m) {
          s_n[j] <- sqrt(1/(n-1)*sum((data[,j] - mean(data[,j]))^2))
        }
        sigma <- s_n/K
      }
      if(var_fun == 2){
        s_n <- c()
        sigma <- list()
        for (k in 1:K) {
          for (j in 1:m) {
            s_n[j] <- sqrt(1/(n-1)*sum((data[,j] - mean(data[,j]))^2))/K
          }
          sigma[[k]] <- s_n}
      }
      if(var_fun == 3){
        s_n <- c()
        sigma <- matrix(0,m,m)
        for (k in 1:K) {
          for (j in 1:m) {
            s_n[j] <- (1/(n-1)*sum((data[,j] - mean(data[,j]))^2))/K
          }
          mat <- matrix(0,m,m)
          diag(mat) <- s_n
          sigma <- mat
        }
      }
      if(var_fun == 4){
        s_n <- c()
        sigma <- vector("list",length = K)
        for (k in 1:K) {
          for (j in 1:m) {
            s_n[j] <- (1/(n-1)*sum((data[,j] - mean(data[,j]))^2))/K
          }
          mat <- matrix(0,m,m)
          diag(mat) <- s_n
          sigma[[k]] <- mat
        }
      }
      #if(var_fun == 3){
      # v3<-em.diag(data,K,steps=steps+20,var_fun = 3)
      #  sigma <- v3$sigma
      #}
      #if(var_fun == 4){
      # v3<-em.diag(data,K,steps=steps+20,var_fun = 4)
      #  sigma <- v3$sigma
      #}
    }
  }


  if(m != 1 && K == 1){
    if (missing(p)){
      p <- rep(1/K, K)
    }
    if (missing(z)){
      z <- rnorm(K,mean = 0, sd=1)
    }
    if (missing(alpha)){
      alpha <- colMeans(data)
    }
    if (missing(beta)){
      beta <- as.numeric(data[sample(nrow(data), 1), ] - alpha)
    }
    if (missing(sigma)){
      if(var_fun == 1){
        s_n <- c()
        for (j in 1:m) {
          s_n[j] <- sqrt(1/(n-1)*sum((data[,j] - mean(data[,j]))^2))
        }
        sigma <- s_n/K
      }
      if(var_fun == 2){
        s_n <- c()
        sigma <- list()
        for (k in 1:K) {
          for (j in 1:m) {
            s_n[j] <- sqrt(1/(n-1)*sum((data[,j] - mean(data[,j]))^2))/K
          }
          sigma[[k]] <- s_n}
      }
      if(var_fun == 3){
        s_n <- c()
        sigma <- matrix(0,m,m)
        for (k in 1:K) {
          for (j in 1:m) {
            s_n[j] <- (1/(n-1)*sum((data[,j] - mean(data[,j]))^2))/K
          }
          mat <- matrix(0,m,m)
          diag(mat) <- s_n
          sigma <- mat
        }
      }
      if(var_fun == 4){
        s_n <- c()
        sigma <- vector("list",length = K)
        for (k in 1:K) {
          for (j in 1:m) {
            s_n[j] <- (1/(n-1)*sum((data[,j] - mean(data[,j]))^2))/K
          }
          mat <- matrix(0,m,m)
          diag(mat) <- s_n
          sigma[[k]] <- mat
        }
      }
      #if(var_fun == 3){
        #v3<-em.diag(data,K,steps=steps+20,var_fun = 3)
       # sigma <- v3$sigma
      #}
      #if(var_fun == 4){
       # v3<-em.diag(data,K,steps=steps+20,var_fun = 4)
      #  sigma <- v3$sigma
     # }
    }
  }

  if(m ==1){
    if (missing(p)){
      p <- rep(1/K, K)
    }
    if (missing(z)){
      z <- rnorm(K,mean = 0, sd=1)
    }
    if (missing(alpha)){
      alpha <- mean(data)
    }
    if (missing(beta)){
      beta <- (sample(data, 1) - alpha)
    }
    if (missing(sigma)){
      sigma <- sd(data)
    }
  }
  ##################################

  s<-1
  while (s <=steps){
    if(K == 1){
      W <- matrix(1,n,1)
      p <- 1
      z <- 0
      beta <- 0
      alpha <- colMeans(data)
      sigma <- colSds(as.matrix(data))
    } else {
      W   <- estep(data,p,z,alpha,beta,sigma,var_fun)
      fit <- mstep(data, W, alpha=alpha, beta=beta,sigma=sigma,var_fun)
      p   <- fit$p
      z  <- fit$z
      beta <- fit$beta
      alpha <- fit$alpha
      sigma <-fit$sigma

    }
    s<-s+1
    setTxtProgressBar(pb,s)
  }

  if(K != 1 && m != 1 && is.na(beta[1]) == FALSE){
    for (l in 1:length(beta)) {
      if(beta[l] == 0 & beta[l+1]<0){
        beta_hat <- beta*(-1)
        z_hat <- z*(-1)
      }
      if(beta[1] < 0){
        beta_hat <- beta*(-1)
        z_hat <- z*(-1)
      }
      if(beta[1] > 0){
        beta_hat <- beta
        z_hat <- z
      }
    }
  }
  if(m == 1 && is.na(beta) == FALSE){
    if(beta < 0)
      beta_hat <- beta*(-1)
    z_hat <- z*(-1)
    if(beta > 0){
      beta_hat <- beta
      z_hat <- z
    }
  }
  if(m != 1 && is.na(beta[1]) != FALSE){
    beta_hat <- beta
    z_hat <- z
  }
  if(m == 1 && is.na(beta) != FALSE){
    beta_hat <- beta
    z_hat <- z
  }
  if(K == 1 && is.na(beta) == FALSE && is.na(z) == FALSE){
    beta_hat <- beta
    z_hat <- z
  }


  if(m != 1){
    if(var_fun == 1){
      m <- length(data[1,])
      dens.all<-matrix(0,length(data[,1]), K)
      for (k in 1:K){

        dens.all[,k] <- p[k]*dmvnorm(x=data, mean=alpha + beta*z[k], sigma=diag(sigma^2,m,m),log=FALSE)
      }
    }
    if(var_fun == 2){
      m <- length(data[1,])
      dens.all<-matrix(0,length(data[,1]), K)
      for (k in 1:K){

        dens.all[,k] <- p[k]*dmvnorm(x=data, mean=alpha + beta*z[k],sigma=diag(unlist(sigma[k])^2,m,m), log=FALSE)
      }
    }
    if(var_fun == 3){
      m <- length(data[1,])
      dens.all<-matrix(0,length(data[,1]), K)
      for (k in 1:K){

        dens.all[,k] <- p[k]*dmvnorm(x=data, mean=alpha + beta*z[k], sigma=matrix(sigma,m,m),log=FALSE)
      }
    }
    if(var_fun == 4){
      m <- length(data[1,])
      dens.all<-matrix(0,length(data[,1]), K)
      for (k in 1:K){

        dens.all[,k] <- p[k]*dmvnorm(x=data, mean=alpha + beta*z[k],sigma=matrix(unlist(sigma[k]),m,m),log=FALSE)
      }
    }
  }
  if(m == 1){
    m = 1
    dens.all<-matrix(0,length(data), K)
    for (k in 1:K){

      dens.all[,k] <- dnorm(x=data, mean=alpha + beta*z[k], sd=sigma)*p[k]
    }
  }


  loglik<- sum(log(apply(dens.all,1,sum)))


  if(var_fun == 1){
    num_parameters <- (K-1) + K + m + m + m }
  if(var_fun == 2){
    num_parameters <- (K-1) + K + m + m + m*K
  }
  if(var_fun == 3){
    num_parameters <- (K-1) + K + m + m + m*(m+1)/2
  }
  if(var_fun == 4){
    num_parameters <- (K-1) + K + m + m + m*(m+1)/2*K
  }


  if (K!=1){
    n <- dim(W)[1]
    K <- dim(W)[2]

    wz_hat <- matrix(0,n,K)
    for (k in 1:K) {
      wz_hat[,k] <- W[,k]*z_hat[k]
    }
    x_proj <- apply(wz_hat,1,sum)
  }
  if(K==1){
    x_proj <- data
  }

  fit<- list("p"=p, "alpha"=alpha,"z"=z,"beta"=beta,"z_hat"=z_hat,"beta_hat"=beta_hat ,"sigma"=sigma,
             "W" =W, "loglikelihood"=loglik, "disparity"=-2*loglik, "number_parameters"=num_parameters,
             "AIC"=-2*loglik+2*num_parameters, "projection" = x_proj, "BIC"= -2*loglik + num_parameters*log(length(data[,1])) )
  class(fit) <- "EM"
  return(fit)
}




####em_covs.diag and em_covs.start(they have the same estep_covs)
estep_covs<- function(data, v, p, z, alpha, beta, gamma, sigma, var_fun){
  if(missing(var_fun)){
    var_fun <- 2
  }

  m <- length(data[1,])

  if(m != 1){
    data <- data.frame(data)
    v <- as.data.frame(v)
    n <- length(data[,1])
    K <- length(p) #
    m <- length(data[1,]) #

    mu <- matrix(0,K,m)
    for (k in 1:K) {
      mu[k,] <- alpha + beta*z[k]
    }

    data_new <- data.frame()
    for (i in 1:n) {
      data_new.1 <- data[i,] - gamma%*%as.numeric(v[i,])
      data_new <- rbind(data_new,data_new.1)
    }
    if(var_fun == 1){
      W_1 <- matrix(0,n,K)
      for (k in 1:K) {
        W_1[,k] <- p[k]*dmvnorm(x=data_new, mean=mu[k,], sigma=diag(sigma^2,m,m), log=FALSE)
      }
    }
    if(var_fun == 2){
      W_1 <- matrix(0,n,K)
      for (k in 1:K) {
        W_1[,k] <- p[k]*dmvnorm(x=data_new, mean=mu[k,], sigma=diag(unlist(sigma[k])^2,m,m), log=FALSE)
      }
    }
    if(var_fun == 3){
      W_1 <- matrix(0,n,K)
      for (k in 1:K) {
        W_1[,k] <- p[k]*dmvnorm(x=data_new, mean=mu[k,], sigma=matrix(sigma,m,m), log=FALSE)
      }
    }
    if(var_fun == 4){
      W_1 <- matrix(0,n,K)
      for (k in 1:K) {
        W_1[,k] <- p[k]*dmvnorm(x=data_new, mean=mu[k,], sigma=matrix(unlist(sigma[k]),m,m), log=FALSE)
      }
    }

    W <- W_1/apply(W_1,1,sum)

    return(W)
  }

  if(m == 1){
    data <- data.frame(as.numeric(unlist(data)))
    v <- data.frame(v)

    n <- length(data[,1])
    K <- length(p)
    mu <- matrix(0,K,m)
    for (k in 1:K) {
      mu[k,] <- alpha + beta*z[k]
    }

    #data_new <- data - gamma*v
    data_new <- data.frame()
    for (i in 1:n) {
      data_new.1 <- data[i,] - gamma%*%as.numeric(v[i,])
      data_new <- rbind(data_new,data_new.1)
    }

    W_1 <- matrix(0,n,K)
    for (k in 1:K) {
      W_1[,k] <- p[k]*dnorm(x=as.numeric(unlist(data_new)), mean=mu[k,], sd=sqrt(sigma[k]), log=FALSE)
    }
    W <- W_1/apply(W_1,1,sum)
    return(W)
  }

}


mstep_covs.start <- function(data, v, W, alpha, beta, gamma, var_fun){
  if(missing(var_fun)){
    var_fun <- 2
  }
  m <- length(data[1,])

  if(m != 1){
    data <- as.data.frame(data)
    v <- as.data.frame(v)
    m <- length(data[1,])
    n <- dim(W)[1]
    K <- dim(W)[2]

    if(var_fun == 1){
      p <- apply(W,2,sum)/n

      counter <- 0
      while (counter <= 10) {
        z_1 <- matrix(0,K,m)
        for (k in 1:K) {
          z_1[k,] <- apply(W[,k]*data,2,sum)
        }
        z_2 <- matrix(0,K,m)
        for (j in 1:m) {
          z_2[,j]  <- beta[j]*z_1[,j]
        }
        z_3 <- apply(z_2,1,sum) ##

        z_4 <- matrix(0,K,m)
        for (j in 1:m) {
          z_4[,j] <- alpha[j]*beta[j]*apply(W,2,sum)
        }
        z_5 <- apply(z_4,1,sum)##

        ##part 3
        z_6 <- matrix(0,K,m)
        for (j in 1:m) {
          z_6[,j] <- ((beta[j])^2)*apply(W,2,sum)
        }
        z_7 <- apply(z_6,1,sum) ##

        ##part 4
        phi <- data.frame()
        for (i in 1:n) {
          phi_current <- t(gamma %*% t(as.matrix(v[i,])))
          phi <- rbind(phi,phi_current)
        }

        z_8 <- matrix(0,K,m)
        for (k in 1:k) {
          z_8[k,] <- apply(W[,k]*phi,2,sum)
        }
        z_9 <- matrix(0,K,m)
        for (j in 1:m) {
          z_9[,j]  <- beta[j]*z_8[,j]
        }
        z_10 <- apply(z_9,1,sum) ##

        z_new <- (z_3 - z_5 - z_10)/z_7
        z_j <- z_new - sum(z_new*p)
        z <- z_j/sqrt(sum((z_j^2)*p))

        ####alpha
        wz <- matrix(0,n,K)
        for (k in 1:K) {
          wz[,k] <- W[,k]*z[k]
        }
        alpha <- (1/n)*(apply(data,2,sum) - beta*sum(apply(wz,1,sum))  - apply(phi,2,sum))

        #####beta
        wz <- matrix(0,n,K)
        for (k in 1:K) {
          wz[,k] <- W[,k]*z[k]
        }

        wz_2 <- matrix(0,n,K)
        for (k in 1:K) {
          wz_2[,k] <- W[,k]*(z[k])^2
        }

        beta_1 <- apply(data*apply(wz,1,sum),2,sum)

        beta_v_1 <- apply(phi*apply(wz,1,sum),2,sum)

        beta_v_2 <- (1/n)*sum(apply(wz,1,sum))* apply(phi,2,sum)

        beta <- (beta_1 - (1/n)*(apply(data,2,sum)*sum(apply(wz,1,sum))) - beta_v_1 + beta_v_2)/(sum(apply(wz_2,1,sum)) - (1/n)*(sum(apply(wz,1,sum)))^2)

        ####gamma
        v_current <- vector("list", length = n)
        for (i in 1:n) {
          v_current[[i]] <- as.numeric(v[i,]) %*% t(as.numeric(v[i,]))
        }
        v_sum <- Reduce("+", v_current)
        v_sum_inverse <- solve(v_sum)

        current_gamma  <- vector("list", length = n)
        Gamma_1 <- vector("list",length = K)
        for (k in 1:K) {
          for (i in 1:n) {
            current_gamma[[i]] <- W[i,k]*(as.numeric(data[i,] - alpha - beta*z[k])%*%t(as.numeric(v[i,])))
            sum_matrix <- Reduce("+", current_gamma)
          }
          Gamma_1[[k]] <- sum_matrix
          gamma_new <- Reduce("+", Gamma_1)
        }
        gamma <- gamma_new%*%v_sum_inverse

        counter = counter+1
      }
      diff<-matrix(0,n, K)
      sigma <- c()
      for (j in 1:m) {
        for (k in 1:K) {
          diff[,k] <- (data[,j] - alpha[j] - beta[j]*z[k] - phi[,j])^2
        }
        sigma[j] <- sqrt(mean(apply(W*diff,1,sum)))
      }
    }

    if(var_fun == 2){
      p <- apply(W,2,sum)/n

      counter <- 0
      while (counter <= 10) {
        z_1 <- matrix(0,K,m)
        for (k in 1:K) {
          z_1[k,] <- apply(W[,k]*data,2,sum)
        }
        z_2 <- matrix(0,K,m)
        for (j in 1:m) {
          z_2[,j]  <- beta[j]*z_1[,j]
        }
        z_3 <- apply(z_2,1,sum) ##

        z_4 <- matrix(0,K,m)
        for (j in 1:m) {
          z_4[,j] <- alpha[j]*beta[j]*apply(W,2,sum)
        }
        z_5 <- apply(z_4,1,sum)

        ##part 3
        z_6 <- matrix(0,K,m)
        for (j in 1:m) {
          z_6[,j] <- ((beta[j])^2)*apply(W,2,sum)
        }
        z_7 <- apply(z_6,1,sum) ##

        ##part 4
        phi <- data.frame()
        for (i in 1:n) {
          phi_current <- t(gamma %*% t(as.matrix(v[i,])))
          phi <- rbind(phi,phi_current)
        }

        z_8 <- matrix(0,K,m)
        for (k in 1:k) {
          z_8[k,] <- apply(W[,k]*phi,2,sum)
        }
        z_9 <- matrix(0,K,m)
        for (j in 1:m) {
          z_9[,j]  <- beta[j]*z_8[,j]
        }
        z_10 <- apply(z_9,1,sum) ##

        z_new <- (z_3 - z_5 - z_10)/z_7
        z_j <- z_new - sum(z_new*p)
        z <- z_j/sqrt(sum((z_j^2)*p))

        ####alpha
        wz <- matrix(0,n,K)
        for (k in 1:K) {
          wz[,k] <- W[,k]*z[k]
        }
        alpha <- (1/n)*(apply(data,2,sum) - beta*sum(apply(wz,1,sum))  - apply(phi,2,sum))

        #####beta
        wz <- matrix(0,n,K)
        for (k in 1:K) {
          wz[,k] <- W[,k]*z[k]
        }

        wz_2 <- matrix(0,n,K)
        for (k in 1:K) {
          wz_2[,k] <- W[,k]*(z[k])^2
        }

        beta_1 <- apply(data*apply(wz,1,sum),2,sum)

        beta_v_1 <- apply(phi*apply(wz,1,sum),2,sum)

        beta_v_2 <- (1/n)*sum(apply(wz,1,sum))* apply(phi,2,sum)

        beta <- (beta_1 - (1/n)*(apply(data,2,sum)*sum(apply(wz,1,sum))) - beta_v_1 + beta_v_2)/(sum(apply(wz_2,1,sum)) - (1/n)*(sum(apply(wz,1,sum)))^2)

        ####gamma
        v_current <- vector("list", length = n)
        for (i in 1:n) {
          v_current[[i]] <- as.numeric(v[i,]) %*% t(as.numeric(v[i,]))
        }
        v_sum <- Reduce("+", v_current)
        v_sum_inverse <- solve(v_sum)

        current_gamma  <- vector("list", length = n)
        Gamma_1 <- vector("list",length = K)
        for (k in 1:K) {
          for (i in 1:n) {
            current_gamma[[i]] <- W[i,k]*(as.numeric(data[i,] - alpha - beta*z[k])%*%t(as.numeric(v[i,])))
            sum_matrix <- Reduce("+", current_gamma)
          }
          Gamma_1[[k]] <- sum_matrix
          gamma_new <- Reduce("+", Gamma_1)
        }
        gamma <- gamma_new%*%v_sum_inverse

        counter = counter+1
      }

      diff<- matrix(0,n,m)
      sigma <- c()
      for (k in 1:K) {
        for (j in 1:m) {
          diff[,j] <- (data[,j] - alpha[j] - beta[j]*z[k] - phi[,j])^2
        }
        sigma[[k]] <- sqrt(apply(W[,k]*diff,2,sum)/sum(W[,k]))
      }
    }

    if (var_fun == 3){
      p <- apply(W,2,sum)/n

      counter <- 0
      while (counter <= 10) {
        z_1 <- matrix(0,K,m)
        for (k in 1:K) {
          z_1[k,] <- apply(W[,k]*data,2,sum)
        }
        z_2 <- matrix(0,K,m)
        for (j in 1:m) {
          z_2[,j]  <- beta[j]*z_1[,j]
        }
        z_3 <- apply(z_2,1,sum) ##

        z_4 <- matrix(0,K,m)
        for (j in 1:m) {
          z_4[,j] <- alpha[j]*beta[j]*apply(W,2,sum)
        }
        z_5 <- apply(z_4,1,sum)##

        ##part 3
        z_6 <- matrix(0,K,m)
        for (j in 1:m) {
          z_6[,j] <- ((beta[j])^2)*apply(W,2,sum)
        }
        z_7 <- apply(z_6,1,sum) ##

        ##part 4
        phi <- data.frame()
        for (i in 1:n) {
          phi_current <- t(gamma %*% t(as.matrix(v[i,])))
          phi <- rbind(phi,phi_current)
        }

        z_8 <- matrix(0,K,m)
        for (k in 1:k) {
          z_8[k,] <- apply(W[,k]*phi,2,sum)
        }
        z_9 <- matrix(0,K,m)
        for (j in 1:m) {
          z_9[,j]  <- beta[j]*z_8[,j]
        }
        z_10 <- apply(z_9,1,sum) ##

        z_new <- (z_3 - z_5 - z_10)/z_7
        z_j <- z_new - sum(z_new*p)
        z <- z_j/sqrt(sum((z_j^2)*p))

        #####
        phi <- data.frame()
        for (i in 1:n) {
          phi_current <- t(gamma %*% t(as.matrix(v[i,])))
          phi <- rbind(phi,phi_current)
        }
        ####alpha
        wz <- matrix(0,n,K)
        for (k in 1:K) {
          wz[,k] <- W[,k]*z[k]
        }
        alpha <- (1/n)*(apply(data,2,sum) - beta*sum(apply(wz,1,sum))  - apply(phi,2,sum))

        #####beta
        wz <- matrix(0,n,K)
        for (k in 1:K) {
          wz[,k] <- W[,k]*z[k]
        }

        wz_2 <- matrix(0,n,K)
        for (k in 1:K) {
          wz_2[,k] <- W[,k]*(z[k])^2
        }

        beta_1 <- apply(data*apply(wz,1,sum),2,sum)

        beta_v_1 <- apply(phi*apply(wz,1,sum),2,sum)

        beta_v_2 <- (1/n)*sum(apply(wz,1,sum))* apply(phi,2,sum)

        beta <- (beta_1 - (1/n)*(apply(data,2,sum)*sum(apply(wz,1,sum))) - beta_v_1 + beta_v_2)/(sum(apply(wz_2,1,sum)) - (1/n)*(sum(apply(wz,1,sum)))^2)

        ####gamma
        v_current <- vector("list", length = n)
        for (i in 1:n) {
          v_current[[i]] <- as.numeric(v[i,]) %*% t(as.numeric(v[i,]))
        }
        v_sum <- Reduce("+", v_current)
        v_sum_inverse <- solve(v_sum)

        current_gamma  <- vector("list", length = n)
        Gamma_1 <- vector("list",length = K)
        for (k in 1:K) {
          for (i in 1:n) {
            current_gamma[[i]] <- W[i,k]*(as.numeric(data[i,] - alpha - beta*z[k])%*%t(as.numeric(v[i,])))
            sum_matrix <- Reduce("+", current_gamma)
          }
          Gamma_1[[k]] <- sum_matrix
          gamma_new <- Reduce("+", Gamma_1)
        }
        gamma <- gamma_new%*%v_sum_inverse

        counter = counter+1
      }
      current  <- vector("list", length = n)
      Sigma_1 <- vector("list",length = K)
      for (k in 1:K) {
        for (i in 1:n) {
          current[[i]] <- W[i,k]*t(as.matrix(data[i,] - alpha - beta*z[k] - phi[i,]))%*%(as.matrix(data[i,] - alpha - beta*z[k] - phi[i,]))
          sum_matrix <- Reduce("+", current)
        }
        Sigma_1[[k]] <- sum_matrix/n
        sigma <- Reduce("+", Sigma_1)
      }
    }

    if (var_fun == 4){
      p <- apply(W,2,sum)/n

      counter <- 0
      while (counter <= 10) {
       # z_1 <- matrix(0,n,K)
      #  for (i in 1:n) {
       #   z_1[i,] <- W[i,]*as.numeric(t(beta) %*% solve(matrix(unlist(sigma),m,m)) %*% t(data[i,] - alpha -  as.numeric(gamma*v[i,])))
      #  }
       # z_2 <- apply(z_1,2,sum)

        #z_3 <- matrix(0,n,K)
        #for (k in 1:K) {
         # z_3[,k] <- W[,k]*as.numeric((t(matrix(beta)) %*% solve(matrix(unlist(sigma),m,m)) %*% matrix(beta)))
        #}
        #z_4 <- apply(z_3,2,sum)

        #z_new <- z_2/z_4
        #z_j <- z_new - sum(z_new*p)
        #z <- z_j/sqrt(sum((z_j^2)*p))
        z_1 <- matrix(0,K,m)
        for (k in 1:K) {
          z_1[k,] <- apply(W[,k]*data,2,sum)
        }
        z_2 <- matrix(0,K,m)
        for (j in 1:m) {
          z_2[,j]  <- beta[j]*z_1[,j]
        }
        z_3 <- apply(z_2,1,sum) ##

        z_4 <- matrix(0,K,m)
        for (j in 1:m) {
          z_4[,j] <- alpha[j]*beta[j]*apply(W,2,sum)
        }
        z_5 <- apply(z_4,1,sum)##

        ##part 3
        z_6 <- matrix(0,K,m)
        for (j in 1:m) {
          z_6[,j] <- ((beta[j])^2)*apply(W,2,sum)
        }
        z_7 <- apply(z_6,1,sum) ##

        ##part 4
        phi <- data.frame()
        for (i in 1:n) {
          phi_current <- t(gamma %*% t(as.matrix(v[i,])))
          phi <- rbind(phi,phi_current)
        }

        z_8 <- matrix(0,K,m)
        for (k in 1:k) {
          z_8[k,] <- apply(W[,k]*phi,2,sum)
        }
        z_9 <- matrix(0,K,m)
        for (j in 1:m) {
          z_9[,j]  <- beta[j]*z_8[,j]
        }
        z_10 <- apply(z_9,1,sum) ##

        z_new <- (z_3 - z_5 - z_10)/z_7
        z_j <- z_new - sum(z_new*p)
        z <- z_j/sqrt(sum((z_j^2)*p))

        #####
        phi <- data.frame()
        for (i in 1:n) {
          phi_current <- t(gamma %*% t(as.matrix(v[i,])))
          phi <- rbind(phi,phi_current)
        }
        ####alpha
        wz <- matrix(0,n,K)
        for (k in 1:K) {
          wz[,k] <- W[,k]*z[k]
        }
        alpha <- (1/n)*(apply(data,2,sum) - beta*sum(apply(wz,1,sum))  - apply(phi,2,sum))

        #####beta
        wz <- matrix(0,n,K)
        for (k in 1:K) {
          wz[,k] <- W[,k]*z[k]
        }

        wz_2 <- matrix(0,n,K)
        for (k in 1:K) {
          wz_2[,k] <- W[,k]*(z[k])^2
        }

        beta_1 <- apply(data*apply(wz,1,sum),2,sum)

        beta_v_1 <- apply(phi*apply(wz,1,sum),2,sum)

        beta_v_2 <- (1/n)*sum(apply(wz,1,sum))* apply(phi,2,sum)

        beta <- (beta_1 - (1/n)*(apply(data,2,sum)*sum(apply(wz,1,sum))) - beta_v_1 + beta_v_2)/(sum(apply(wz_2,1,sum)) - (1/n)*(sum(apply(wz,1,sum)))^2)

        ####gamma
        v_current <- vector("list", length = n)
        for (i in 1:n) {
          v_current[[i]] <- as.numeric(v[i,]) %*% t(as.numeric(v[i,]))
        }
        v_sum <- Reduce("+", v_current)
        v_sum_inverse <- solve(v_sum)

        current_gamma  <- vector("list", length = n)
        Gamma_1 <- vector("list",length = K)
        for (k in 1:K) {
          for (i in 1:n) {
            current_gamma[[i]] <- W[i,k]*(as.numeric(data[i,] - alpha - beta*z[k])%*%t(as.numeric(v[i,])))
            sum_matrix <- Reduce("+", current_gamma)
          }
          Gamma_1[[k]] <- sum_matrix
          gamma_new <- Reduce("+", Gamma_1)
        }
        gamma <- gamma_new%*%v_sum_inverse

        counter = counter+1
      }
      ##sigma
      current  <- vector("list", length = n)
      sigma <- vector("list",length = K)
      for (k in 1:K) {
        for (i in 1:n) {
          current[[i]] <- W[i,k]*t(as.matrix(data[i,] - alpha - beta*z[k] - phi[i,] ))%*%(as.matrix(data[i,] - alpha - beta*z[k] - phi[i,] ))
          add <- function(x) Reduce("+", x)
          sum_matrix <- add(current)
        }
        sigma[[k]] <- sum_matrix/sum(W[,k])
      }
    }

  }

  if (m == 1){
    data <- data.frame(as.numeric(unlist(data)))
    v <- as.data.frame(v)
    m <- length(data[1,])
    n <- dim(W)[1]
    K <- dim(W)[2]

    p <- apply(W,2,sum)/n
    counter <- 0
    while (counter <= 10) {

      z_1_pre <- data.frame()
      for (i in 1:n) {
        data_new.1 <- data[i,] - alpha - gamma%*%as.numeric(v[i,])
        z_1_pre <- rbind(z_1_pre, data_new.1)
      }

      z_1 <- matrix(0,K,m)
      for (k in 1:K) {
        z_1[k,] <- sum(W[,k]*z_1_pre)
      }

      z_2 <- c()
      for (j in 1:m) {
        z_2 <- (beta)*apply(W,2,sum)
      }
      z_new <- z_1/z_2
      z_j <- z_new - sum(z_new*p)
      z <- z_j/sqrt(sum((z_j^2)*p))

      ##beta
      phi <- data.frame()
      for (i in 1:n) {
        phi_current <- gamma%*%as.numeric(v[i,])
        phi <- rbind(phi,phi_current)
      }

      wz <- matrix(0,n,K)
      for (k in 1:K) {
        wz[,k] <- W[,k]*z[k]
      }

      wz_2 <- matrix(0,n,K)
      for (k in 1:K) {
        wz_2[,k] <- W[,k]*(z[k])^2
      }

      beta_1 <- apply(data*apply(wz,1,sum),2,sum)

      beta_v_1 <- apply(phi*apply(wz,1,sum),2,sum)

      beta_v_2 <- (1/n)*sum(apply(wz,1,sum))* apply(phi,2,sum)

      beta <- (beta_1 - (1/n)*(apply(data,2,sum)*sum(apply(wz,1,sum))) - beta_v_1 + beta_v_2)/(sum(apply(wz_2,1,sum)) - (1/n)*(sum(apply(wz,1,sum)))^2)

      ##alpha
      wz <- matrix(0,n,K)
      for (k in 1:K) {
        wz[,k] <- W[,k]*z[k]
      }
      alpha <- (1/n)*(apply(data,2,sum) - beta*sum(apply(wz,1,sum))  - apply(phi,2,sum))

      #gamma
      v_current <- vector("list", length = n)
      for (i in 1:n) {
        v_current[[i]] <- as.numeric(v[i,]) %*% t(as.numeric(v[i,]))
      }
      v_sum <- Reduce("+", v_current)
      v_sum_inverse <- solve(v_sum)

      current_gamma  <- vector("list", length = n)
      Gamma_1 <- vector("list",length = K)
      for (k in 1:K) {
        for (i in 1:n) {
          current_gamma[[i]] <- W[i,k]*(as.numeric(data[i,] - alpha - beta*z[k])%*%t(as.numeric(v[i,])))
          sum_matrix <- Reduce("+", current_gamma)
        }
        Gamma_1[[k]] <- sum_matrix
        gamma_new <- Reduce("+", Gamma_1)
      }
      gamma <- gamma_new%*%v_sum_inverse

      counter <- counter + 1
    }

    diff<- matrix(0,n,1)
    sigma <- c()
    for (k in 1:K) {
      diff <- (as.numeric(unlist(data)) - alpha - beta*z[k] - as.numeric(unlist(phi)))^2
      sigma[k] <- sqrt(sum(W[,k]*diff)/sum(W[,k]))
    }

  }


  return(list("p"=p, "z"=z, "alpha"=alpha, "beta"=beta, "gamma"=gamma,"sigma"=sigma))
}



em_covs.start <-  function(data, v, K, p, alpha, z, beta, sigma, gamma, steps, var_fun){
  pb<-txtProgressBar(min=0, max=steps, style=3)

  data <- as.data.frame(data)
  v <- as.data.frame(v)
  q <- as.numeric(length(v[1,]))

  if(missing(var_fun)){
    var_fun <- 2
  }
  m <- as.numeric(length(data[1,]))
  if(m != 1){
    m <- length(data[1,])
    n <- length(data[,1])
  }
  if(m == 1){
    n <- length(data)
  }

  if(m != 1){
    if(missing(gamma)){
      m <- length(data[1,])
      q <- length(v[1,])
      dep_vars <- colnames(data)
      ind_vars <- colnames(v)
      var_comb <- expand.grid(dep_vars, ind_vars )
      formula_vec <- sprintf("%s ~ %s", var_comb$Var1, var_comb$Var2)
      lm_res <- lapply( formula_vec, function(f) {
        fit1 <- lm( f, data = data.frame(data,v))
        return(fit1$coefficients[2])
      })
      gamma_new <- as.numeric(unlist(lm_res))
      gamma <- matrix(gamma_new,m,q)
    }
    if (missing(p)){
      p <- rep(1/K, K)
    }
    if (missing(z)){
      z <- rnorm(K,mean = 0, sd=1)
    }
    if (missing(alpha)){
      alpha <- colMeans(data)
    }
    if (missing(beta)){
      beta <- as.numeric(data[sample(nrow(data), 1), ] - alpha)
    }
    if (missing(sigma)){
      if(var_fun == 1){
        s_n <- c()
        for (j in 1:m) {
          s_n[j] <- sqrt(1/(n-1)*sum((data[,j] - mean(data[,j]))^2))
        }
        sigma <- s_n/K
      }
      if(var_fun == 2){
        s_n <- c()
        sigma <- list()
        for (k in 1:K) {
          for (j in 1:m) {
            s_n[j] <- sqrt(1/(n-1)*sum((data[,j] - mean(data[,j]))^2))/K
          }
          sigma[[k]] <- s_n}
      }
      if(var_fun == 3){
        s_n <- c()
        sigma <- matrix(0,m,m)
        for (k in 1:K) {
          for (j in 1:m) {
           # s_n[j] <- (1/(n-1)*sum((data[,j] - mean(data[,j]))^2))/K
            s_n[j] <- (1/(n-1)*sum((data[,j] - mean(data[,j]))^2))
          }
          mat <- matrix(0,m,m)
          diag(mat) <- s_n
          sigma <- mat
        }
      }
      if(var_fun == 4){
        s_n <- c()
        sigma <- vector("list",length = K)
        for (k in 1:K) {
          for (j in 1:m) {
            s_n[j] <- (1/(n-1)*sum((data[,j] - mean(data[,j]))^2))/K
          }
          mat <- matrix(0,m,m)
          diag(mat) <- s_n
          sigma[[k]] <- mat
        }
      }
    }
  }

  if(m == 1){
    if (missing(p)){
      p <- rep(1/K, K)
    }
    if (missing(z)){
      z <- rnorm(K,mean = 0, sd=1)
    }
    if (missing(alpha)){
      alpha <- mean(unlist(data))
    }
    if (missing(beta)){
      beta <- (sample(unlist(data), 1) - alpha)
    }
    if (missing(sigma)){
      sigma <- sd(unlist(data))
    }
    if (missing(gamma)){
      m <- as.numeric(length(data[1,]))
      q <- as.numeric(length(v[1,]))
      dep_vars <- colnames(data)
      ind_vars <- colnames(v)
      var_comb <- expand.grid(dep_vars, ind_vars )
      formula_vec <- sprintf("%s ~ %s", var_comb$Var1, var_comb$Var2)
      lm_res <- lapply( formula_vec, function(f) {
        fit1 <- lm( f, data = data.frame(data,v))
        return(fit1$coefficients[2])
      })
      gamma_new <- as.numeric(unlist(lm_res))
      gamma <- matrix(gamma_new,m,q)
    }
  }


  s<-1
  while (s <=steps){
    W   <- estep_covs(data,v,p,z,alpha,beta,gamma,sigma,var_fun)

    fit <- mstep_covs.start(data, v, W, alpha=alpha, beta=beta, gamma=gamma,var_fun)
    p   <- fit$p
    z  <- fit$z
    beta <- fit$beta
    alpha <- fit$alpha
    gamma <- fit$gamma
    sigma <-fit$sigma
    s<-s+1
    setTxtProgressBar(pb,s)
  }
  #################################
  if (K != 1 && m != 1 && is.na(beta[1]) == FALSE) {
    for (l in 1:length(beta)) {
      if(beta[l] == 0 & beta[l+1]<0){
        beta_hat <- beta*(-1)
        z_hat <- z*(-1)
      }
      if(beta[1] < 0){
        beta_hat <- beta*(-1)
        z_hat <- z*(-1)
      }
      if(beta[1] > 0){
        beta_hat <- beta
        z_hat <- z
      }
    }
  }
  if(m != 1 && is.na(beta[1]) != FALSE){
    beta_hat <- beta
    z_hat <- z
  }

  if(m == 1 && is.na(beta) == FALSE){
    if(beta < 0)
      beta_hat <- beta*(-1)
    z_hat <- z*(-1)
    if(beta > 0){
      beta_hat <- beta
      z_hat <- z
    }
  }

  if(m == 1 && is.na(beta[1]) != FALSE){
    beta_hat <- beta
    z_hat <- z
  }


  if(m != 1){
    data_new <- data.frame()
    for (i in 1:n) {
      data_new.1 <- data[i,] - gamma%*%as.numeric(v[i,])
      data_new <- rbind(data_new,data_new.1)
    }
    if(var_fun == 1){
      m <- length(data[1,])
      dens.all<-matrix(0,length(data[,1]), K)
      for (k in 1:K){
        dens.all[,k]<- p[k]*dmvnorm(x=data_new, mean=alpha + beta*z[k], sigma=diag(sigma^2,m,m), log=FALSE)
      }
    }
    if(var_fun == 2){
      m <- length(data[1,])
      dens.all<-matrix(0,length(data[,1]), K)
      for (k in 1:K){
        dens.all[,k]<- p[k]*dmvnorm(x=data_new, mean=alpha + beta*z[k], sigma=diag(unlist(sigma[k])^2,m,m), log=FALSE)
      }
    }
    if(var_fun == 3){
      m <- length(data[1,])
      dens.all<-matrix(0,length(data[,1]), K)
      for (k in 1:K){
        dens.all[,k]<- p[k]*dmvnorm(x=data_new, mean=alpha + beta*z[k], sigma=matrix(sigma,m,m), log=FALSE)
      }
    }
    if(var_fun == 4){
      m <- length(data[1,])
      dens.all<-matrix(0,length(data[,1]), K)
      for (k in 1:K){
        dens.all[,k]<- p[k]*dmvnorm(x=data_new, mean=alpha + beta*z[k], sigma=matrix(unlist(sigma[k]),m,m), log=FALSE)
      }
    }
  }

  if(m == 1){
    n <- length(data[,1])
    K <- length(p)
    mu <- matrix(0,K,m)
    for (k in 1:K) {
      mu[k,] <- alpha + beta*z[k]
    }

    data_new <- data.frame()
    for (i in 1:n) {
      data_new.1 <- data[i,] - gamma%*%as.numeric(v[i,])
      data_new <- rbind(data_new,data_new.1)
    }

    dens.all <- matrix(0,n,K)
    for (k in 1:K) {
      dens.all[,k] <- p[k]*dnorm(x=as.numeric(unlist(data_new)), mean=mu[k,], sd=sqrt(sigma[k]), log=FALSE)
    }
  }

  loglik<- sum(log(apply(dens.all,1,sum)))

  if(m != 1){
    if(var_fun == 1){
      num_parameters <- (K-1) + K + m + m + m + m*q
    }
    if(var_fun == 2){
      num_parameters <- (K-1) + K + m + m + m*K + m*q
    }
    if(var_fun == 3){
      num_parameters <- (K-1) + K + m + m + m*(m+1)/2 + m*q
    }
    if(var_fun == 4){
      num_parameters <- (K-1) + K + m + m + m*(m+1)/2*K + m*q
    }
  }

  if(m == 1){
    num_parameters <- (K-1) + K + m + m + m*K + m*q
  }

  fit<- list("p"=p, "z"=z_hat,"alpha"=alpha,"beta"=beta_hat, "gamma"=gamma,"sigma"=sigma,  "W" =W, "loglikelihood"=loglik, "disparity"=-2*loglik, "number_parameters"=num_parameters, "AIC"=-2*loglik+2*num_parameters)#,"BIC"= -2*loglik + num_parameters*log(length(data[,1])))

  class(fit) <- "EM"
  return(fit)
}




####em_2level_fun.start().
estep_mlevel_se<- function(data, p, z, alpha, beta, sigma, var_fun){

  if(missing(var_fun)){
    var_fun <- 2
  }
  data_numeric <- data[sapply(data, is.numeric)]
  data_factor <- as.factor(data[,names(Filter(is.factor,data))])

  m <- length(data_numeric[1,])
  K <- length(p)
  nr <- nlevels(data_factor)

  mu <- matrix(0,K,m)
  for (k in 1:K) {
    mu[k,] <- alpha + beta*z[k]
  }

  if(var_fun == 1) {
    X <- split(data_numeric, data_factor)
    m_ik <- data.frame()
    fik_j <- matrix(0,1,K)
    for (i in X) {
      for (k in 1:K) {
        for (j in 1:length(i[,1])) {
          fik_j[,k] <- prod(dmvnorm(x=data.frame(i)[j,], mean=mu[k,], sigma=diag(sigma^2,m,m), log=FALSE))
        }
      }
      m_ik <- rbind(m_ik,fik_j)
    }
    W_1 <- matrix(0,nr,K)
    for (k in 1:K) {
      W_1[,k] <- p[k]*m_ik[,k]
    }
  }

  if(var_fun == 2) {
    X <- split(data_numeric, data_factor)
    m_ik <- data.frame()
    fik_j <- matrix(0,1,K)
    for (i in X) {
      for (k in 1:K) {
        for (j in 1:length(i[,1])) {
          fik_j[,k] <- prod(dmvnorm(x=data.frame(i)[j,], mean=mu[k,], sigma=diag(unlist(sigma[k])^2,m,m), log=FALSE))
        }
      }
      m_ik <- rbind(m_ik,fik_j)
    }
    W_1 <- matrix(0,nr,K)
    for (k in 1:K) {
      W_1[,k] <- p[k]*m_ik[,k]
    }
  }

  if(var_fun == 3) {
    X <- split(data_numeric, data_factor)
    m_ik <- data.frame()
    fik_j <- matrix(0,1,K)
    for (i in X) {
      for (k in 1:K) {
        for (j in 1:length(i[,1])) {
          fik_j[,k] <- prod(dmvnorm(x=data.frame(i)[j,], mean=mu[k,], sigma=matrix(sigma,m,m), log=FALSE))
        }
      }
      m_ik <- rbind(m_ik,fik_j)
    }
    W_1 <- matrix(0,nr,K)
    for (k in 1:K) {
      W_1[,k] <- p[k]*m_ik[,k]
    }
  }

  if(var_fun == 4) {
    X <- split(data_numeric, data_factor)
    m_ik <- data.frame()
    fik_j <- matrix(0,1,K)
    for (i in X) {
      for (k in 1:K) {
        for (j in 1:length(i[,1])) {
          fik_j[,k] <- prod(dmvnorm(x=data.frame(i)[j,], mean=mu[k,], sigma= matrix(unlist(sigma[k]),m,m), log=FALSE))
        }
      }
      m_ik <- rbind(m_ik,fik_j)
    }
    W_1 <- matrix(0,nr,K)
    for (k in 1:K) {
      W_1[,k] <- p[k]*m_ik[,k]
    }
  }
  W <- W_1/apply(W_1,1,sum)
  return(W)

}



mstep_mlevel_se <- function(data, W, alpha, beta, sigma, var_fun){

  data_numeric <- data[sapply(data, is.numeric)]
  data_factor <- as.factor(data[,names(Filter(is.factor,data))])
  n_i_name_numeric <- as.data.frame(summary(data_factor,maxsum = length(unique(data_factor))))
  n_i <- as.vector(n_i_name_numeric[,1])

  if(missing(var_fun)){
    var_fun <- 2
  }

  n <- length(data_numeric[,1])
  m <- length(data_numeric[1,])
  K <- dim(W)[2]
  nr <- nlevels(data_factor)

  if(var_fun == 1){

    p <- apply(W,2,sum)/nr

    counter <- 0
    while (counter <= 5) {
      x_i <- data.frame()
      X <- split(data_numeric, data_factor)
      for (i in X) {
        for (j in 1:length(i[,1])) {
          xij <- colSums(data.frame(i))
        }
        x_i <- rbind(x_i,xij)
      }

      ##z
      z_1 <- matrix(0,K,m)
      for (k in 1:K) {
        z_1[k,] <- apply(W[,k]*x_i,2,sum)
      }
      z_2 <- matrix(0,K,m)
      for (j in 1:m) {
        z_2[,j]  <- beta[j]*z_1[,j]
      }
      z_3 <- apply(z_2,1,sum) ##
      #part 2
      z_4 <- matrix(0,K,m)
      for (j in 1:m) {
        z_4[,j] <- alpha[j]*beta[j]*apply(W*n_i,2,sum)
      }
      z_5 <- apply(z_4,1,sum)##

      ##part 3
      z_6 <- matrix(0,K,m)
      for (j in 1:m) {
        z_6[,j] <- ((beta[j])^2)*apply(W*n_i,2,sum)
      }
      z_7 <- apply(z_6,1,sum) ##

      z_new <- (z_3 - z_5)/z_7
      z_j <- z_new - sum(z_new*p)
      z <- z_j/sqrt(sum((z_j^2)*p))

      x_i <- data.frame()
      X <- split(data_numeric, data_factor)
      for (i in X) {
        for (j in 1:length(i[,1])) {
          xij <- colSums(data.frame(i))
        }
        x_i <- rbind(x_i,xij)
      }

      wz <- matrix(0,nr,K)
      for (k in 1:K) {
        wz[,k] <- W[,k]*z[k]
      }
      wz_2 <- matrix(0,nr,K)
      for (k in 1:K) {
        wz_2[,k] <- W[,k]*(z[k])^2
      }

      beta_1 <- apply(x_i*apply(wz,1,sum),2,sum)
      beta <- (beta_1 - (1/n)*(apply(data_numeric,2,sum)*sum(n_i*apply(wz,1,sum)))) /
        (sum(n_i*apply(wz_2,1,sum)) - (1/n)*(sum(n_i*apply(wz,1,sum)))^2)

      ##alpha
      alpha <- (1/n)*(apply(data_numeric,2,sum) - beta*sum(n_i*apply(wz,1,sum)))

      counter = counter+1
    }

    ##sigma
    diff<-matrix(0,n, K)
    sigma <- c()
    for (l in 1:m) {
      for (k in 1:K) {
        diff[,k] <- (data_numeric[,l] - alpha[l] - beta[l]*z[k])^2
      }
      diff_i <- data.frame()
      diff_df <- data.frame(cbind(diff,data_factor))
      Y <- split(diff_df, data_factor)
      for (i in Y) {
        for (j in 1:length(i[,1])) {
          diffij <- colSums(data.frame(i))
        }
        diff_i <- rbind(diff_i,diffij)
      }
      diff_i <- as.matrix(diff_i[,-length(diff_i[1,])])
      sigma[l] <- sqrt(mean(apply(W*diff_i,1,sum)))
    }
  }

  if(var_fun == 2){

    p <- apply(W,2,sum)/nr

    counter <- 0
    while (counter <= 5) {
      x_i <- data.frame()
      X <- split(data_numeric, data_factor)
      for (i in X) {
        for (j in 1:length(i[,1])) {
          xij <- colSums(data.frame(i))
        }
        x_i <- rbind(x_i,xij)
      }

      ###z
      z_1 <- matrix(0,K,m)
      for (k in 1:K) {
        z_1[k,] <- apply(W[,k]*x_i,2,sum)
      }
      z_2 <- matrix(0,K,m)
      for (j in 1:m) {
        z_2[,j]  <- beta[j]*z_1[,j]
      }
      z_3 <- apply(z_2,1,sum) ##
      #part 2
      z_4 <- matrix(0,K,m)
      for (j in 1:m) {
        z_4[,j] <- alpha[j]*beta[j]*apply(W*n_i,2,sum)
      }
      z_5 <- apply(z_4,1,sum)##

      ##part 3
      z_6 <- matrix(0,K,m)
      for (j in 1:m) {
        z_6[,j] <- ((beta[j])^2)*apply(W*n_i,2,sum)
      }
      z_7 <- apply(z_6,1,sum) ##

      z_new <- (z_3 - z_5)/z_7
      z_j <- z_new - sum(z_new*p)
      z <- z_j/sqrt(sum((z_j^2)*p))

      x_i <- data.frame()
      X <- split(data_numeric, data_factor)
      for (i in X) {
        for (j in 1:length(i[,1])) {
          xij <- colSums(data.frame(i))
        }
        x_i <- rbind(x_i,xij)
      }

      wz <- matrix(0,nr,K)
      for (k in 1:K) {
        wz[,k] <- W[,k]*z[k]
      }
      wz_2 <- matrix(0,nr,K)
      for (k in 1:K) {
        wz_2[,k] <- W[,k]*(z[k])^2
      }

      beta_1 <- apply(x_i*apply(wz,1,sum),2,sum)
      beta <- (beta_1 - (1/n)*(apply(data_numeric,2,sum)*sum(n_i*apply(wz,1,sum)))) /
        (sum(n_i*apply(wz_2,1,sum)) - (1/n)*(sum(n_i*apply(wz,1,sum)))^2)

      ##############alpha
      alpha <- (1/n)*(apply(data_numeric,2,sum) - beta*sum(n_i*apply(wz,1,sum)))

      counter = counter+1
    }

    ##sigma
    diff<- matrix(0,n,m)
    sigma <- c()
    for (k in 1:K) {
      for (l in 1:m) {
        diff[,l] <- (data_numeric[,l] - alpha[l] - beta[l]*z[k])^2
      }
      diff_i <- data.frame()
      diff_df <- data.frame(cbind(diff,data_factor))
      Y <- split(diff_df, data_factor)
      for (i in Y) {
        for (j in 1:length(i[,1])) {
          diffij <- colSums(data.frame(i))
        }
        diff_i <- rbind(diff_i,diffij)
      }
      diff_i <- as.matrix(diff_i[,-length(diff_i[1,])])###
      sigma[[k]] <- sqrt(apply(W[,k]*diff_i,2,sum)/sum(n_i*W[,k]))
    }

  }


  if(var_fun == 3){
    p <- apply(W,2,sum)/nr

    counter <- 0
    while (counter <= 5) {
      x_i <- data.frame()
      X <- split(data_numeric, data_factor)
      for (i in X) {
        for (j in 1:length(i[,1])) {
          xij <- colSums(data.frame(i))
        }
        x_i <- rbind(x_i,xij)
      }

      z_1 <- matrix(0,nr,K)
      for (i in 1:nr) {
        z_1[i,] <- W[i,]*as.numeric((t(beta)  %*% solve(matrix(unlist(sigma),m,m)) %*% matrix(matrix(unlist(x_i), nr, m)[i,] - alpha)))
      }
      z_2 <- apply(z_1,2,sum)
      z_3 <- apply(W*n_i,2,sum)*as.numeric((t(matrix(beta)) %*% solve(matrix(unlist(sigma),m,m)) %*% matrix(beta)))

      z_new <- z_2/z_3
      z_j <- z_new - sum(z_new*p)
      z <- z_j/sqrt(sum((z_j^2)*p))

      ##beta
      wz <- matrix(0,nr,K)
      for (k in 1:K) {
        wz[,k] <- W[,k]*z[k]
      }
      wz_2 <- matrix(0,nr,K)
      for (k in 1:K) {
        wz_2[,k] <- W[,k]*(z[k])^2
      }

      beta_1 <- apply(x_i*apply(wz,1,sum),2,sum)
      beta <- (beta_1 - (1/n)*(apply(data_numeric,2,sum)*sum(n_i*apply(wz,1,sum)))) /
        (sum(n_i*apply(wz_2,1,sum)) - (1/n)*(sum(n_i*apply(wz,1,sum)))^2)

      ##alpha
      alpha <- (1/n)*(apply(data_numeric,2,sum) - beta*sum(n_i*apply(wz,1,sum)))

      counter = counter+1
    }

    ##sigma
    x_part_pre <- vector("list")
    Z <- split(data_numeric, data_factor)
    for (i in Z) {
      for (k in 1:K) {
        for (j in 1:length(i[,1])) {
          x_part_pre[[j]] <- t(as.matrix(data_numeric[j,] - alpha - beta*z[k]))%*%(as.matrix(data_numeric[j,] - alpha - beta*z[k]))
          sum_x_part_pre <- Reduce("+", x_part_pre)
        }
      }
    }
    current  <- vector("list", length = nr)
    Sigma_1 <- vector("list",length = K)
    for (k in 1:K) {
      for (i in 1:nr) {
        current[[i]] <- W[i,k]*sum_x_part_pre
        sum_matrix <- Reduce("+", current)
      }
      Sigma_1[[k]] <- sum_matrix/n
      sigma <- Reduce("+", Sigma_1)
    }
  }

  if(var_fun == 4){
    p <- apply(W,2,sum)/nr

    counter <- 0
    while (counter <= 5) {
      x_i <- data.frame()
      X <- split(data_numeric, data_factor)
      for (i in X) {
        for (j in 1:length(i[,1])) {
          xij <- colSums(data.frame(i))
        }
        x_i <- rbind(x_i,xij)
      }

      z_1 <- matrix(0,nr,K)
      for (i in 1:nr) {
        for (k in 1:K) {
          z_1[,k] <- W[,k]*as.numeric((t(beta)  %*% solve(matrix(unlist(sigma[[k]]),m,m)) %*%
                                         matrix(matrix(unlist(x_i), nr, m)[i,] - alpha)))
        }
      }
      z_2 <- apply(z_1,2,sum)

      z_3 <- apply(W*n_i,2,sum)
      z_4 <- matrix(0,1,K)
      for (k in 1:K) {
        z_4[,k] <- z_3[k]*as.numeric((t(matrix(beta)) %*% solve(matrix(unlist(sigma[[k]]),m,m)) %*% matrix(beta)))
      }
      z_new <- z_2/z_4
      z_j <- z_new - sum(z_new*p)
      z <- z_j/sqrt(sum((z_j^2)*p))

      ##beta
      wz <- matrix(0,nr,K)
      for (k in 1:K) {
        wz[,k] <- W[,k]*z[k]
      }
      wz_2 <- matrix(0,nr,K)
      for (k in 1:K) {
        wz_2[,k] <- W[,k]*(z[k])^2
      }

      beta_1 <- apply(x_i*apply(wz,1,sum),2,sum)
      beta <- (beta_1 - (1/n)*(apply(data_numeric,2,sum)*sum(n_i*apply(wz,1,sum)))) /
        (sum(n_i*apply(wz_2,1,sum)) - (1/n)*(sum(n_i*apply(wz,1,sum)))^2)

      ##alpha
      alpha <- (1/n)*(apply(data_numeric,2,sum) - beta*sum(n_i*apply(wz,1,sum)))

      counter = counter+1
    }

    ##sigma
    x_part_pre <- vector("list")
    Z <- split(data_numeric, data_factor)
    for (i in Z) {
      for (k in 1:K) {
        for (j in 1:length(i[,1])) {
          x_part_pre[[j]] <- t(as.matrix(data_numeric[j,] - alpha - beta*z[k]))%*%(as.matrix(data_numeric[j,] - alpha - beta*z[k]))
          sum_x_part_pre <- Reduce("+", x_part_pre)
        }
      }
    }

    current  <- vector("list", length = nr)
    sigma <- vector("list",length = K)
    for (k in 1:K) {
      for (i in 1:nr) {
        current[[i]] <- W[i,k]*sum_x_part_pre
        add <- function(x) Reduce("+", x)
        sum_matrix <- add(current)
      }
      sigma[[k]] <- sum_matrix/sum(n_i*W[,k])
    }

  }

  return(list("p"=p, "z"=z, "alpha"=alpha, "beta"=beta, "sigma"=sigma))
}


em_2level_fun.start <- function(data, K, steps, p, z, beta, alpha, sigma, var_fun){
  pb<-txtProgressBar(min=0, max=steps, style=3)

  #data_numeric <- select_if(data, is.numeric)
  data_numeric <- data[sapply(data, is.numeric)]
  data_factor <- as.factor(data[,names(Filter(is.factor,data))])

  if(missing(var_fun)){
    var_fun <- 2
  }

  n <- length(data_numeric[,1])
  m <- length(data_numeric[1,])
  nr <- nlevels(data_factor)


  if (missing(p)){
    p <- rep(1/K, K)

  }
  if (missing(z)){
    z <- rnorm(K, mean = 0, sd=1)
  }
  if (missing(alpha)){
    alpha <- colMeans(data_numeric)
  }
  if (missing(beta)){
    beta <- as.numeric(data_numeric[sample(nrow(data_numeric), 1), ] - alpha)
  }
  if (missing(sigma)){
    if(var_fun == 1){
      s_n <- c()
      for (j in 1:m) {
        s_n[j] <- sqrt(1/(n-1)*sum((data_numeric[,j] - mean(data_numeric[,j]))^2))
      }
      sigma <- s_n/K
    }
    if(var_fun == 2){
      s_n <- c()
      sigma <- list()
      for (k in 1:K) {
        for (j in 1:m) {
          s_n[j] <- sqrt(1/(n-1)*sum((data_numeric[,j] - mean(data_numeric[,j]))^2))/K
        }
        sigma[[k]] <- s_n}
    }
    if(var_fun == 3){
      s_n <- c()
      sigma <- matrix(0,m,m)
      for (k in 1:K) {
        for (j in 1:m) {
          s_n[j] <- (1/(n-1)*sum((data[,j] - mean(data[,j]))^2))/K
        }
        mat <- matrix(0,m,m)
        diag(mat) <- s_n
        sigma <- mat
      }
    }
    if(var_fun == 4){
      s_n <- c()
      sigma <- vector("list",length = K)
      for (k in 1:K) {
        for (j in 1:m) {
          s_n[j] <- (1/(n-1)*sum((data[,j] - mean(data[,j]))^2))/K
        }
        mat <- matrix(0,m,m)
        diag(mat) <- s_n
        sigma[[k]] <- mat
      }
    }
  }

  s<-1
  while (s <=steps){
    W   <- estep_mlevel_se(data, p, z, alpha, beta, sigma, var_fun)
    fit <- mstep_mlevel_se(data, W, alpha=alpha, beta=beta, sigma = sigma,var_fun)
    p   <- fit$p
    z  <- fit$z
    beta <- fit$beta
    alpha <- fit$alpha
    sigma <-fit$sigma
    s<-s+1
    setTxtProgressBar(pb,s)
  }

  if(is.na(beta[1]) == FALSE){
    for (l in 1:length(beta)) {
      if(beta[l] == 0 & beta[l+1]<0){
        beta_hat <- beta*(-1)
        z_hat <- z*(-1)
      }
      if(beta[1] < 0){
        beta_hat <- beta*(-1)
        z_hat <- z*(-1)
      }
      if(beta[1] > 0){
        beta_hat <- beta
        z_hat <- z
      }
    }
  }

  if(is.na(beta[1]) != FALSE){
    beta_hat <- beta
    z_hat <- z
  }

  ######log-likelihood
  mu <- matrix(0,K,m)
  for (k in 1:K) {
    mu[k,] <- alpha + beta*z[k]
  }

  if(var_fun == 1){
    X <- split(data_numeric, data_factor)
    m_ik <- data.frame()
    fik_j <- matrix(0,1,K)
    for (i in X) {
      for (k in 1:K) {
        for (j in 1:length(i[,1])) {
          fik_j[,k] <- prod(dmvnorm(x=data.frame(i)[j,], mean=mu[k,], sigma=diag(sigma^2,m,m), log=FALSE))
        }
      }
      m_ik <- rbind(m_ik,fik_j)
    }
    dens.all <- matrix(0,nr,K)
    for (k in 1:K) {
      dens.all[,k] <- p[k]*m_ik[,k]
    }
  }
  if(var_fun == 2){
    X <- split(data_numeric, data_factor)
    m_ik <- data.frame()
    fik_j <- matrix(0,1,K)
    for (i in X) {
      for (k in 1:K) {
        for (j in 1:length(i[,1])) {
          fik_j[,k] <- prod(dmvnorm(x=data.frame(i)[j,], mean=mu[k,], sigma=diag(unlist(sigma[k])^2,m,m), log=FALSE))
        }
      }
      m_ik <- rbind(m_ik,fik_j)
    }
    dens.all <- matrix(0,nr,K)
    for (k in 1:K) {
      dens.all[,k] <- p[k]*m_ik[,k]
    }
  }
  if(var_fun == 3){
    X <- split(data_numeric, data_factor)
    m_ik <- data.frame()
    fik_j <- matrix(0,1,K)
    for (i in X) {
      for (k in 1:K) {
        for (j in 1:length(i[,1])) {
          fik_j[,k] <- prod(dmvnorm(x=data.frame(i)[j,], mean=mu[k,], sigma=matrix(sigma,m,m), log=FALSE))
        }
      }
      m_ik <- rbind(m_ik,fik_j)
    }
    dens.all <- matrix(0,nr,K)
    for (k in 1:K) {
      dens.all[,k] <- p[k]*m_ik[,k]
    }
  }
  if(var_fun == 4){
    X <- split(data_numeric, data_factor)
    m_ik <- data.frame()
    fik_j <- matrix(0,1,K)
    for (i in X) {
      for (k in 1:K) {
        for (j in 1:length(i[,1])) {
          fik_j[,k] <- prod(dmvnorm(x=data.frame(i)[j,], mean=mu[k,], sigma= matrix(unlist(sigma[k]),m,m), log=FALSE))
        }
      }
      m_ik <- rbind(m_ik,fik_j)
    }
    dens.all <- matrix(0,nr,K)
    for (k in 1:K) {
      dens.all[,k] <- p[k]*m_ik[,k]
    }
  }

  loglik<- sum(log(apply(dens.all,1,sum)))


  if(var_fun == 1){
    num_parameters <- (K-1) + K + m + m + m
  }
  if(var_fun == 2){
    num_parameters <- (K-1) + K + m + m + m*K
  }
  if(var_fun == 3){
    num_parameters <- (K-1) + K + m + m + m*(m+1)/2
  }
  if(var_fun == 4){
    num_parameters <- (K-1) + K + m + m + m*(m+1)/2*K
  }

  fit<- list("p"=p, "alpha"=alpha,"z"=z,"beta"=beta,"z_hat"=z_hat,"beta_hat"=beta_hat ,"sigma"=sigma,
             "W" =W, "loglikelihood"=loglik, "disparity"=-2*loglik, "number_parameters"=num_parameters,
             "AIC"=-2*loglik+2*num_parameters)

  class(fit) <- "EM"

  return(fit)
}



##################em_2level_covs.start()

estep_mlevel_covS <- function(data, v, p, z, alpha, beta, gamma, sigma, var_fun) {
  if(missing(var_fun)){
    var_fun <- 2
  }

  data <- data.frame(data)
  v <- data.frame(v)

  #data_numeric <- select_if(data, is.numeric)
  data_numeric <- data[sapply(data, is.numeric)]
  data_factor <- as.factor(data[,names(Filter(is.factor,data))])

  n <- length(data_numeric[,1])
  m <- length(data_numeric[1,])
  K <- length(p)
  r <- nlevels(data_factor)

  mu <- matrix(0,K,m)
  for (k in 1:K) {
    mu[k,] <- alpha + beta*z[k]
  }

  data_new <- data.frame()
  for (i in 1:n) {
    data_new.1 <- data_numeric[i,] - gamma%*%as.numeric(v[i,])
    data_new <- rbind(data_new,data_new.1)
  }
  if(var_fun == 1) {
    X <- split(data_new, data_factor)
    m_ik <- data.frame()
    fik_j <- matrix(0,1,K)
    for (i in X) {
      for (k in 1:K) {
        for (j in 1:length(i[,1])) {
          fik_j[,k] <- prod(dmvnorm(x=data.frame(i)[j,], mean=mu[k,], sigma=diag(sigma^2,m,m), log=FALSE))
        }
      }
      m_ik <- rbind(m_ik,fik_j)
    }
    W_1 <- matrix(0,r,K)
    for (k in 1:K) {
      W_1[,k] <- p[k]*m_ik[,k]
    }
  }

  if(var_fun == 2) {
    X <- split(data_new, data_factor)
    m_ik <- data.frame()
    fik_j <- matrix(0,1,K)
    for (i in X) {
      for (k in 1:K) {
        for (j in 1:length(i[,1])) {
          fik_j[,k] <- prod(dmvnorm(x=data.frame(i)[j,], mean=mu[k,], sigma=diag(unlist(sigma[k])^2,m,m), log=FALSE))
        }
      }
      m_ik <- rbind(m_ik,fik_j)
    }
    W_1 <- matrix(0,r,K)
    for (k in 1:K) {
      W_1[,k] <- p[k]*m_ik[,k]
    }
  }


  W <- W_1/apply(W_1,1,sum)
  return(W)

}


mstep_mlevel_covS <- function(data, v, W, alpha, beta, gamma, var_fun){
  data <- data.frame(data)
  v <- data.frame(v)
  #data_numeric <- select_if(data, is.numeric)
  data_numeric <- data[sapply(data, is.numeric)]
  data_factor <- as.factor(data[,names(Filter(is.factor,data))])
  n_i_name_numeric <- as.data.frame(summary(data_factor,maxsum = length(unique(data_factor))))
  n_i <- as.vector(n_i_name_numeric[,1])
  data_numeric_v <- cbind(data_numeric,v)

  if(missing(var_fun)){
    var_fun <- 2
  }

  n <- length(data_numeric[,1])
  m <- length(data_numeric[1,])
  nr <- nlevels(data_factor)
  K <- dim(W)[2]
  q <- as.numeric(length(v[1,]))


  p <- apply(W,2,sum)/nr

  counter <- 0
  while (counter <= 20) {
    x_i <- data.frame()
    X <- split(data_numeric, data_factor)
    for (i in X) {
      for (j in 1:length(i[,1])) {
        xij <- colSums(data.frame(i))
      }
      x_i <- rbind(x_i,xij)
    }

    v_i <- data.frame()
    V <- split(data.frame(v), data_factor)
    for (i in V) {
      for (j in 1:length(i[,1])) {
        vij <- colSums(data.frame(i))
      }
      v_i <- rbind(v_i,vij)
    }

    phi <- data.frame()
    for (i in 1:n) {
      phi_current <- t(gamma %*% t(as.matrix(v[i,])))
      phi <- rbind(phi,phi_current)
    }

    phi_i <- data.frame()
    P <- split(data.frame(phi), data_factor)
    for (i in P) {
      for (j in 1:length(i[,1])) {
        phiij <- colSums(data.frame(i))
      }
      phi_i <- rbind(phi_i,phiij)
    }


    z_1 <- matrix(0,K,m)
    for (k in 1:K) {
      z_1[k,] <- apply(W[,k]*x_i,2,sum)
    }
    z_2 <- matrix(0,K,m)
    for (j in 1:m) {
      z_2[,j]  <- beta[j]*z_1[,j]
    }
    z_3 <- apply(z_2,1,sum) ##

    z_4 <- matrix(0,K,m)
    for (j in 1:m) {
      z_4[,j] <- alpha[j]*beta[j]*apply(W*n_i,2,sum)
    }
    z_5 <- apply(z_4,1,sum)##

    ##part 3
    z_6 <- matrix(0,K,m)
    for (j in 1:m) {
      z_6[,j] <- ((beta[j])^2)*apply(W*n_i,2,sum)
    }
    z_7 <- apply(z_6,1,sum) ##

    ##part 4
    z_8 <- matrix(0,K,m)
    for (k in 1:k) {
      z_8[k,] <- apply(W[,k]*phi_i,2,sum)
    }
    z_9 <- matrix(0,K,m)
    for (j in 1:m) {
      z_9[,j]  <- beta[j]*z_8[,j]
    }
    z_10 <- apply(z_9,1,sum)

    z_new <- (z_3 - z_5 - z_10)/z_7
    z_j <- z_new - sum(z_new*p)
    z <- z_j/sqrt(sum((z_j^2)*p))

    ################beta
    wz <- matrix(0,nr,K)
    for (k in 1:K) {
      wz[,k] <- W[,k]*z[k]
    }
    wz_2 <- matrix(0,nr,K)
    for (k in 1:K) {
      wz_2[,k] <- W[,k]*(z[k])^2
    }

    beta_1 <- apply(x_i*apply(wz,1,sum),2,sum)

    beta_v_1 <- apply(phi_i*apply(wz,1,sum),2,sum)#apply(data.frame(v_i)*apply(wz,1,sum),2,sum)#
    beta_v_2 <- (1/n)*sum(n_i*apply(wz,1,sum))* apply(phi,2,sum) #sum(n_i*apply(wz,1,sum))* sum(v)#

    beta <- (beta_1 - (1/n)*(apply(data_numeric,2,sum)*sum(n_i*apply(wz,1,sum))) - beta_v_1 + beta_v_2 )/
      (sum(n_i*apply(wz_2,1,sum)) - (1/n)*(sum(n_i*apply(wz,1,sum)))^2)

    ##############alpha

    alpha <- (1/n)*(apply(data_numeric,2,sum) - beta*sum(n_i*apply(wz,1,sum)) - t(gamma %*% as.matrix(colSums(v))))

    #########gamma
    #####inverse part
    v_current <- vector("list", length = n)
    for (l in 1:n) {
      v_current[[l]] <- as.numeric(v[l,]) %*% t(as.numeric(v[l,]))
    }
    v_sum <- Reduce("+", v_current)
    v_sum_inverse <- solve(v_sum)

    ######part 1
    x_v_ij <- vector("list", length = n)
    for (i in 1:n) {
      x_v_ij[[i]] <- t(as.matrix(data_numeric[i,]))%*%as.matrix(v[i,])
    }
    sum_x_v <- Reduce("+", x_v_ij)

    ######part 2
    alpha_v_ij <- vector("list", length = n)
    for (i in 1:n) {
      alpha_v_ij[[i]] <- t(alpha)%*%as.matrix(v[i,])
    }
    sum_alpha_v <- Reduce("+", alpha_v_ij)

    ######part 3
    beta_v_ij <- vector("list", length = n)
    for (i in 1:n) {
      beta_v_ij[[i]] <- as.matrix(beta)%*%as.matrix(v[i,])
    }

    beta_v_i <-  vector("list", length = nr)
    B <- split(beta_v_ij, data_factor)
    for (b in B) {
      for (i in 1:nr) {
        beta_v_i[[i]] <-   Reduce("+", b)
      }
    }

    current_gamma <- vector("list", length = nr)
    Gamma_1 <- vector("list",length = K)
    for (k in 1:K) {
      for (i in 1:nr) {
        current_gamma[[i]] <- W[i,k]*z[k]*matrix(unlist(beta_v_i[i]),m,q)
        sum_matrix <- Reduce("+", current_gamma)
      }
      Gamma_1[[k]] <- sum_matrix
      gamma_new <- Reduce("+", Gamma_1)
    }

    gamma_7 <- sum_x_v - sum_alpha_v - gamma_new

    gamma <- gamma_7%*%v_sum_inverse




    counter = counter+1
  }

  #######sigma
  if(var_fun == 1){
    diff<-matrix(0,n, K)
    sigma <- c()
    for (l in 1:m) {
      for (k in 1:K) {
        diff[,k] <- (data_numeric[,l] - alpha[l] - beta[l]*z[k] - phi[,l])^2
      }
      diff_i <- data.frame()
      diff_df <- data.frame(cbind(diff,data_factor))
      Y <- split(diff_df, data_factor)
      for (i in Y) {
        for (j in 1:length(i[,1])) {
          diffij <- colSums(data.frame(i))
        }
        diff_i <- rbind(diff_i,diffij)
      }
      diff_i <- as.matrix(diff_i[,-length(diff_i[1,])])
      sigma[l] <- sqrt(mean(apply(W*diff_i,1,sum)))
    }
  }

  if(var_fun == 2){
    diff<- matrix(0,n,m)
    sigma <- c()
    for (k in 1:K) {
      for (l in 1:m) {
        diff[,l] <- (data_numeric[,l] - alpha[l] - beta[l]*z[k] - phi[,l])^2
      }
      diff_i <- data.frame()
      diff_df <- data.frame(cbind(diff,data_factor))
      Y <- split(diff_df, data_factor)
      for (i in Y) {
        for (j in 1:length(i[,1])) {
          diffij <- colSums(data.frame(i))
        }
        diff_i <- rbind(diff_i,diffij)
      }
      diff_i <- as.matrix(diff_i[,-length(diff_i[1,])])
      sigma[[k]] <- sqrt(apply(W[,k]*diff_i,2,sum)/sum(n_i*W[,k]))
    }
  }


  return(list("p"=p, "z"=z, "alpha"=alpha, "beta"=beta, "sigma"=sigma,"gamma"=gamma))
}


em_2level_covs.start <- function(data,v,K,steps,p,z,beta,alpha,gamma,sigma,var_fun){
  pb<-txtProgressBar(min=0, max=steps, style=3)

  data <- data.frame(data)
  v <- data.frame(v)
  data_numeric <- data[sapply(data, is.numeric)]
  data_factor <- as.factor(data[,names(Filter(is.factor,data))])

  if(missing(var_fun)){
    var_fun <- 2
  }

  n <- length(data_numeric[,1])
  m <- length(data_numeric[1,])
  nr <- nlevels(data_factor)
  q <- as.numeric(length(v[1,]))
  #K <- dim(W)[2]

  if(missing(gamma)){
    m <- as.numeric(length(data_numeric[1,]))
    q <- as.numeric(length(v[1,]))
    dep_vars <- colnames(data_numeric)
    ind_vars <- colnames(v)
    var_comb <- expand.grid(dep_vars, ind_vars )
    formula_vec <- sprintf("%s ~ %s", var_comb$Var1, var_comb$Var2)
    lm_res <- lapply( formula_vec, function(f) {
      fit1 <- lm( f, data = data.frame(data_numeric,v))
      return(fit1$coefficients[2])
    })
    gamma_new <- as.numeric(unlist(lm_res))
    gamma <- matrix(gamma_new,m,q)
  }

  if (missing(p)){
    p <- rep(1/K, K)

  }
  if (missing(z)){
    z <- rnorm(K, mean = 0, sd=1)
  }
  if (missing(alpha)){
    alpha <- colMeans(data_numeric)
  }
  if (missing(beta)){
    beta <- as.numeric(data_numeric[sample(nrow(data_numeric), 1), ] - alpha)
  }
  if (missing(sigma)){
    if(var_fun == 1){
      s_n <- c()
      for (j in 1:m) {
        s_n[j] <- sqrt(1/(n-1)*sum((data_numeric[,j] - mean(data_numeric[,j]))^2))
      }
      sigma <- s_n/K
    }
    if(var_fun == 2){
      s_n <- c()
      sigma <- list()
      for (k in 1:K) {
        for (j in 1:m) {
          s_n[j] <- sqrt(1/(n-1)*sum((data_numeric[,j] - mean(data_numeric[,j]))^2))/K
        }
        sigma[[k]] <- s_n}
    }
  }


  s<-1
  while (s <=steps){
    W   <- estep_mlevel_covS(data, v, p, z, alpha, beta, gamma, sigma, var_fun)
    fit <- mstep_mlevel_covS(data, v, W, alpha=alpha, beta=beta, gamma=gamma, var_fun)
    p   <- fit$p
    z  <- fit$z
    beta <- fit$beta
    alpha <- fit$alpha
    gamma <- fit$gamma
    sigma <-fit$sigma
    s<-s+1
    setTxtProgressBar(pb,s)
  }

  if(is.na(beta[1]) == FALSE){
    for (l in 1:length(beta)) {
      if(beta[l] == 0 & beta[l+1]<0){
        beta_hat <- beta*(-1)
        z_hat <- z*(-1)
      }
      if(beta[1] < 0){
        beta_hat <- beta*(-1)
        z_hat <- z*(-1)
      }
      if(beta[1] > 0){
        beta_hat <- beta
        z_hat <- z
      }
    }
  }

  if(is.na(beta[1]) != FALSE){
    beta_hat <- beta
    z_hat <- z
  }


  mu <- matrix(0,K,m)
  for (k in 1:K) {
    mu[k,] <- alpha + beta*z[k]
  }

  data_new <- data.frame()
  for (i in 1:n) {
    data_new.1 <- data_numeric[i,] -  gamma%*%as.numeric(v[i,])#gamma*v[i]
    data_new <- rbind(data_new,data_new.1)
  }
  ######log-likelihood
  if(var_fun == 1){
    X <- split(data_new, data_factor)
    m_ik <- data.frame()
    fik_j <- matrix(0,1,K)
    for (i in X) {
      for (k in 1:K) {
        for (j in 1:length(i[,1])) {
          fik_j[,k] <- prod(dmvnorm(x=data.frame(i)[j,], mean=mu[k,], sigma=diag(sigma^2,m,m), log=FALSE))
        }
      }
      m_ik <- rbind(m_ik,fik_j)
    }
    dens.all <- matrix(0,nr,K)
    for (k in 1:K) {
      dens.all[,k] <- p[k]*m_ik[,k]
    }
  }
  if(var_fun == 2){
    X <- split(data_new, data_factor)
    m_ik <- data.frame()
    fik_j <- matrix(0,1,K)
    for (i in X) {
      for (k in 1:K) {
        for (j in 1:length(i[,1])) {
          fik_j[,k] <- prod(dmvnorm(x=data.frame(i)[j,], mean=mu[k,], sigma=diag(unlist(sigma[k])^2,m,m), log=FALSE))
        }
      }
      m_ik <- rbind(m_ik,fik_j)
    }
    dens.all <- matrix(0,nr,K)
    for (k in 1:K) {
      dens.all[,k] <- p[k]*m_ik[,k]
    }
  }



  loglik<- sum(log(apply(dens.all,1,sum)))

  m <- as.numeric(m)
  q <- as.numeric(q)

  if(var_fun == 1){
    num_parameters <- (K-1) + K + m + m + m + m*q
  }
  if(var_fun == 2){
    num_parameters <- (K-1) + K + m + m + m*K + m*q
  }


  fit<- list("p"=p, "alpha"=alpha,"z"=z,"beta"=beta,"z_hat"=z_hat,"beta_hat"=beta_hat ,"sigma"=sigma,"gamma"=gamma,
             "W" =W, "loglikelihood"=loglik, "disparity"=-2*loglik, "number_parameters"=num_parameters,
             "AIC"=-2*loglik+2*num_parameters)

  class(fit) <- "EM"

  return(fit)

}

