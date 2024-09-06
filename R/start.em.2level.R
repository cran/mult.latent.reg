#' @title EM algorithm for multivariate two level model with covariates
#' @name mult.em_2level
#' @description This function extends the one-level version \link{mult.em_1level},
#' and it is designed to obtain Maximum Likelihood Estimates (MLE) using the EM algorithm for nested (structured) multivariate data,
#' e.g. multivariate test scores (such as on numeracy, literacy)  of students nested in different classes or schools.
#' The resulting estimates can be applied for clustering or constructing league tables (ranking of observations).
#' With the inclusion of covariates, the model allows fitting a multivariate response model for further regression analysis.
#' Detailed information about the model used in this function can be found in Zhang et al. (2023).
#' @param data A data set object; we denote the dimension to be \eqn{m}.
#' @param v Covariate(s).
#' @param K Number of mixture components, the default is \code{K = 2}. Note that when \code{K = 1}, \code{z} and \code{beta} will be 0.
#' @param steps Number of iterations, the default is \code{steps = 20}.
#' @param start Containing parameters involved in the proposed model (\code{p}, \code{alpha}, \code{z}, \code{beta}, \code{sigma}, \code{gamma}) in a list,
#'              the starting values can be obtained through the use of \link{start_em}. More details can be found in \link{start_em}.
#' @param option Four options for selecting the starting values for the parameters in the model. The default is \code{option = 1}.
#'               More details can be found in \link{start_em}.
#' @param var_fun There are two types of variance specifications; \code{var_fun = 1}, the same diagonal variance specification to all \code{K} components of the mixture;
#'                \code{var_fun = 2}, different diagonal variance matrices for different components;
#'                The default is \code{var_fun = 2}.
#' @references Zhang, Y., Einbeck, J. and Drikvandi, R. (2023). A multilevel multivariate response model for data with latent structures.
#'             In: Proceedings of the 37th International Workshop on Statistical Modelling, pages 343-348.
#'             Link on RG: \url{https://www.researchgate.net/publication/375641972_A_multilevel_multivariate_response_model_for_data_with_latent_structures}
#' @return The estimated parameters in the model \eqn{x_{ij} = \alpha + \beta z_k + \Gamma v_{ij} + \varepsilon_{ij}} obtained through the EM algorithm,
#'         where the upper-level unit is indexed by \eqn{i}, and the lower-level unit is indexed by \eqn{j}.
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
#'  \item{starting_values}{A list of starting values for parameters used in the EM algorithm.}
#' @seealso \code{\link{mult.reg_2level}}.
#'
#' @examples
#' \donttest{
#' ##examples for data without covariates.
#' data(trading_data)
#' set.seed(49)
#' trade_res <- mult.em_2level(trading_data, K=4, steps = 10, var_fun = 2)
#'
#' i_1 <- apply(trade_res$W, 1, which.max)
#'ind_certain <- rep(as.vector(i_1),c(4,5,5,3,5,5,4,4,5,5,5,5,5,5,5,5,5,5,
#'3,5,5,5,5,4,4,5,5,5,4,5,4,5,5,5,3,5,5,5,5,5,5,4,5,4))
#'colors <- c("#FF6600","#66BD63", "lightpink","purple")
#'plot(trading_data[,-3],pch = 1, col = colors[factor(ind_certain)])
#'legend("topleft", legend=c("Mass point 1", "Mass point 2","Mass point 3","Mass point 4"),
#' col=c("#FF6600","purple","#66BD63","lightpink"),pch = 1, cex=0.8)
#'
#'###The Twins data
#'library(lme4)
#'set.seed(26)
#'twins_res <- mult.em_2level(twins_data[,c(1,2,3)],v=twins_data[,c(4,5,6)],
#'K=2, steps = 20, var_fun = 2)
#'coeffs <- twins_res$gamma
#' ##Compare to the estimated coefficients obtained using individual two-level models (lmer()).
#'summary(lmer(SelfTouchCodable ~ Depression + PSS + Anxiety + (1 | id) ,
#'data=twins_data, REML = TRUE))$coefficients[2,1]
#'}
#' @import "mvtnorm"
#' @import "stats"
#' @import "utils"
#' @import "matrixStats"
#' @import "lme4"
#' @export
library(matrixStats)
library(mvtnorm)
mult.em_2level <- function(data, v, K, start, steps = 10, var_fun = 2, option = 1){


  if (missing(v)) {

    result <- em_2level_fun(data, K, start, steps,  var_fun, option)

  } else {

    result <- em_2level_covs(data, v, K, start, steps, var_fun, option)

  }
  return(result)
}




##########2level functions
estep_mlevel_se <- function(data, p, z, alpha, beta, sigma, var_fun){

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


em_2level_fun <- function(data, K, start, steps, var_fun, option){
  pb<-txtProgressBar(min=0, max=steps, style=3)
  data_numeric <- data[sapply(data, is.numeric)]
  data_factor <- as.factor(data[,names(Filter(is.factor,data))])

  if(missing(var_fun)){
    var_fun <- 2
  }
  if(missing(option)){
    option <- 1
  }
  n <- length(data_numeric[,1])
  m <- length(data_numeric[1,])
  nr <- nlevels(data_factor)

  if(missing(start)){
    start <- start_em(data = data, K = K, var_fun = var_fun, steps = steps, option = option)
    p <- start$p
    alpha <- start$alpha
    beta <- start$beta
    z <- start$z
    sigma <- start$sigma
  } else {
    p <- start$p
    alpha <- start$alpha
    beta <- start$beta
    z <- start$z
    sigma <- start$sigma
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

  fit<- list("p"=p, "alpha"=alpha,"z"=z_hat,"beta"=beta_hat ,"sigma"=sigma,
             "W" =W, "loglikelihood"=loglik, "disparity"=-2*loglik, "number_parameters"=num_parameters,
             "AIC"=-2*loglik+2*num_parameters, "starting_values" = start)

  class(fit) <- "EM"

  return(fit)
}


#####two-level+covs function
estep_mlevel_covS <- function(data, v, p, z, alpha, beta, gamma, sigma, var_fun) {
  if(missing(var_fun)){
    var_fun <- 2
  }

  data <- data.frame(data)
  v <- as.data.frame(v)

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
  v <- as.data.frame(v)

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
  q <- length(v[1,])


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

    B <- split(beta_v_ij, data_factor)
    beta_v_i <- vector("list", length = length(B))
    for (i in seq_along(B)) {
      matrices_list <- B[[i]]
      beta_v_i[[i]] <- Reduce("+", matrices_list)
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

    phi <- data.frame()
    for (i in 1:n) {
      phi_current <- t(gamma %*% t(as.matrix(v[i,])))
      phi <- rbind(phi,phi_current)
    }

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

    phi <- data.frame()
    for (i in 1:n) {
      phi_current <- t(gamma %*% t(as.matrix(v[i,])))
      phi <- rbind(phi,phi_current)
    }


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




em_2level_covs <- function(data, v, K, start, steps, var_fun, option){
  pb<-txtProgressBar(min=0, max=steps, style=3)
  data <- data.frame(data)
  v <- as.data.frame(v)
  data_numeric <- data[sapply(data, is.numeric)]
  data_factor <- as.factor(data[,names(Filter(is.factor,data))])

  if(missing(var_fun)){
    var_fun <- 2
  }
  if(missing(option)){
    option <- 1
  }
  n <- length(data_numeric[,1])
  m <- length(data_numeric[1,])
  nr <- nlevels(data_factor)
  q <- as.numeric(length(v[1,]))

  if(missing(start)){
    start <- start_em(data = data, v = v, K = K, var_fun = var_fun, steps = steps, option = option)
    p <- start$p
    alpha <- start$alpha
    beta <- start$beta
    z <- start$z
    sigma <- start$sigma
    gamma <- start$gamma
  } else {
    p <- start$p
    alpha <- start$alpha
    beta <- start$beta
    z <- start$z
    sigma <- start$sigma
    gamma <- start$gamma
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

  if(var_fun == 1){
    num_parameters <- (K-1) + K + m + m + m + m*q
  }
  if(var_fun == 2){
    num_parameters <- (K-1) + K + m + m + m*K + m*q
  }


  fit<- list("p"=p, "alpha"=alpha,"z"=z_hat,"beta"=beta_hat ,"sigma"=sigma,"gamma"=gamma,
             "W" =W, "loglikelihood"=loglik, "disparity"=-2*loglik, "number_parameters"=num_parameters,
             "AIC"=-2*loglik+2*num_parameters, "starting_values" = start)

  class(fit) <- "EM"

  return(fit)

}

