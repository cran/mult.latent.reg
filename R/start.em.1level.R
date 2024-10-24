#' @title EM algorithm for multivariate one level model with covariates
#' @name mult.em_1level
#' @description This function is used to obtain the Maximum Likelihood Estimates (MLE) using the EM algorithm for one-level multivariate data.
#'              The estimates enable users to conduct clustering, ranking, and simultaneous dimension reduction on the multivariate dataset.
#'              Furthermore, when covariates are included, the function supports the fitting of multivariate response models, expanding its utility for regression analysis.
#'              The details of the model used in this function can be found in Zhang and Einbeck (2024). Note that this function is designed for multivariate data.
#'              When the dimension of the data is 1, please use \link[npmlreg]{alldist} as an alternative. A warning message will also be displayed when the input data is a univariate dataset.
#' @param data A data set object; we denote the dimension to be \eqn{m}.
#' @param v Covariate(s).
#' @param K Number of mixture components, the default is \code{K = 2}. Note that when \code{K = 1}, \code{z} and \code{beta} will be 0.
#' @param steps Number of iterations, the default is \code{steps = 20}.
#' @param start Containing parameters involved in the proposed model (\code{p}, \code{alpha}, \code{z}, \code{beta}, \code{sigma}, \code{gamma}) in a list,
#'              the starting values can be obtained through the use of \link{start_em}. More details can be found in \link{start_em}.
#' @param option Four options for selecting the starting values for the parameters in the model. The default is option = 1.
#'                More details can be found in \link{start_em}.
#' @param var_fun There are four types of variance specifications;
#'                \code{var_fun = 1}, the same diagonal variance specification to all \code{K} components of the mixture;
#'                \code{var_fun = 2}, different diagonal variance matrices for different components.
#'                \code{var_fun = 3}, the same full (unrestricted) variance for all components.
#'                \code{var_fun = 4}, different full (unrestricted) variance matrices for different components.
#'                The default is \code{var_fun = 2}.
#' @return The estimated parameters in the model \eqn{x_{i} = \alpha + \beta z_k + \Gamma v_i + \varepsilon_i} obtained through the EM algorithm at the convergence.
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
#'  \item{starting_values}{A list of starting values for parameters used in the EM algorithm.}
#' @seealso \code{\link{mult.reg_1level}}.
#' @note It is worth noting that due to the sequential nature of the updates within the M-step,
#'      this algorithm can be considered an ECM algorithm.
#' @references Zhang, Y. and Einbeck, J. (2024). A Versatile Model for Clustered and Highly Correlated Multivariate Data. J Stat Theory Pract 18(5).\doi{10.1007/s42519-023-00357-0}
#' @examples
#' ##example for data without covariates.
#' data(faithful)
#' res <- mult.em_1level(faithful,K=2,steps = 10,var_fun = 1)
#'
#'
#' ## Graph showing the estimated one-dimensional space with cluster centers in red and alpha in green.
#' x <- res$alpha[1]+res$beta[1]*res$z
#' y <- res$alpha[2]+res$beta[2]*res$z
#' plot(faithful,col = 8)
#' points(x=x[1],y=y[1],type = "p",col = "red",pch = 17)
#' points(x=x[2],y=y[2],type = "p",col = "red",pch = 17)
#' points(x=res$alpha[1],y=res$alpha[2],type = "p",col = "darkgreen",pch = 4)
#' slope <- (y[2]-y[1])/(x[2]-x[1])
#' intercept <- y[1]-slope*x[1]
#' abline(intercept, slope, col="red")
#'
#' ##Graph showing the originaldata points being assigned to different
#'  ##clusters according to the Maximum a posterior (MAP) rule.
#' index <- apply(res$W, 1, which.max)
#' faithful_grouped <- cbind(faithful,index)
#' colors <- c("#FDAE61", "#66BD63")
#' plot(faithful_grouped[,-3], pch = 1, col = colors[factor(index)])
#'
#' \donttest{
#' ##example for data with covariates.
#' data(fetal_covid_data)
#' set.seed(2)
#' covid_res <- mult.em_1level(fetal_covid_data[,c(1:5)],v=fetal_covid_data$status_bi, K=3, steps = 20,
#'              var_fun = 2)
#' coeffs <- covid_res$gamma
#' ##compare with regression coefficients from fitting individual linear models.
#' summary(lm( UpperFaceMovements ~ status_bi,data=fetal_covid_data))$coefficients[2,1]
#' summary(lm( Headmovements ~ status_bi,data=fetal_covid_data))$coefficients[2,1]
#'}
#' @import "mvtnorm"
#' @import "stats"
#' @import "utils"
#' @import "matrixStats"
#' @export
library(matrixStats)
library(mvtnorm)
library(stats)
library(utils)
mult.em_1level <- function(data, v, K, start, steps = 10, var_fun = 2, option = 1) {
  #data_numeric <- data[sapply(data, is.numeric)]
  #m <- length(data_numeric[1, ])
  data <- as.data.frame(data)
  m <- length(data[1, ])

  if (m == 1 ) {
    warning("Please use alldist() available in npmlreg in the case of m = 1")
    return()
  }


  if (missing(v)) {

    result <- em_fun(data, K, start, steps, var_fun, option)
  } else {

    result <- em_covs(data, v, K, start, steps, var_fun, option)

  }
  return(result)


}


###em_fun
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
    n <- length(data[,1]) #n is added
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
      while (counter <= 5) {
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
      while (counter <= 5) {
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
      while (counter <= 5) {
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
      while (counter <= 5) {
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


##################################

em_fun <- function(data, K, start, steps, var_fun, option){
  pb<-txtProgressBar(min=0, max=steps, style=3)
  data <- data.frame(data)
  m <- length(data[1,])
  n <- length(data[,1])

  if(missing(var_fun)){
    var_fun <- 2
  }
  if(missing(option)){
    option <- 1
  }
  #if(m != 1 && K == 1){
    if(missing(start)){
      start <- start_em(data = data, K = K ,var_fun = var_fun, steps = steps,option = option)
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
  #}


 # if(m != 1 && K == 1){
  #  if (missing(p)){
   #   p <- rep(1/K, K)
  #  }
   # if (missing(z)){
    #  z <- rnorm(K,mean = 0, sd=1)
    #}
    #if (missing(alpha)){
     # alpha <- colMeans(data)
    #}
    #if (missing(beta)){
     # beta <- as.numeric(data[sample(nrow(data), 1), ] - alpha)
    #}
    #if (missing(sigma)){
     # if(var_fun == 1){
      #  s_n <- c()
       # for (j in 1:m) {
        #  s_n[j] <- sqrt(1/(n-1)*sum((data[,j] - mean(data[,j]))^2))
        #}
        #sigma <- s_n/K
     # }
      #if(var_fun == 2){
       # s_n <- c()
        #sigma <- list()
        #for (k in 1:K) {
         # for (j in 1:m) {
          #  s_n[j] <- sqrt(1/(n-1)*sum((data[,j] - mean(data[,j]))^2))/K
          #}
          #sigma[[k]] <- s_n}
      #}
      #if(var_fun == 3){
       # v3<-em.diag(data,K,steps=steps+20,var_fun = 3)
      #  sigma <- v3$sigma
      #}
      #if(var_fun == 4){
       # v3<-em.diag(data,K,steps=steps+20,var_fun = 4)
        #sigma <- v3$sigma
      #}
   # }
  #}

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
  ##"projection" = x_proj
  fit<- list("p"=p, "alpha"=alpha,"z"=z_hat,"beta"=beta_hat ,"sigma"=sigma,
             "W" =W, "loglikelihood"=loglik, "disparity"=-2*loglik, "number_parameters"=num_parameters,
             "AIC"=-2*loglik+2*num_parameters, "BIC"= -2*loglik + num_parameters*log(length(data[,1])) , "starting_values" = start)
  class(fit) <- "EM"
  return(fit)
}


###############em_covs
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


mstep_covs <- function(data, v, W, alpha, beta, gamma, var_fun){
  if(missing(var_fun)){
    var_fun <- 2
  }
  m <- length(data[1,])
  K <- dim(W)[2] ##newly added

  if(m != 1 && K != 1){
    data <- data.frame(data)
    v <- as.data.frame(v)
    m <- length(data[1,])
    n <- dim(W)[1]
    K <- dim(W)[2]

    if(var_fun == 1){
      p <- apply(W,2,sum)/n

      counter <- 0
      while (counter <= 5) {
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

      phi <- data.frame() ##newly added
      for (i in 1:n) {
        phi_current <- t(gamma %*% t(as.matrix(v[i,])))
        phi <- rbind(phi,phi_current)
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
      while (counter <= 5) {
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

      phi <- data.frame()
      for (i in 1:n) {
        phi_current <- t(gamma %*% t(as.matrix(v[i,])))
        phi <- rbind(phi,phi_current)
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
      while (counter <= 5) {
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

      phi <- data.frame()
      for (i in 1:n) {
        phi_current <- t(gamma %*% t(as.matrix(v[i,])))
        phi <- rbind(phi,phi_current)
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
      while (counter <= 5) {
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
        ###alpha
        wz <- matrix(0,n,K)
        for (k in 1:K) {
          wz[,k] <- W[,k]*z[k]
        }
        alpha <- (1/n)*(apply(data,2,sum) - beta*sum(apply(wz,1,sum))  - apply(phi,2,sum))

        ####beta
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

        ###gamma
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

      phi <- data.frame()
      for (i in 1:n) {
        phi_current <- t(gamma %*% t(as.matrix(v[i,])))
        phi <- rbind(phi,phi_current)
      }

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


  if(K == 1){
    n <- as.numeric(length(data[,1]))
    W <- matrix(1,n,1)
    m <- as.numeric(length(data[1,]))
    data <- data.frame(data)
    v <- as.data.frame(v)
    K <- dim(W)[2]

    p <- 1
    z <- 0
    beta <- rep(0, m)
    if(var_fun == 1){
      p <- 1
      beta <- rep(0, m)
      z <- 0
      counter <- 0
      while (counter <= 5) {
        ####alpha
        phi <- data.frame() ##newly added
        for (i in 1:n) {
          phi_current <- t(gamma %*% t(as.matrix(v[i,])))
          phi <- rbind(phi,phi_current)
        }
        alpha <- (1/n)*(apply(data,2,sum)  - apply(phi,2,sum))

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

        counter = counter + 1
      }
      ####sigma
      phi <- data.frame() ##newly added
      for (i in 1:n) {
        phi_current <- t(gamma %*% t(as.matrix(v[i,])))
        phi <- rbind(phi,phi_current)
      }

      diff<-matrix(0,n, K)
      sigma <- c()
      for (j in 1:m) {
        for (k in 1:K) {
          diff[,k] <- (data[,j] - alpha[j] - phi[,j])^2
        }
        sigma[j] <- sqrt(mean(apply(W*diff,1,sum)))
      }
    }

    if (var_fun == 3){
      p <- 1
      beta <- rep(0, m)
      z <- 0
      counter <- 0
      while (counter <= 5) {
        ####alpha
        phi <- data.frame() ##newly added
        for (i in 1:n) {
          phi_current <- t(gamma %*% t(as.matrix(v[i,])))
          phi <- rbind(phi,phi_current)
        }
        alpha <- (1/n)*(apply(data,2,sum)  - apply(phi,2,sum))

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

        counter = counter + 1
      }
      ####sigma
      phi <- data.frame()
      for (i in 1:n) {
        phi_current <- t(gamma %*% t(as.matrix(v[i,])))
        phi <- rbind(phi,phi_current)
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

  }


  return(list("p"=p, "z"=z, "alpha"=alpha, "beta"=beta, "gamma"=gamma,"sigma"=sigma))
}



em_covs <-  function(data, v, K, start, steps, var_fun, option){
  pb<-txtProgressBar(min=0, max=steps, style=3)

  data <- data.frame(data)
  v <- as.data.frame(v)
  q <- as.numeric(length(v[1,]))

  if(missing(var_fun)){
    var_fun <- 2
  }
  if(missing(option)){
    option <- 1
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
    if(missing(start)){
      start <- start_em(data = data, v = v, K = K , var_fun = var_fun , steps = steps, option = option)
      p <- start$p
      alpha <- start$alpha
      beta <- start$beta
      z <- start$z
      gamma <- start$gamma
      sigma <- start$sigma
    } else {
      p <- start$p
      alpha <- start$alpha
      beta <- start$beta
      z <- start$z
      gamma <- start$gamma
      sigma <- start$sigma
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
    if(K==1){
      W <- matrix(1,n,1)
      p <- 1
      z <- 0
      beta <- rep(0, m)
      fit <- mstep_covs(data, v, W, alpha=alpha, beta=beta, gamma=gamma,var_fun)
      alpha <- fit$alpha
      gamma <- fit$gamma
      sigma <-fit$sigma

    } else {

    W   <- estep_covs(data,v,p,z,alpha,beta,gamma,sigma,var_fun)
    fit <- mstep_covs(data, v, W, alpha=alpha, beta=beta, gamma=gamma,var_fun)
    p   <- fit$p
    z  <- fit$z
    beta <- fit$beta
    alpha <- fit$alpha
    gamma <- fit$gamma
    sigma <-fit$sigma
    }

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

  if(K == 1){
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

  fit<- list("p"=p, "alpha"=alpha, "z"=z_hat,"beta"=beta_hat, "gamma"=gamma,"sigma"=sigma,
             "W" =W, "loglikelihood"=loglik, "disparity"=-2*loglik, "number_parameters"=num_parameters,
             "AIC"=-2*loglik+2*num_parameters,"BIC"= -2*loglik + num_parameters*log(length(data[,1])),
             "starting_values" = start)

  class(fit) <- "EM"
  return(fit)
}



