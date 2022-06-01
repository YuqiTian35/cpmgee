#' Calculate CDFs conditional on covariates
#'
#' This functions calculates the CDFs conditional on covariates based on the fitted model
#' and new data.
#'
#' @param mod the model
#' @param cont_y the continuous response data
#' @param new.data the new data
#' @param at.y a numeric vector of cut-off points P(y <= at.y | new.data)
#' @param se if confidence intervals needed (default = TRUE)
#' @return A list containing the following components:
#' @return \item{est}{a vector of estimated condtional CDFs}
#' @return \item{se}{a vector of estimated standard errors}
#' @return \item{lb}{a vector of estimated lower bounds of 95\% confidence intervals}
#' @return \item{ub}{a vector of estimated upper bounds of 95\% confidence intervals}
#' @export
cdf_cpmgee <- function(mod, cont_y, new.data, at.y=0, se=TRUE){
  # inputs
  coef <- mod$coefficients
  cov.par <- mod$robust.var
  order.y <- sort(unique(cont_y))
  alpha <- coef[1:(length(order.y)-1)]
  beta <- coef[length(order.y):length(coef)]
  cumprob <- function(x) 1 / (1 + exp(-x))
  
  # x %*% b
  xb <- as.matrix(new.data) %*% beta
  
  index <- sapply(at.y, FUN=function(x) {
    if(x<min(order.y)[1]) result <- Inf
    else if (x==min(order.y)[1]) result <- 1
    else if(x >= max(order.y)[1]) result <- -Inf
    else which(order.y>=x)[1]-1})
  
  m.alpha <-  coef[index]
  m.alpha <- ifelse(is.infinite(index), index, m.alpha)
  if(length(at.y)==1){
    lb <- as.matrix(outer(m.alpha, xb, "+")[,,1])
  } else lb <- t(outer(m.alpha, xb, "+")[,,1])
  m.cdf <- cumprob(lb)
  
  
  if(se){
    
    deriv <- function(x) exp(-x)/(1+exp(-x))^2
    
    cdf.se <- matrix(NA, ncol=length(at.y), nrow=dim(new.data)[1])
    lb.se <- matrix(NA, ncol=length(at.y), nrow=dim(new.data)[1])
    
    var <- as.matrix(cov.par)
    
    for(i in 1:length(at.y)) {
      
      var.i <- var[c(index[i], length(order.y):length(coef)),
                   c(index[i], length(order.y):length(coef))]
      
      dcdf.dtheta <- cbind(-deriv(lb[,i]),
                           -deriv(lb[,i]) * as.matrix(new.data) )
      dlb.dtheta <- as.matrix(cbind(1, new.data))
      cdf.se[,i] <- sqrt(diag(dcdf.dtheta %*% var.i%*% t(dcdf.dtheta)))
      lb.se[, i] <- sqrt(diag(dlb.dtheta%*%var.i%*% t(dlb.dtheta)))
      
    }
    ci.lb <- sapply(1:length(at.y), FUN=function(i)
    { cumprob(lb[, i]  - qnorm(0.975)*lb.se[, i])})
    ci.ub <- sapply(1:length(at.y), FUN=function(i)
    { cumprob(lb[, i]  + qnorm(0.975)*lb.se[, i])})
    
    
    result <- list(est=m.cdf,
                   se=cdf.se,
                   lb=ci.lb,
                   ub=ci.ub)
  } else{
    result <- list(est = m.cdf)
  }
  
  return(result)
}


#' Calculate quantiles conditional on covariates
#'
#' This functions calculates the quantiles conditional on covariates based on the fitted model
#' and new data.
#'
#' @param mod the model
#' @param cont_y the continuous response data
#' @param new.data the new data
#' @param probs a numeric vector of pth quantiles
#' @param se if confidence intervals needed (default = TRUE)
#' @return A list containing the following components:
#' @return \item{est}{a vector of estimated condtional quantiles}
#' @return \item{lb}{a vector of estimated lower bounds of 95\% confidence intervals}
#' @return \item{ub}{a vector of estimated upper bounds of 95\% confidence intervals}
#' @export
quantile_cpmgee <- function(mod, cont_y, new.data, probs=0.5, se=TRUE){
  # inputs
  coef <- mod$coefficients
  cov.par <- mod$robust.var
  order.y <- sort(unique(cont_y))
  alpha <- coef[1:(length(order.y)-1)]
  beta <- coef[length(order.y):length(coef)]
  cumprob <- function(x) 1 / (1 + exp(-x))
  
  # type 2
  order.y_left <- c(order.y, order.y[length(order.y)])
  # type 1
  order.y_right <- c(order.y[1], order.y)
  
  quantile <- matrix(NA, nrow=dim(new.data)[1], ncol=length(probs))
  xb <- as.matrix(new.data) %*% beta
  m.alpha <- alpha
  lb <- t(outer(m.alpha, xb, "+")[,,1])
  
  m.cdf <- cumprob(lb)
  m.cdf <- cbind(0, m.cdf, 1)
  
  for(i in 1: length(probs)){
    try({
      index.1 <- apply(m.cdf, 1, FUN=function(x){ max(which(x<=probs[i]))[1]} )
      index.2 <- apply(m.cdf, 1, FUN=function(x){ min(which(x>=probs[i]))[1]} )
      
      # index - type 2
      index.y1_left <- ifelse(index.1 > length(order.y_left), Inf, order.y_left[index.1])
      index.y2_left <- ifelse(index.2 > length(order.y_left), Inf, order.y_left[index.2])
      # index - type 1
      index.y1_right <- ifelse(index.1 > length(order.y_right), Inf, order.y_right[index.1])
      index.y2_right <- ifelse(index.2 > length(order.y_right), Inf, order.y_right[index.2])
      
      # CDF - type 2
      index.y1.cdf_left <- ifelse(index.1 == 0, 0, m.cdf[cbind(1:dim(new.data)[1], index.1)])
      index.y2.cdf_left <- ifelse(index.2 > length(order.y_left), 1, m.cdf[cbind(1:dim(new.data)[1], index.2)])
      # CDF - type 1
      index.y1.cdf_right <- ifelse(index.1 == 0, 0, m.cdf[cbind(1:dim(new.data)[1], index.1)])
      index.y2.cdf_right <- ifelse(index.2 > length(order.y_right), 1, m.cdf[cbind(1:dim(new.data)[1], index.2)])
      
      # quantile - type 2
      quantile_left <- ifelse(index.1 == index.2, index.y1_left,
                              (index.y2_left - index.y1_left) / (index.y2.cdf_left - index.y1.cdf_left) *
                                (probs[i] - index.y1.cdf_left) + index.y1_left)
      quantile_left <- ifelse(is.na(quantile_left), max(order.y), quantile_left)
      # quantile - type 1
      quantile_right <- ifelse(index.1 == index.2, index.y1_right,
                               (index.y2_right - index.y1_right) / (index.y2.cdf_right - index.y1.cdf_right) *
                                 (probs[i] - index.y1.cdf_right) + index.y1_right)
      quantile_right <- ifelse(is.na(quantile_right), max(order.y), quantile_right)
      
      # weight
      all_weights <- c(0, (m.cdf[-length(m.cdf)] - m.cdf[1]) / (m.cdf[length(m.cdf) - 1] - m.cdf[1]), 1)
      weight <- all_weights[index.1]
      # weighted average
      quantile[,i] <- weight * quantile_left + (1 - weight) * quantile_right
    })
  }
  
  result <- quantile
  
  if(se){
    deriv <- function(x) exp(-x)/(1+exp(-x))^2
    
    quantile.lb <- quantile.ub <- matrix(NA, nrow=dim(new.data)[1], ncol=length(probs))
    lb.se <- matrix(NA, ncol=dim(lb)[2], nrow=dim(new.data)[1])
    var <- as.matrix(cov.par)
    
    for(i in 1:dim(lb)[2]){
      var.i <- var[c(i, length(order.y):length(coef)),
                   c(i, length(order.y):length(coef))]
      
      dcdf.dtheta <- cbind(deriv(lb[,i]),
                           deriv(lb[,i])*as.matrix(new.data) )
      dlb.dtheta <- as.matrix(cbind(1, new.data))
      lb.se[,i] <- sqrt(diag(dlb.dtheta %*% var.i %*% t(dlb.dtheta)))
    }
    
    ci.lb <- sapply(1:dim(lb)[2], FUN=function(i) { cumprob(lb[, i] - qnorm(0.975) * lb.se[, i])})
    ci.ub <- sapply(1:dim(lb)[2], FUN=function(i) { cumprob(lb[, i] + qnorm(0.975) * lb.se[, i])})
    ci.lb <- matrix(ci.lb, nrow=dim(new.data)[1])
    ci.ub <- matrix(ci.ub, nrow=dim(new.data)[1])
    
    ci.lb <- cbind(0, ci.lb, 1)
    ci.ub <- cbind(0, ci.ub, 1)
    
    for(i in 1: length(probs)){
      try({
        # lower bound
        index.1 <- apply(ci.lb, 1, FUN=function(x){ max(which( x<= probs[i]))[1]} )
        index.2 <- apply(ci.lb, 1, FUN=function(x){ min(which(x >= probs[i]))[1]} )
        
        # index - type 2
        index.y1_left <- ifelse(index.1>length(order.y_left), Inf, order.y_left[index.1])
        index.y2_left <- ifelse(index.2>length(order.y_left), Inf, order.y_left[index.2])
        # index - type 1
        index.y1_right <- ifelse(index.1>length(order.y_right), Inf, order.y_right[index.1])
        index.y2_right <- ifelse(index.2>length(order.y_right), Inf, order.y_right[index.2])
        
        # CDF - type 2
        index.y1.cdf_left <- ifelse(index.1 == 0, 0, ci.lb[cbind(1:dim(new.data)[1], index.1)])
        index.y2.cdf_left <- ifelse(index.2 > length(order.y_left), 1, ci.lb[cbind(1:dim(new.data)[1], index.2)])
        # CDF - type 1
        index.y1.cdf_right <- ifelse(index.1 == 0, 0, ci.lb[cbind(1:dim(new.data)[1], index.1)])
        index.y2.cdf_right <- ifelse(index.2 > length(order.y_right), 1, ci.lb[cbind(1:dim(new.data)[1], index.2)])
        
        # quantile - type 2
        quantile_left <- ifelse(index.1==index.2, index.y1_left,
                                (index.y2_left - index.y1_left) / (index.y2.cdf_left - index.y1.cdf_left) *
                                  (probs[i] - index.y1.cdf_left) + index.y1_left)
        quantile_left <- ifelse(is.infinite(quantile_left), max(order.y), quantile_left)
        # quantile - type 1
        quantile_right <- ifelse(index.1==index.2, index.y1_right,
                                 (index.y2_right - index.y1_right)/(index.y2.cdf_right - index.y1.cdf_right) *
                                   (probs[i] - index.y1.cdf_right) + index.y1_right)
        quantile_right <- ifelse(is.infinite(quantile_right), max(order.y), quantile_right)
        
        # weight
        all_weights <- c(0, (ci.lb[-length(ci.lb)] - ci.lb[1]) / (ci.lb[length(ci.lb) - 1] - ci.lb[1]), 1)
        weight <- all_weights[index.1]
        # weighted average
        quantile.lb[,i] <- weight * quantile_left + (1 - weight) * quantile_right
        
        
        # upper bound
        index.1 <- apply(ci.ub, 1, FUN=function(x){ max(which(x<=probs[i]))[1]} )
        index.2 <- apply(ci.ub, 1, FUN=function(x){ min(which(x>=probs[i]))[1]} )
        
        # index - type 2
        index.y1_left <- ifelse(index.1>length(order.y_left), Inf, order.y_left[index.1])
        index.y2_left <- ifelse(index.2>length(order.y_left), Inf, order.y_left[index.2])
        # index - type 1
        index.y1_right <- ifelse(index.1>length(order.y_right), Inf, order.y_right[index.1])
        index.y2_right <- ifelse(index.2>length(order.y_right), Inf, order.y_right[index.2])
        
        # CDF - type 2
        index.y1.cdf_left <- ifelse(index.1 == 0, 0, ci.ub[cbind(1:dim(new.data)[1], index.1)])
        index.y2.cdf_left <- ifelse(index.2 > length(order.y_left), 1, ci.ub[cbind(1:dim(new.data)[1], index.2)])
        # CDF - type 1
        index.y1.cdf_right <- ifelse(index.1 == 0, 0, ci.ub[cbind(1:dim(new.data)[1], index.1)])
        index.y2.cdf_right <- ifelse(index.2 > length(order.y_right), 1, ci.ub[cbind(1:dim(new.data)[1], index.2)])
        
        # quantile - type 2
        quantile_left <- ifelse(index.1==index.2, index.y1_left,
                                (index.y2_left - index.y1_left) / (index.y2.cdf_left - index.y1.cdf_left) *
                                  (probs[i] - index.y1.cdf_left) + index.y1_left)
        quantile_left <- ifelse(is.infinite(quantile_left), max(order.y), quantile_left)
        # quantile - type 1
        quantile_right <- ifelse(index.1==index.2, index.y1_right,
                                 (index.y2_right - index.y1_right)/(index.y2.cdf_right - index.y1.cdf_right) *
                                   (probs[i] - index.y1.cdf_right) + index.y1_right)
        quantile_right <- ifelse(is.infinite(quantile_right), max(order.y), quantile_right)
        
        # weight
        all_weights <- c(0, (ci.ub[-length(ci.ub)] - ci.ub[1]) / (ci.ub[length(ci.ub) - 1] - ci.ub[1]), 1)
        weight <- all_weights[index.1]
        # weighted average
        quantile.ub[,i] <- weight * quantile_left + (1 - weight) * quantile_right
        
        
      })
      
    }
    
    result <- list(est = quantile,
                   lb = quantile.ub,
                   ub = quantile.lb)
  }
  return(result)
}



#' Calculate expectations conditional on covariates
#'
#' This functions calculates the expectations conditional on covariates based on the fitted model
#' and new data.
#'
#' @param mod the model
#' @param cont_y the continuous response data
#' @param new.data the new data
#' @param se if confidence intervals needed (default = TRUE)
#' @return A list containing the following components:
#' @return \item{est}{a vector of estimated condtional expectations}
#' @return \item{se}{a vector of estimated standard errors}
#' @return \item{lb}{a vector of estimated lower bounds of 95\% confidence intervals}
#' @return \item{ub}{a vector of estimated upper bounds of 95\% confidence intervals}
#' @export
mean_cpmgee <- function(mod, cont_y, new.data, se=TRUE){
  # inputs
  coef <- mod$coefficients
  cov.par <- mod$robust.var
  order.y <- sort(unique(cont_y))
  alpha <- coef[1:(length(order.y)-1)]
  beta <- coef[length(order.y):length(coef)]
  cumprob <- function(x) 1 / (1 + exp(-x))
  
  n.alpha <- length(order.y) - 1
  xb <- as.matrix(new.data) %*% beta
  m.alpha <- alpha
  lb <- t(outer(m.alpha, xb, "+")[,,1])
  m.s <- cumprob(lb)
  m.f <- t(apply(m.s, 1, FUN=function(x) c(x[1:n.alpha], 1) - c(0, x[1:n.alpha])))
  m.mean <- apply(m.f, 1, FUN=function(x) sum(x*order.y))
  
  if(se){
    deriv <- function(x) exp(-x)/(1+exp(-x))^2
    
    dmean.dalpha <- t(apply(deriv(lb),
                            1, FUN=function(x) x*(order.y[2:length(order.y)] - order.y[1:n.alpha])))
    dmean.dbeta <- apply(dmean.dalpha, 1, sum)*as.matrix(new.data)
    dmean.dtheta <- cbind(dmean.dalpha, dmean.dbeta)
    mean.var <- diag(dmean.dtheta %*% cov.par %*% t(dmean.dtheta))
    mean.se <- sqrt(mean.var)
    result <- cbind(m.mean, mean.se)
    ci <- t(apply(result, 1, FUN=function(x) c(x[1] - qnorm(0.975) * x[2],
                                               x[1] + qnorm(0.975) * x[2])))
    result <- as.data.frame(cbind(result, ci))
    colnames(result) <- c("est", "se", "lb", "ub")
    result <- list(est = matrix(result$est, ncol=1),
                   se = matrix(result$se, ncol=1),
                   lb = matrix(result$lb, ncol=1),
                   ub = matrix(result$ub, ncol=1))
  } else{
    result <- list(est = m.mean)
  }
  
  return(result)
}

