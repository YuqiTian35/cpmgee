#' CPMs for clustered continuous response data based on GEE methods for ordinal data
#'
#' This function function fits cumulative probability models (CPMs) clustered/longitudinal continuous response data
#' based on GEE methods for ordinal data.
#'
#' @param formula an R formula object
#' @param subjects a character string specifying the name of the subject variable
#' @param data a data frame including response data and covariates
#' @param corr.mod a character string specifying the working correlation structure
#'  ("independence", "exchangeable", and "ar1")
#' @param alpha an initial value for the association parameter in the range of 0.05 to 0.95.
#' @param fit.opt a vector of options to control the behavior of the fitting algorithm
#' @return A list containing the following components:
#' @return \item{max.id}{number of clusters}
#' @return \item{fitted.values}{a vector of the fitted values}
#' @return \item{linear.predictors}{a vector of linear predictors}
#' @return \item{coefficients}{a vector of interecept and regression parameters}
#' @return \item{robust.var}{the robust (sandwich) variance matrix}
#' @return \item{alpha}{the estimate of the association parameter}
#' @export
#'
#' @details CPMs are useful for the analysis of continuous response data which may need to be transformed 
#' prior to fitting standard regression models. CPMs are semi-parametric linear transformation models; 
#' they nonparametrically estimate the appropriate transformation as part of the fitting procedure.
#'
#' We propose two feasible and computationally efficient approaches to
#' fit CPMs for clustered continuous response variables with different working correlation structures
#' (independence, exchangeable and AR1 working correlation structures).
#'
#' CPMs with independence working correlation can be efficiently fit to clustered continuous responses
#' with thousands of distinct values based on CPMs and sandwich estimator for variances.
#'
#' To improve efficiency, CPMs with more complex working correlation structures (exchangeable and AR1)
#' can be fit with a one-step GEE estimator for repolr (repeated measures proportional odds logistic regression
#' proposed by Parsons). The number of distinct response values can be further reduced by
#' equal-quantile binning or rounding.
#'
#' Estimates of the mean, quantiles, and exceedance probabilities conditional on covariates (new data)
#' can be derived from the model fit.
#'
#' @seealso \code{\link{cdf_cpmgee}, \link{quantile_cpmgee}, \link{mean_cpmgee}}
#'
#'
#' @references
#' Tian et al. "Analyzing clustered continuous response variables with ordinal
#'  regression models." (2022) (to be submitted)
#' @references
#' Parsons, N. (2017). repolr: an R package for fitting proportional-odds models to repeated
#' ordinal scores. R package version 3.4 https://CRAN.R-project.org/package=repolr
#' @references
#' Harrell, F. (2020). rms: Regression modeling strategies. R package version 6.1.0. https://CRAN.R-project.org/package=rms
#' @references
#' Liu, Q., Shepherd, B. E., Li, C., & Harrell Jr, F. E. (2017). Modeling continuous response variables using ordinal regression. 
#' Statistics in Medicine, 36(27), 4316-4335.
#'
#'
#' @import Matrix
#' @import rms
#' @import diagonals
#' @importFrom Rcpp evalCpp
#' @importFrom methods as
#' @importFrom stats as.formula coef dnorm model.matrix qnorm terms terms.formula
#' @useDynLib cpmgee, .registration=TRUE
#' @exportPattern "^[[:alpha:]]+"
#'
#' @examples
#' data(data)
#'
#' # independence working correlation structure
#' mod_cpmgee_ind <- cpmgee(formula = y ~ x + t, data = data, 
#' subjects = 'id', corr.mod = 'independence')
#'
#' # exchangeable working correlation structure
#' mod_cpmgee_ex <- cpmgee(formula = y ~ x + t, data = data,
#' subjects = 'id', corr.mod = 'exchangeable', alpha = 0.5)
#'
#' # new data
#' new_data <- data.frame(x = c(0,1), t = 0.2)
#'
#' # conditional quantities for independence working correlation structure
#' mean_ind <- mean_cpmgee(mod_cpmgee_ind, data$y, new_data)
#' median_ind <- quantile_cpmgee(mod_cpmgee_ind, data$y, new_data, 0.5)
#' cdf_ind <- cdf_cpmgee(mod_cpmgee_ind, data$y, new_data, 5)
#'
#' # conditional quantities for exchangeable working correlation structure
#' mean_ex <- mean_cpmgee(mod_cpmgee_ex, data$y, new_data)
#' median_ex <- quantile_cpmgee(mod_cpmgee_ex, data$y, new_data, 0.5)
#' cdf_ex <- cdf_cpmgee(mod_cpmgee_ex, data$y, new_data, 5)

cpmgee <- function(formula, subjects, data, corr.mod = "independence",
                   alpha = 0.5, fit.opt = rep(NA, 5)){
  # start
  call <- match.call()

  # set-up
  corr.mods <- c("ar1", "exchangeable", "independence")
  icorr.mod <- as.integer(match(corr.mod, corr.mods, -1))
  if (icorr.mod < 1){stop("corr.mod: set to be independence, ar1 or exchangeable")}
  rcorr.mod <- corr.mod
  diffmeth <- 'analytic'
  initial <- 'orm'
  alpha <- as.double(alpha)
  set.fit.opt <- c(cmaxit = 10, omaxit = 5, ctol = 0.001, otol = 0.00001, h = 0.01)
  set.fit.opt[which(is.na(fit.opt) == FALSE)] <- fit.opt[which(is.na(fit.opt) == FALSE)]
  if(set.fit.opt[1] < 4){set.fit.opt[1] <- 4}
  if(alpha < 0.05 | alpha > 0.95){stop("alpha: invalid correlation parameter")}
  subjects <- as.character(subjects)
  isubject <- as.integer(match(subjects, names(as.data.frame(data)), -1))
  if (isubject < 1){stop("subjects: unknown subject name")}
  orig.formula <- as.formula(formula)
  
  # convert id to numeric
  data$subjects.ord <- as.numeric(factor(data[,subjects], levels = unique(data[,subjects])))
  data <- data[order(data$subjects.ord),]
  
  # max cluster size
  times <- seq(1, max(table(data$subjects.ord)))
  
  # convert response to ordinal variable
  y <- data[,all.vars(formula)[1]]
  if(is.factor(y)){
    ylevels <- levels(y)
    y.ord <- unclass(y)
  }else{
    ylevels <- sort(unique(y))
    y.ord <- match(y, ylevels)
  }
  data$y.ord <- y.ord
  formula.ord <- as.formula(paste('y.ord', '~', paste(all.vars(formula)[-1], collapse = ' + ')))

  # categories 
  categories <- max(y.ord)
  categories1 <- categories - 1
  
  # model matrix
  exdata <- ord.expand(formula = formula.ord, times = times,
                       data = data, subjects = 'subjects.ord', categories = categories)
  formula <- exdata$formula
  dat <- list(data = exdata$data)
  covariate_part <- model.matrix(formula.ord, data = exdata$data)[,-1]
  covariate_names <- colnames(covariate_part)
  intercept_names <- paste0('cuts', exdata$data$cuts[1:categories1])
  # intercept part for one observation
  i_intercept <- Matrix::Diagonal(categories1)
  intercept_names <- paste0('cuts', exdata$data$cuts[1:categories1])
  # number of replicates
  v <- rep(1, nrow(data))
  intercept_part <- kronecker(v, i_intercept)
  Xmat <- list(design = cbind(intercept_part, covariate_part))
  var.names <- c(intercept_names, covariate_names)

  mod <- rms::orm(formula = formula.ord, data = data, x = TRUE, y = TRUE) # default logit link
  inv.logit <- function(x) 1 / (1 + exp(-x))

  # independence working correlation
  if(corr.mod == 'independence' | length(times) == 1){
    mod_robust <- rms::robcov(fit = mod, cluster = data[,subjects])
    coeffs <- -mod_robust$coefficients
    robust.var <- mod_robust$var
    rownames(robust.var) <- colnames(robust.var) <- var.names
    names(coeffs) <- var.names

    return(list(
      title = "CPM GEE",
      call = call,
      data = call[["data"]],
      subjects = subjects,
      formula = mod_robust$sformula,
      orig.formula = mod_robust$sformula,
      corr.mod = 'independence',
      times = times,
      categories = categories,
      id = data[,subjects],
      max.id = length(unique(data[,subjects])),
      y = exdata$data[,all.vars(formula)[1]],
      linear.predictors = as.vector(Xmat$design %*% mod[['coefficients']]),
      fitted.values = inv.logit(mod$linear.predictors),
      coefficients = coeffs,
      robust.var = robust.var
    ))
  }

  mod[['coefficients']] <- -mod[['coefficients']]
  mod[['linear.predictors']] <- as.vector(Xmat$design %*% mod[['coefficients']])
  mod[['fitted.values']] <- inv.logit(mod$linear.predictors)
  mod[['data']] <- exdata$data
  mod[['y']] <- exdata$data[,all.vars(formula)[1]]
  xsmat <- smat(coeff = mod$coefficients[1:categories1])

  # 2. update alpha
  xhgmat <- hgmat(mod = mod, smat = xsmat, X = Xmat,
                  modtype = "glm", diffmeth = diffmeth, alpha = alpha,
                  corrmod = corr.mod, h = set.fit.opt[5])

  xupalpha <- upalpha(hgmat = xhgmat, alpha = alpha, diffmeth = diffmeth, h = set.fit.opt[5])
  alpha <- xupalpha$alpha

  # 3. update beta
  xicormat <- list(irmat = icormat(mod=mod, smat=xsmat, modtype='glm', diffmeth = diffmeth,
                                   alpha=alpha, corrmod = corr.mod, h = set.fit.opt[5]))
  mod.ordgee <- ordgee(mod = mod, icormat = xicormat, X = Xmat,
                          ctimes = times, categories = categories,
                          omaxit = as.integer(set.fit.opt[2]), otol = as.double(set.fit.opt[4]))
  coeffs <- mod.ordgee$coefficients

  polycuts <- NA
  # force increasing for intercepts
  decreasing_ind <- which(diff(coeffs[1:categories1]) <= 0)
  if(length(decreasing_ind) > 0){
    for(ind in decreasing_ind){
      # if the differences are extremely small -> use the mid value
      if(coeffs[ind+1] - coeffs[ind-1] < 1e-5){
        coeffs[ind] <- (coeffs[ind+1] + coeffs[ind-1]) / 2
      }else{
        coeffs[ind] <- coeffs[ind-1] + 1e-5
      }
    }
    ## putting the cofficients back to model (recalculate in fitted/linpred etc.)
    mod.ordgee <- fixmod(mod.ordgee, coeffs, Xmat)
    coeffs <- mod.ordgee$coefficients
  }

  xsmat <- smat(coeff = mod.ordgee$coefficients[1:categories1])

  # covariance matrices
  xhgmat <- hgmat_cov(mod = mod.ordgee, smat = xsmat,
                      X = Xmat, modtype = "gee", diffmeth = diffmeth,
                      alpha = alpha, corrmod = corr.mod, h = set.fit.opt[5])
  inv_hmat <- Matrix::solve(xhgmat$hmat)
  robust.var <- inv_hmat %*% xhgmat$gmat %*% inv_hmat

  inv_hmat <- Matrix::solve(xhgmat$hmat)
  robust.var <- inv_hmat %*% xhgmat$gmat %*% inv_hmat
  polycuts <- poly_robust.var <- NA

  rownames(robust.var) <- colnames(robust.var) <- var.names
  coeffs <- as.numeric(coeffs); names(coeffs) <- var.names

  # output
  cpmgee.mod <- list(title = "CPM GEE",
                     call = call,
                     data = call[["data"]],
                     subjects = subjects,
                     formula = formula,
                     orig.formula = orig.formula,
                     corr.mod = rcorr.mod,
                     times = times,
                     categories = categories,
                     max.id = as.numeric(mod.ordgee$max.id),
                     id = as.numeric(mod.ordgee$id),
                     y = as.numeric(mod.ordgee$y),
                     linear.predictors = as.numeric(mod.ordgee$linear.predictors),
                     fitted.values = as.numeric(mod.ordgee$fitted.values),
                     coefficients = coeffs,
                     robust.var = as.matrix(robust.var),
                     alpha = as.numeric(alpha),
                     fit.opt = set.fit.opt)
}
