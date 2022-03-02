# CPM GEE

cpmgee <- function(formula, subjects, data, times, categories, corr.mod = "independence", 
                              alpha = 0.5, po.test = FALSE, fixed = FALSE, poly = NULL, space = NULL, fit.opt = rep(NA, 5)){
  # start
  call <- match.call()
  
  # set-up
  corr.mods <- c("ar1", "uniform", "independence")
  icorr.mod <- as.integer(match(corr.mod, corr.mods, -1))
  if (icorr.mod < 1){stop("corr.mod: set to be independence, ar1 or uniform")}
  rcorr.mod <- corr.mod
  # diffmeths <- c("analytic", "numeric")
  # idiffmeth <- as.integer(match(diffmeth, diffmeths, -1))
  # if (idiffmeth < 1){stop("diffmeth: set to be numeric or analytic")}
  diffmeth <- 'analytic'
  initial <- 'orm'
  if(times[1] != 1){stop("times: times should be vector with first value set to 1")}
  alpha <- as.double(alpha)
  po.test <- as.logical(po.test)
  fixed <- as.logical(fixed)
  set.fit.opt <- c(cmaxit = 10, omaxit = 5, ctol = 0.001, otol = 0.00001, h = 0.01)
  set.fit.opt[which(is.na(fit.opt) == FALSE)] <- fit.opt[which(is.na(fit.opt) == FALSE)]
  if(set.fit.opt[1] < 4){set.fit.opt[1] <- 4}
  if(alpha < 0.05 | alpha > 0.95){stop("alpha: invalid correlation parameter")}
  categories <- as.integer(categories)
  categories1 <- categories - 1
  subjects <- as.character(subjects)
  isubject <- as.integer(match(subjects, names(as.data.frame(data)), -1))
  if (isubject < 1){stop("subjects: unknown subject name")}
  orig.formula <- as.formula(formula)
  # if(corr.mod == "independence" | length(times) == 1){alpha <- 0;   fixed <- TRUE;   corr.mod <- "uniform"}
  # if(length(times) == 1){diffmeth <- "analytic"}
  # if(is.null(poly) == FALSE){po.test <- FALSE; poly <- as.integer(poly);  initial = "glm"}
  if(is.null(poly) == FALSE){po.test <- FALSE; poly <- as.integer(poly)}
  if(sum(diff(times) > 0) != (length(times) - 1)){stop("times: invalid vector of times")}
  if(is.null(space) != TRUE){if(sum(diff(space) > 0) != categories1){stop("space: invalid vector of spacings")}}
  if(is.null(space) != TRUE){if(space[1] != 1){stop("space: space should be vector with first value set to 1")}}
  
  # convert id to numeric
  data$subjects.ord <- as.numeric(factor(data[,subjects], levels = unique(data[,subjects])))
  data <- data[order(data$subjects.ord),]
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

  # model matrix
  exdata <- ord.expand(space = space, formula = formula.ord, times = times, poly = poly,
                       data = data, subjects = 'subjects.ord', categories = categories)
  formula <- exdata$formula
  dat <- list(data = exdata$data)
  # Xmat <- list(design = Matrix::Matrix(model.matrix(formula, data = exdata$data), sparse = TRUE)) 
  covariate_part <- model.matrix(formula.ord, data = exdata$data)[,-1]
  covariate_names <- colnames(covariate_part)
  intercept_names <- paste0('cuts', exdata$data$cuts[1:categories1])
  # intercept part for one observation
  if(is.null(poly) == FALSE){
    i_poly <- poly(exdata$data$pcuts, poly)
    poly_names <- paste0('poly(pcuts, ', poly, ')', 1:poly)
    Xmat <- list(design = as(cbind(1, i_poly, covariate_part), 'dgCMatrix'))
    # var.names <- c('(Intercept)', poly_names, covariate_names)
  }else{
    i_intercept <- Matrix::Diagonal(categories1)
    intercept_names <- paste0('cuts', exdata$data$cuts[1:categories1])
    # number of replicates
    v <- rep(1, nrow(data))
    intercept_part <- kronecker(v, i_intercept)
    Xmat <- list(design = cbind(intercept_part, covariate_part))
  }
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
      corr.mod <- 'independence',
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
  if(is.null(poly) == FALSE){
    polycuts <- as.numeric(polydesign %*% coeffs[1:(poly1)])
    coeffs <- c(polycuts, mod.ordgee$coefficients[(poly1+1):length(mod.ordgee$coefficients)])
    xsmat <- smat(coeff=polycuts)
  } else {
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
  }
  
  xsmat <- smat(coeff = mod.ordgee$coefficients[1:categories1])
  
  # covariance matrices
  xhgmat <- hgmat_cov(mod = mod.ordgee, smat = xsmat, 
                      X = Xmat, modtype = "gee", diffmeth = diffmeth, 
                      alpha = alpha, corrmod = corr.mod, h = set.fit.opt[5])
  inv_hmat <- Matrix::solve(xhgmat$hmat) 
  robust.var <- inv_hmat %*% xhgmat$gmat %*% inv_hmat
  if(is.null(poly) == FALSE){
    poly_robust.var <- robust.var
    polycuts.robust_int <- polydesign %*% robust.var[1:(poly1), 1:(poly1)] %*% Matrix::t(polydesign)
    polycuts.robust_cov <- robust.var[(poly1+1):nrow(robust.var), (poly1+1):nrow(robust.var)] 
    polycuts.robust_intcov <- polydesign %*% robust.var[1:(poly1), (poly1+1):nrow(robust.var)]
    robust.var <- as.matrix(rbind(cbind(polycuts.robust_int, polycuts.robust_intcov),
                                  cbind(t(polycuts.robust_intcov), polycuts.robust_cov)))
  } else{
    inv_hmat <- Matrix::solve(xhgmat$hmat) 
    robust.var <- inv_hmat %*% xhgmat$gmat %*% inv_hmat
    polycuts <- poly_robust.var <- NA
  }
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
                     poly.mod = list(poly = poly, polycuts = list(coeff = polycuts, robust.var = poly_robust.var)),
                     max.id = as.numeric(mod.ordgee$max.id),
                     id = as.numeric(mod.ordgee$id),
                     y = as.numeric(mod.ordgee$y),
                     linear.predictors = as.numeric(mod.ordgee$linear.predictors),
                     fitted.values = as.numeric(mod.ordgee$fitted.values),
                     coefficients = coeffs,
                     robust.var = as.matrix(robust.var),
                     fixed = fixed,
                     alpha = as.numeric(alpha),
                     fit.opt = set.fit.opt,
                     grad1 = as.numeric(xupalpha$gvb),
                     grad2 = as.numeric(xupalpha$ggvb))
}
