# hgmat_cov R version 
# can handle missing data

hgmat_cov <- function(mod, smat, X, modtype, diffmeth, alpha, corrmod, h){
  # input
  Xmat <- X$design
  linpred <- mod$linear.predictors
  fitted <- mod$fitted.values
  y <- mod$y
  residuals <- y - fitted
  # cmat, gicmat, ggicmat
  ismat <- smat$ismat
  ssmat <- smat$smat
  varmat <- sqrt(fitted * (1 - fitted))
  
  if (modtype == "glm"){
    id <- mod$data$subjects
  } else if (modtype == "gee") {
    id = c(mod$id)
  }  
  
  # start and end index for each subject
  end_ind <- cumsum(table(id))
  start_ind <- c(1, end_ind[-length(end_ind)]+1)
  names(start_ind) <- names(end_ind)
  
  # calculate new variables 
  n_id <- length(unique(id))
  cats1 <- nrow(ismat)
  table_id <- table(id)
  ntimes <- table_id / cats1

  # create icormat for differrent ntimes[i]
  distinct_time <- unique(ntimes)
  icormat_dict <- lapply(1:length(distinct_time), function(i){
    cmat <- cmat(ctimes = 1:distinct_time[i], alpha = alpha, corrmod = corrmod, 
                 diffmeth = diffmeth, h = h)
    icmat <- cmat$icmat
    icor <- as(kronecker(icmat, ismat), 'dgCMatrix')
    return(icor)
  })
  names(icormat_dict) <- distinct_time
  
  hmat <- gmat <- Matrix::Matrix(0, nrow=ncol(Xmat), ncol=ncol(Xmat), sparse = TRUE)
  
  for(i in 1:n_id){
    # new variables
    linpred_i <- linpred[start_ind[i]:end_ind[i]]
    dmat_i <- Matrix::Diagonal(n = length(linpred_i), exp(linpred_i) / (1 + exp(linpred_i))^2)
    Xmat_i <- Xmat[start_ind[i]:end_ind[i],]
    ddmat_i <- Matrix::crossprod(Xmat_i, dmat_i)
    
    icor <- icormat_dict[[as.character(ntimes[i])]]
    
    # residuals
    residuals_i <- residuals[start_ind[i]:end_ind[i]]
    res_sqr <- Matrix::tcrossprod(residuals_i, residuals_i)
    
    # between category covariance blocks
    varmat_i <- varmat[start_ind[i]:end_ind[i]]
    vmat_i <- Matrix::Diagonal(n=length(varmat_i), 1 / varmat_i)
    start_i <- ((1:ntimes[i]) - 1) * cats1 + 1
    end_i <- (1:ntimes[i]) * cats1
    insvcov_i <- rep(0, cats1 * cats1 * ntimes[i])
    for(t in 1:ntimes[i]){
      varmat_it <- Matrix::Diagonal(n = cats1, varmat_i[start_i[t]:end_i[t]])
      insvcov_i[((t-1)*cats1^2+1) : (t*cats1^2)] <- as.vector(varmat_it %*% ssmat %*% varmat_it)
    }
    
    qfullresmat <- res_sqr
    diagonals::fatdiag(qfullresmat, size = cats1) <- insvcov_i
    
    vcovmat_i <- vmat_i %*% icor %*% vmat_i
    
    dv_i <- ddmat_i %*% vcovmat_i
    hmat_i <- Matrix::tcrossprod(dv_i, ddmat_i)
    gmat_i <- Matrix::tcrossprod(dv_i %*% qfullresmat, dv_i)
    
    hmat <- mat_add(hmat, hmat_i)
    gmat <- mat_add(gmat, as(gmat_i, "dgCMatrix"))
  }
  
  return(list('hmat' = hmat,
              'gmat' = gmat))
  
  
}

