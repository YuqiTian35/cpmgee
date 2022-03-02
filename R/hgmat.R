# hgmat R version 
# can handle missing data

hgmat <- function(mod, smat, X, modtype, diffmeth, alpha, corrmod, h){
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
    id <- c(mod$id)
  }  
  
  
  # calculate new variables
  n_id <- length(unique(id))
  cats1 <- nrow(ismat)
  table_id <- table(id)
  ntimes <- table_id / cats1
  
  # start and end index for each subject
  end_ind <- cumsum(table_id)
  start_ind <- end_ind - table_id + 1
  
  # create icormat for differrent ntimes[i]
  distinct_time <- unique(ntimes)
  icormat_dict <- lapply(1:length(distinct_time), function(i){
    cmat <- cmat(ctimes = 1:distinct_time[i], alpha = alpha, corrmod = corrmod, 
                 diffmeth = diffmeth, h = h)
    icmat <- cmat$icmat
    if(diffmeth == "analytic"){
      gicmat <- cmat$gicmat
      ggicmat <- cmat$ggicmat
    } else if(diffmeth == "numeric"){
      gicmat <- cmat$licmat
      ggicmat <- cmat$uicmat
    }
    
    # icor, gicor, ggicor
    icor <- as(kronecker(icmat, ismat), 'dgCMatrix')
    gicor <- as(kronecker(gicmat, ismat), 'dgCMatrix')
    ggicor <- as(kronecker(ggicmat, ismat), 'dgCMatrix')
    
    return(list(icor = icor, gicor = gicor, ggicor = ggicor))
  })
  names(icormat_dict) <- distinct_time
  
  hmat <- ghmat <- gghmat <- gmat <- ggmat <- gggmat <- Matrix::Matrix(0, nrow=ncol(Xmat), ncol=ncol(Xmat), sparse = TRUE)
  icors <- vector("list", length = n_id)
  for(i in 1:n_id){
    # new variables
    linpred_i <- linpred[start_ind[i]:end_ind[i]]
    dmat_i <- Matrix::Diagonal(n = length(linpred_i), exp(linpred_i) / (1 + exp(linpred_i))^2)
    Xmat_i <- Xmat[start_ind[i]:end_ind[i],]
    ddmat_i <- Matrix::crossprod(Xmat_i, dmat_i)
    
    # icor, gicor, ggicor
    cor_list <- icormat_dict[[as.character(ntimes[i])]]
    icor <- cor_list$icor
    gicor <- cor_list$gicor
    ggicor <- cor_list$ggicor
    
    # residuals
    residuals_i <- residuals[start_ind[i]:end_ind[i]]
    res_sqr <- Matrix::tcrossprod(residuals_i, residuals_i)
    
    # between category covariance blocks
    varmat_i <- varmat[start_ind[i]:end_ind[i]]
    vmat_i <- Matrix::Diagonal(n=length(varmat_i), 1 / varmat_i)
    start_i <- ((1:ntimes[i]) - 1) * cats1 + 1
    end_i <- (1:ntimes[i]) * cats1
    # }))
    insvcov_i <- rep(0, cats1 * cats1 * ntimes[i])
    for(t in 1:ntimes[i]){
      varmat_it <- Matrix::Diagonal(n = cats1, varmat_i[start_i[t]:end_i[t]])
      insvcov_i[((t-1)*cats1^2+1) : (t*cats1^2)] <- as.vector(varmat_it %*% ssmat %*% varmat_it)
    }
    
    qfullresmat <- res_sqr
    diagonals::fatdiag(qfullresmat, size = cats1) <- insvcov_i
    
    if(diffmeth == 'analytic'){
      # hmat, ghmat, gghmat
      vcovmat_i <- vmat_i %*% icor %*% vmat_i
      gvcovmat_i <- vmat_i %*% gicor %*% vmat_i
      ggvcovmat_i <- vmat_i %*% ggicor %*% vmat_i
      
      dv_i <- ddmat_i %*% vcovmat_i
      dgv_i <- ddmat_i %*% gvcovmat_i
      dggv_i <- ddmat_i %*% ggvcovmat_i
      
      hmat_i <- Matrix::tcrossprod(dv_i, ddmat_i)
      ghmat_i <- Matrix::tcrossprod(dgv_i, ddmat_i)
      gghmat_i <- Matrix::tcrossprod(dggv_i, ddmat_i)
      
      # gmat, ggmat, gggmat
      dvq_i <- dv_i %*% qfullresmat
      dgvq_i <- dgv_i %*% qfullresmat
      dggvq_i <- dggv_i %*% qfullresmat
      
      gmat_i <- Matrix::tcrossprod(dvq_i, dv_i)
      ggmat_i <- Matrix::tcrossprod(dgvq_i, dv_i) +  Matrix::tcrossprod(dvq_i, dgv_i)
      gggmat_i <- Matrix::tcrossprod(dggvq_i, dv_i) + 2 * Matrix::tcrossprod(dgvq_i, dgv_i) + Matrix::tcrossprod(dvq_i, dggv_i)
    } else if(diffmeth == "numeric"){
      # hmat, ghmat, gghmat
      vcovmat_i <- vmat_i %*% icor %*% vmat_i
      gvcovmat_i <- vmat_i %*% gicor %*% vmat_i
      ggvcovmat_i <- vmat_i %*% ggicor %*% vmat_i
      
      dv_i <- ddmat_i %*% vcovmat_i
      dgv_i <- ddmat_i %*% gvcovmat_i
      dggv_i <- ddmat_i %*% ggvcovmat_i
      
      hmat_i <- dv_i %*% t(ddmat_i)
      ghmat_i <- dgv_i %*% t(ddmat_i)
      gghmat_i <- dggv_i %*% t(ddmat_i)

      # gmat, ggmat, gggmat
      gmat_i <- dv_i %*% qfullresmat %*% t(dv_i)
      ggmat_i <- dgv_i %*% qfullresmat %*% t(dgv_i)
      gggmat_i <- dggv_i %*% qfullresmat %*% t(dggv_i)

    }
    hmat <- mat_add(hmat, hmat_i)
    ghmat <- mat_add(ghmat, ghmat_i)
    gghmat <- mat_add(gghmat, gghmat_i)
    gmat <- mat_add(gmat, as(gmat_i, "dgCMatrix"))
    ggmat <- mat_add(ggmat, as(ggmat_i, "dgCMatrix"))
    gggmat <- mat_add(gggmat, as(gggmat_i, "dgCMatrix"))
    icors[[i]] <- icor
  }
  
  if(diffmeth == 'analytic'){
    return(list('icormat' = Matrix::bdiag(icors),
                'hmat' = hmat,
                'ghmat' = ghmat,
                'gghmat' = gghmat,
                'gmat' = gmat,
                'ggmat' = ggmat,
                'gggmat' = gggmat))
  }else if(diffmeth == 'numeric'){
    return(list('icormat' = Matrix::bdiag(icors),
                'hmat' = hmat,
                'lhmat' = ghmat,
                'uhmat' = gghmat,
                'gmat' = gmat,
                'lgmat' = ggmat,
                'ugmat' = gggmat))
  }else{
    return(0)
  }
  
}
