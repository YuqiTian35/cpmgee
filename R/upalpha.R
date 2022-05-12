#' Update the association parameter
#' 
#' This function updates alpha, the association parameter
#' 
#' @param hgmat a list of H and G matrices (see reference for repolr for more details)
#' @param alpha a numeric value for current alpha value
#' @param diffmeth a character string specifying the method used for estimation of alpha
#' @param h a numeric value of finite differencing
#' @return A list containing the following components:
#' @return \item{gvb}{the first derivative of generalized variance at convergence}
#' @return \item{ggvb}{the second derivative of generalized variance at convergence}
#' @return \item{alpha}{the updated assocation parameter}


upalpha <- function(hgmat, alpha, diffmeth, h){
  if(diffmeth == 'analytic'){
    # inputs
    hmat <- hgmat$hmat
    ghmat <- hgmat$ghmat
    gghmat <- hgmat$gghmat
    gmat <- hgmat$gmat
    ggmat <- hgmat$ggmat
    gggmat <- hgmat$gggmat
    
    # new variables
    ncoeff <- ncol(hmat)
    
    # eigen analysis
    h_eigen <- eigen(hmat)
    heigval <- h_eigen$values
    heigvec <- h_eigen$vectors
    g_eigen <- eigen(gmat)
    geigval <- g_eigen$values
    geigvec <- g_eigen$vectors
    
    # consturct sparse matrix
    cprodh <- Matrix::crossprod(heigvec, ghmat %*% heigvec)
    gldethmat <- sum(Matrix::diag(cprodh))
    cprodg <- Matrix::crossprod(geigvec, ggmat %*% geigvec)
    gldetgmat <- sum(Matrix::diag(cprodg))
    
    heigdenom <- matrix(rep(heigval, ncoeff), nrow=ncoeff) - matrix(rep(heigval, each=ncoeff), nrow=ncoeff)
    geigdenom <- matrix(rep(geigval, ncoeff), nrow=ncoeff) - matrix(rep(geigval, each=ncoeff), nrow=ncoeff)
    
    # hmat - bdiag
    fhnum <- Matrix::bdiag(lapply(1:ncoeff, function(i){
      vec_i <- heigvec[,i]
      return(Matrix::tcrossprod(vec_i, vec_i))
    }))
    # gmat - bdiag
    fgnum <- Matrix::bdiag(lapply(1:ncoeff, function(i){
      vec_i <- geigvec[,i]
      return(Matrix::tcrossprod(vec_i, vec_i))
    }))
    
    # structure to store intermediate results
    intmat <- do.call(cbind, replicate(ncoeff, Matrix::Diagonal(ncoeff), simplify=FALSE))
    dd <- sapply(1:ncoeff, function(k){
      ## hmat
      vkheigdenom <- rep(heigdenom[,k], each=ncoeff)
      qh <- which(vkheigdenom !=0 )
      vkheigdenom[qh] <- 1 / vkheigdenom[qh]
      bhmat <- Matrix::Diagonal(n=ncoeff^2, vkheigdenom)
      fhmat <- Matrix::tcrossprod(intmat %*% bhmat %*% fhnum, intmat)
      dheig <- mat_minus(gghmat, 2 * ghmat %*% fhmat %*% ghmat)
      
      ## gmat
      vkgeigdenom <- rep(geigdenom[,k], each=ncoeff)
      qg <- which(vkgeigdenom !=0 )
      vkgeigdenom[qg] <- 1 / vkgeigdenom[qg]
      bgmat <- Matrix::Diagonal(n=ncoeff^2, vkgeigdenom)
      fgmat <- Matrix::tcrossprod(intmat %*% bgmat %*% fgnum, intmat)
      dgeig <- mat_minus(gggmat, 2 * ggmat %*% fgmat %*% ggmat)
      
      return(c((t(heigvec[,k]) %*% dheig %*% heigvec[,k])[1,1], 
               (t(geigvec[,k]) %*% dgeig %*% geigvec[,k])[1,1]))
    })
    ddheig <- dd[1,]
    ddgeig <- dd[2,]
  
    
    # calculate gvb and ggvb
    ggldethmat <- sum(ddheig / heigval) - sum((Matrix::diag(cprodh) / heigval)^2)
    ggldetgmat <- sum(ddgeig / geigval) - sum((Matrix::diag(cprodg) / geigval)^2)
    gvb <- gldetgmat - 2 * gldethmat
    ggvb <- ggldetgmat - 2 * ggldethmat
    
    # update alpha
    phi <- log(alpha) - log(1 - alpha)
    gphi <- exp(phi) / (1 + exp(phi))^2
    ggphi <- exp(phi) * (1 - exp(phi)) / (1 + exp(phi))^2
    gvbphi <- gphi * gvb
    ggvbphi <- ggphi * gvb + gphi^2 * ggvb
    phiinc <- gvbphi / abs(ggvbphi)
    if(abs(phiinc) > 1){
      if(phiinc > 1){
        nphi <- phi - 1
      }else{
        nphi <- phi + 1
      }
    }else{
      nphi <- phi - phiinc
    }
    
    nalpha <- exp(nphi) / (1 + exp(nphi))
    if(nalpha > 0.95){
      nalpha <- 0.95
    }
    if(nalpha < 0.05){
      nalpha <- 0.05
    }
    
    # output
    return(list(gvb = gvbphi, ggvb = ggvbphi, alpha = nalpha))
    
  } else if(diffmeth == 'numeric'){
    # inputs
    hmat <- hgmat$hmat
    hmatl <- hgmat$lhmat
    hmatu <- hgmat$uhmat
    gmat <- hgmat$gmat
    gmatl <- hgmat$lgmat
    gmatu <- hgmat$ugmat
    
    # varance matrices
    ldetvcov <- log(prod(eigen(gmat)$values)) - 2 * log(prod(eigen(hmat)$values))
    ldetvcovl <- log(prod(eigen(gmatl)$values)) - 2 * log(prod(eigen(hmatl)$values))
    ldetvcovu <- log(prod(eigen(gmatu)$values)) - 2 * log(prod(eigen(gmatu)$values))
    
    ## 
    gvbphi <- (ldetvcovu - ldetvcovl) / (2 * h)
    ggvbphi <- (ldetvcovu - 2 * ldetvcov + ldetvcovl) / h^2
    phi <- log(alpha) - log(1 - alpha)
    phiinc <- gvbphi / abs(ggvbphi)
    
    if(abs(phiinc) > 1){
      if(phiinc > 1){
        nphi <- phi - 1
      }else{
        nphi <- phi + 1
      }
    }else{
      nphi <- phi - phiinc
    }
    
    nalpha <- exp(nphi) / (1 + exp(nphi))
    if(nalpha > 0.95){
      nalpha <- 0.95
    }
    if(nalpha < 0.05){
      nalpha <- 0.05
    }
    
    # output
    return(list(gvb = gvbphi, ggvb = ggvbphi, alpha = nalpha))
    
  }else{
    return(0)
  }
}
