#' Obtain the working correlation matrix
#' 
#' 
#' @param mod the fitted model with independence working correlation structure
#' @param smat the matrix of correlations among elements in $Z_{it}$
#' @param modtype a character string specifying the model type (CPM or GLM)
#' @param diffmeth a character string specifying the method used for estimation of alpha
#' @param alpha the association parameter
#' @param corrmod the specified working correlation structure
#' @param h a numeric value of options to control the fitting algorithm
#' @return the working correlation matrix

icormat <- function(mod, smat, modtype = "gee", diffmeth, alpha, corrmod, h){
  ismat <- smat$ismat
  ssmat <- smat$smat
  
  if (modtype == "glm"){
    id <- mod$data$subjects
  } else if (modtype == "gee") {
    id = mod$id
  }  
  
  # number of times / 
  n_id <- length(unique(id)) # number of subjects
  cats1 <- nrow(ismat)
  table_id <- table(id)
  ntimes <- table_id / cats1
  
  end_ind <- cumsum(table_id)
  start_ind <- end_ind - table_id + 1
  
  
  # create icormat for differrent ntimes[i]
  distinct_time <- unique(ntimes)
  icormat_dict <- lapply(1:length(distinct_time), function(i){
    cmat <- cmat(ctimes = 1:distinct_time[i], alpha = alpha, corrmod = corrmod, 
                 diffmeth = diffmeth, h = h)
    icmat <- cmat$icmat
    icor <- kronecker(icmat, ismat)
    return(icor)
  })
  names(icormat_dict) <- distinct_time
  
  # create icormat
  icormat <- Matrix::bdiag(lapply(1:n_id, function(i){
    icor <- icormat_dict[[as.character(ntimes[i])]]
    return(icor)
  }))
  
  return(icormat)
}


