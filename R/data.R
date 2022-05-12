#' Simulated clustered data
#'
#' A dataset with 100 clusters and each cluster has 6 observations.
#'
#' @format A data frame with 600 rows and 4 variables:
#' \describe{
#'   \item{y}{response data}
#'   \item{x}{time-invariant covariate (normally distributed)}
#'   \item{t}{time-varying covariate (0 to 1 for each cluster)}
#'   \item{id}{id for each cluster}
#' }
"data"