\name{cpmgee-package}
\alias{cpmgee-package}
\docType{package}
\title{CPMs for Clustered Continuous Response Variables Based on GEE Methods}
\description{We extend a widely used ordinal regression model, the cumulative probability model (CPM), to fit clustered continuous response variables based on generalized estimating equation (GEE) methods for ordinal responses. With our approaches, estimates of marginal parameters, cumulative distribution functions (CDFs), expectations, and quantiles conditional on covariates can be obtained without pre-transformation of the potentially skewed continuous response data.
}
\details{aaa}
\author{Yuqi Tian
Maintainer: Chun Li <cli77199@usc.edu>
}
\references{
Tian et al. "Analyzing clustered continuous response variables with ordinal regresson models." (2022) (to be submitted)

Parsons, N. (2017). repolr: an R package for fitting proportional-odds models to repeated ordinal scores. R package version 3.4 https://CRAN.R-project.org/package=repolr

Parsons, N., et al. "Repeated measures proportional odds logistic regression analysis of ordinal score data in the statistical software package R." Computational Statistics & Data Analysis 53.3 (2009): 632-641.


Harrell, F. (2020). rms: Regression modeling strategies. R package version 6.1.0. https://CRAN.R-project.org/package=rms
}
\keyword{package}
\seealso{
\code{\link{cdf_cpmgee}}, \code{\link{quantile_cpmgee}}, \code{\link{mean_cpmgee}}
}
\examples{
data(data)
# independence working correlation structure
mod_cpmgee_ind <- cpmgee(formula = y ~ x + t, data = data, categories = length(unique(data$y)),
subjects = 'id', times = 1:6, corr.mod = 'independence', alpha = 0.5)
# exchangeable working correlation structure
mod_cpmgee_ex <- cpmgee(formula = y ~ x + t, data = data, categories = length(unique(data$y)),
subjects = 'id', times = 1:6, corr.mod = 'exchangeable', alpha = 0.5)
# new data
new.data <- data.frame(x = c(0,1), t = 0.2)
# conditional quantities
mean_ind <- mean_cpmgee(mod_cpmgee_ind, data$y, new.data)
median_ind <- quantile_cpmgee(mod_cpmgee_ind, data$y, new.data, 0.5)
cdf_ind <- cdf_cpmgee(mod_cpmgee_ind, data$y, new.data, 0.5)
}
