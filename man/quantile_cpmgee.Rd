% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/conditional_functions.R
\name{quantile_cpmgee}
\alias{quantile_cpmgee}
\title{Calculate quantiles conditional on covariates}
\usage{
quantile_cpmgee(mod, cont_y, new.data, probs = 0.5, se = TRUE)
}
\arguments{
\item{mod}{the model}

\item{cont_y}{the continuous response data}

\item{new.data}{the new data}

\item{probs}{a numeric vector of pth quantiles}

\item{se}{if confidence intervals needed (default = TRUE)}
}
\value{
A list containing the following components:

\item{est}{a vector of estimated condtional quantiles}

\item{lb}{a vector of estimated lower bounds of 95\% confidence intervals}

\item{ub}{a vector of estimated upper bounds of 95\% confidence intervals}
}
\description{
This functions calculates the quantiles conditional on covariates based on the fitted model
and new data.
}
