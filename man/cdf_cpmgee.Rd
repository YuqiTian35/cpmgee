% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/conditional_functions.R
\name{cdf_cpmgee}
\alias{cdf_cpmgee}
\title{Calculate CDFs conditional on covariates}
\usage{
cdf_cpmgee(mod, cont_y, new.data, at.y = 0, se = TRUE)
}
\arguments{
\item{mod}{the model}

\item{cont_y}{the continuous response data}

\item{new.data}{the new data}

\item{at.y}{a numeric vector of cut-off points P(y <= at.y | new.data)}

\item{se}{if confidence intervals needed (default = TRUE)}
}
\value{
A list containing the following components:

\item{est}{a vector of estimated condtional CDFs}

\item{se}{a vector of estimated standard errors}

\item{lb}{a vector of estimated lower bounds of 95\% confidence intervals}

\item{ub}{a vector of estimated upper bounds of 95\% confidence intervals}
}
\description{
This functions calculates the CDFs conditional on covariates based on the fitted model
and new data.
}
