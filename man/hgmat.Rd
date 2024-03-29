% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hgmat.R
\name{hgmat}
\alias{hgmat}
\title{Calculate H and G matrix}
\usage{
hgmat(mod, smat, X, modtype, diffmeth, alpha, corrmod, h)
}
\arguments{
\item{mod}{the fitted model with independence working correlation structure}

\item{smat}{the matrix of correlations among elements in $Z_{it}$}

\item{X}{the design matrix}

\item{modtype}{a character string specifying the model type (CPM or GLM)}

\item{diffmeth}{a character string specifying the method used for estimation of alpha}

\item{alpha}{the association parameter}

\item{corrmod}{the specified working correlation structure}

\item{h}{a numeric value of options to control the fitting algorithm}
}
\value{
A list including the working correlation matrix, H, G, and derivatives of H and G
}
\description{
See details in Parsons N., et al. "A generalized estimating 
equation method for fitting autocorrelated ordinal score data with an application in 
horticultural research." Journal of the Royal Statistical Society: Series C 
(Applied Statistics) 55.4 (2006): 507-524.
}
