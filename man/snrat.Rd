% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/adjust.R
\name{snrat}
\alias{snrat}
\title{snrat}
\usage{
snrat(dat, mdl, unc)
}
\arguments{
\item{mdl}{matrix of MDLs with dimension the same as the data}

\item{unc}{matrix of uncertainties with dimension the same as the data}

\item{data}{data frame of daily constituent concentrations (ndays X ncons)}
}
\description{
\code{snrat} adjusts censored concentrations below the MDL using signal
to noise ratio approach for inflating uncertainty in PMF
}
