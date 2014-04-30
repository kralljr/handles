\name{adjust}
\alias{adjust}
\title{Adjusting data below the MDL}
\usage{
adjust(x, ...)
}
\arguments{
  \item{data}{data frame of daily constituent
  concentrations with date as first column}

  \item{mdl}{vector or matrix of MDLs for each variable of
  data}

  \item{method}{one of three methods.  'substitute'
  substitutes all censored constant proportion of the mdl
  (sub * mdl).  'exclude' excludes constituents with more
  than experc censored data.  'likelihood' multiply imputes
  censored data with a likelihood-based method.}

  \item{sub}{proportion of MDL to substitute censored
  concentrations.  Default is 0.5.}

  \item{experc}{If a constituent exceeds this value, using
  the exclude method, the constituent will be dropped from
  the analysis.}

  \item{N}{number of draws from posterior for
  likelihood-based method}

  \item{burnin}{number of samples to discard for
  likelihood-based method #param ... other arguments}
}
\description{
\code{adjust} adjusts censored concentrations below the MDL
}
\details{
These are functions to create complete PM2.5 constitutent
concentrations by substituting, excluding, or imputing
censored data. Works on dataframe where the first column is
date and all subsequent columns are concentrations of
chemical constituents.
}
\examples{
adjust(nycdat, )
apca(nycdat, mdl = mdls, adjust = "substitute")
}
\seealso{
Other adjust1: \code{\link{print.adjust}}
}

