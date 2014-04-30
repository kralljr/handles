\name{apca}
\alias{apca}
\title{Absolute Principal Component Analysis (APCA)}
\usage{
apca(x, ...)
}
\arguments{
  \item{data}{data frame of daily constituent
  concentrations with date as first column}

  \item{tots}{vector of total concentrations (total PM2.5)
  for each day.  If null, uses \code{rowSums(data)}}

  \item{nsources}{number of sources.  If null, uses number
  of eigenvalues of the correlation matrix greater than
  one.}

  \item{adjust}{method to adjust censored concentrations
  below minimum detection limits (MDLs).  If null, uses
  complete case data}

  \item{mdl}{either a vector of mdls corresponding to each
  constituent or a matrix of mdls corresponding to each
  constituent and each day. #param ... other arguments}
}
\description{
\code{apca} performs APCA on daily PM2.5 constituent
concentrations
}
\details{
This is a function to estimate PM2.5 source profiles and
daily PM2.5 source concentrations from PM2.5 constituent
concentrations observed at ambient monitors (e.g. EPA
Chemical Speciation Network). Works on dataframe where the
first column is date and all subsequent columns are
concentrations of chemical constituents. See Thurston and
Spengler (1985, Atmospheric Environment)
}
\examples{
data(nycdat)
data(nycmdl)
apca(nycdat)
apca(nycdat, mdl = nycmdl, adjust = "substitute")
}

