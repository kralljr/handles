% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sim_APCA.R
\name{rmdups}
\alias{rmdups}
\title{\code{rmdups} reassigns classification for duplicate classifications}
\usage{
rmdups(profile, class, probs, k, sprofiles = speciate)
}
\arguments{
\item{profile}{Source apportionment profile matrix with number of columns equal to number of chemical constituents and number of rows equal to number of sources}

\item{class}{character vector of classifications}

\item{probs}{vector of probabilities corresponding to classifications}
}
\description{
\code{rmdups} reassigns classification for duplicate classifications
}
