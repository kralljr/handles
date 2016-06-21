#' HANdles Detection Limits when Estimating Sources
#'
#' handles helps users estimate sources of PM2.5 with censored data below
#' minimum detection limits (MDL) by first adjusting censored PM2.5 
#' constituent concentrations and then performing source apportionment
#' using Absolute Principal Component Analysis (APCA). 
#' Methods for adjusting censored data include a constant proportion
#' of the MDL (e.g. 1/2 MDL), removing constituents with many 
#' censored observations, and a likelihood-based method. 
#' 
#' Source apportionment
#' 
#' To perform source apportionment, users should have their 
#' data as a dataframe where the first column is the date and
#' other columns correspond to constituent concentrations
#' \code{\link{apca}}
#'
#' Adjusting censored data below MDLs
#'
#' To adjust censored data, users can select a substitution method,
#' which substitutes censored data with a constant proportion of the MDL,
#' an exclude method, which excludes constituents with many censored
#' concentrations, and a likelihood-based multiple imputation approach
#' \code{\link{adjust}}
#'
#' @references George D. Thurston, John D. Spengler (1985).  A quantitative 
#' assessment of source contributions to inhalable particulate matter pollution
#' in metropolitan Boston, 19(1) 9-25. 
#' @docType package
#' @name handles
#' @import nlme
#' @import class
#' @import MCMCpack
#' @import truncnorm
NULL
