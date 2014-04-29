handles
=======

HANdles Detection Limits when Estimating Sources


Functions to adjust censored PM2.5 constituent concentrations that fall below minimum detection limits (MDL).  With handles you can:

  1. Adjust censored data using 3 methods:
    * Substitute censored concentrations with a constant proportion of the MDL (e.g. 1/2 MDL)
    * Exclude constituents with many missing observations
    * Multiply impute the censored data using a likelihood-based method
  2. Apply Absolute Principal Component Analysis (APCA, Thurston and Spengler (1985) *Atmospheric Environment*)
  
Many researchers use reported data from ambient PM2.5 constituent monitors directly in source apportionment analyses.  The reported data frequently fall below MDLs and these reported values may have no association with the true concentrations below the MDLs.  handles allows users to compare source apportionment results between several methods for adjusting censored data. handles implements two methods common in source apportionment literature: Substitute and Exclude.  Additionally, handles implements a likelihood-based method, which is preferred for imputing censored multivariate data (see Helsel (2010) *Occupational Hygiene*, Hopke et al. (2001) *Biometrics*).  
