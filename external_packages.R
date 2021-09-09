check_Rcpp = require("Rcpp")

if (check_Rcpp == FALSE){ install.packages("Rcpp")}
check_RcppArmadillo = require("RcppArmadillo")
if (check_RcppArmadillo == FALSE){ install.packages("RcppArmadillo")}
check_combinat = require("combinat")
if (check_combinat == FALSE){ install.packages("combinat")}
check_RColorBrewer = require("RColorBrewer")
if (check_RColorBrewer == FALSE){ install.packages("RColorBrewer")}
check_gplots = require("gplots")
if (check_gplots == FALSE){ install.packages("gplots")}
