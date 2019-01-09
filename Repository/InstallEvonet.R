if (!require("devtools")) {
  install.packages("devtools")
  library(devtools)}
devtools::install_github( "statnet/tergmLite",force=T)
devtools::install_github( "statnet/EpiModel", ref ="fast_edgelist",force=T)
install_github("EvoNetHIV/EvoNet",subdir="pkg",ref="absdiffby",force=TRUE)
library(evonet)