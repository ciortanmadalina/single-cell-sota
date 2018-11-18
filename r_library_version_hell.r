# This will upgrade R lang and break other libraries
#source("http://cf.10xgenomics.com/supp/cell-exp/rkit-install-2.0.0.R")


install.packages("dplyr")
require(devtools)
devtools::install_url("http://cran.r-project.org/src/contrib/Archive/dplyr/dplyr_0.1.2.tar.gz")
# we need rlang version 0.3.0 not higher, otherwise outlier removal function won't work
packageVersion('rlang')
devtools::install_version("rlang", version = "0.3.0", dep = FALSE)

install.packages("clues")

devtools::install_github("VCCRI/CIDR")

install.packages("plot3D")
install.packages("rgl")

