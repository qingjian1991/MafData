#create a new
#http://tinyheero.github.io/jekyll/update/2015/07/26/making-your-first-R-package.html


#Setup
#install.packages("devtools")
#install.packages("roxygen2")

#1) Creating the Framework for your First Package

devtools::create("MafData")

##2) add dependence
devtools::use_package("dplyr") # Defaults to imports
#> Adding dplyr to Imports
#> Refer to functions with dplyr::fun()
devtools::use_package("dplyr", "Suggests")
devtools::use_package("tidyverse")
#> Adding dplyr to Suggests
#> Use requireNamespace("dplyr", quietly = TRUE) to test if package is
#>  installed, then use dplyr::fun() to refer to functions.

#3) add addiction Data

#R data
x <- c(1:10)
devtools::use_data(x)

# raw data
devtools::use_data_raw()
devtools::use_data()


#6) Documnet
usethis::use_vignette("MafProcess")

#5) How do I Document My Functions?
devtools::document()

#4) Edit your code and load your code
devtools::load_all()

#devtools::use_gpl3_license()

# 7) install your packages
devtools::install()

devtools::install(build_vignettes = T)




