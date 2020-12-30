# This R scrpt details the step required to build the CRAN .tar.gz file for submission to CRAN and how to build the manual PDF.
#------------------------------------------------

# Build PDF. If AWAPer.pdf already exists, then delete before running.
path <- find.package("AWAPer")
system(paste(shQuote(file.path(R.home("bin"), "R")),"CMD", "Rd2pdf", shQuote(path)))

# Pre-build the vignettes to avod their execuation during package installation.
# NOTE: Must manually move image files from vignettes/figure/ to vignettes/ after knit
#       and remove figure/ from .Rmd
library(knitr)
knitr::knit("vignettes/Point_rainfall.Rmd.orig", output = "vignettes/Point_rainfall.Rmd")
knitr::knit("vignettes/Catchment_avg_ET_rainfall.Rmd.orig", output = "vignettes/Catchment_avg_ET_rainfall.Rmd")
browseVignettes("AWAPer")

# Build the pavkage for CRAN
devtools::build(".")
