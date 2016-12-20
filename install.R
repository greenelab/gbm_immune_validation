# Gregory Way 2016 - GBM Immune Validation
# install.R
#
# The script installs all R dependencies required for the pipeline
#
# Usage:
# Run once before implementing pipeline: "Rscript install.R"

library("methods")

mirror <- "http://cran.us.r-project.org"
install.packages("checkpoint", repos = mirror)

library("checkpoint")
dir.create(".checkpoint")
checkpoint("2016-08-16", checkpointLocation = ".")

library("optparse")
library("readr")
library("ggplot2")
library("dplyr")
library("data.table")
library("survival")
library("survminer")
library("car")
library("Rcpp")

# Load bioconductor packages
source("https://bioconductor.org/biocLite.R")
biocLite("preprocessCore", suppressUpdates = TRUE)
biocLite("GSVA", suppressUpdates = TRUE)

library("preprocessCore")
library("GSVA")

sink("sessionInfo.txt")
sessionInfo()
sink()
