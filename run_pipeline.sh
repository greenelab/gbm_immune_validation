#!/bin/bash

# Install R Dependencies
Rscript --vanilla install.R

# Download Gene Expression and Clinical Information
./download_gbm.sh

# Describe data and output subset clinical matrix
Rscript --vanilla describe_data.R

# Run ssGSEA on TCGA GBM data using LM22 gene sets
Rscript --vanilla ssGSEA_lm22.R

# Output final results figures compared to IHC validation
Rscript --vanilla generate_validation_figures.R

# Perform survival analysis using ssGSEA enrichment scores
Rscript --vanilla immune_survival.R

