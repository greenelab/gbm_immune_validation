# Gregory Way 2016 - GBM Immune Profiles
# describe_data.R
#
# Usage:
# Run in command line:
#
#      Rscript --vanilla describe_data.R
#
# Output:
# 1) Subtype distribution .png plots for age, survival, and gender.
# 2) Clinical dataframe used in downstream analyses

options(warn = -1)

suppressMessages(library(checkpoint))
suppressMessages(checkpoint("2016-08-16", checkpointLocation = "."))

library(ggplot2)
library(dplyr)
source(file.path("util", "base_theme.R"))  # Default plotting parameters

# 1. Load Data
# Load full clinical matrix
clinical_df <- readr::read_tsv(file.path("data", "GBM_clinicalMatrix"))

# Load gene expression matrix column names
gene_exp <- readr::read_lines(file.path("data", "HT_HG-U133A"), n_max = 1)
gene_exp <- unlist(strsplit(gene_exp, "\t"))
gene_exp <- gene_exp[2:length(gene_exp)]

# 2. Subset clinical matrix to only samples with gene expression measurements
clin_gene_df <- clinical_df %>% dplyr::filter(sampleID %in% gene_exp)

# How many samples does this remove?
# nrow(clinical_df) - nrow(clin_gene_df) # 90

# How many samples are primary tumor (01), duplicate (02), normal (11)?

# table(clin_gene_df$sample_type, useNA = "ifany")
# Primary Tumor Solid Tissue Normal 
#           529                  10 

# What is the distribution of subtypes?
# table(clin_gene_df$GeneExp_Subtype, useNA = "ifany")
#Classical Mesenchymal      Neural   Proneural        <NA> 
#      145         158          87         139          10 

# For plots, filter only to primary tumors with subtype classification
clin_subtype_df <- clin_gene_df %>%
  dplyr::filter(sampleID %in% gene_exp & sample_type == "Primary Tumor")

# 3. Plot subtype distributions
dir.create("figures", showWarnings = FALSE)

# What is the distribution of age across subtypes?
ggplot(clin_subtype_df, aes(x = age_at_initial_pathologic_diagnosis)) +
  base_density + xlim(c(0, 110)) + xlab("Age at Initial Diagnosis") +
  ylab("Density") + theme_gbm()
ggsave(file.path("figures", "age_distribution.png"), units = "in",
                 width = 3, height = 2)

# What is the distribution of survival times across subtypes?
ggplot(clin_subtype_df, aes(x = as.numeric(paste(days_to_death)))) +
  base_density + xlab("Days to Death") + ylab("Density") + xlim(c(0, 5000)) +
  theme_gbm()
ggsave(file.path("figures", "survival_distribution.png"), units = "in",
       width = 3, height = 2)

# What is the distribution of genders across subtypes?
ggplot(clin_subtype_df, aes(x = gender)) +
  geom_bar(aes(group = GeneExp_Subtype, fill = GeneExp_Subtype), na.rm = TRUE) +
  xlab("Gender") + ylab("Count") + theme_gbm() +
  theme(axis.text.x = element_text(size = rel(0.7), angle = 0),
        legend.position = "right")
ggsave(file.path("figures", "gender_distribution.png"), units = "in",
       width = 3, height = 2)

# 4. Output clinical dataframe to use in downstream analysis
clinical_columns <- c("sampleID", "sample_type", "GeneExp_Subtype",
                      "age_at_initial_pathologic_diagnosis", "_EVENT", "gender",
                      "days_to_death")
clinical_output <- clin_gene_df %>% dplyr::select(one_of(clinical_columns))
readr::write_tsv(clinical_output, file.path("data", "clinical_dataframe.tsv"))
