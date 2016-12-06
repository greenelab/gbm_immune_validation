# Gregory Way 2016 - GBM Immune Validation
# ssGSEA_lm22.R
#
# Runs single sample gene set enrichment analysis (ssGSEA) (Barbie et al. 2009)
# using LM22 defined immune cell gene sets (Newman et al. 2015).
#
# Usage:
# Run in command line:
#
#       Rscript --vanilla ssGSEA_lm22.R
#
# Output:
# single ssGSEA results for lm22 based genes

suppressMessages(library(checkpoint))
suppressMessages(checkpoint("2016-08-16", checkpointLocation = "."))

gene_exp <- readr::read_tsv(file.path("data", "HT_HG-U133A"))
gene_input <- as.matrix(gene_exp[, 1])
gene_exp <- as.matrix(gene_exp[, 2:ncol(gene_exp)])
rownames(gene_exp) <- gene_input

# This file is a processed subset of Sup Table S1_DEGs in Newman et al. 2015
# and is processed by `process_supplemental_data.R`
lm22_genes <- read.table(file.path("data", "ssGSEA_lm22_genes.tsv"),
                         sep = "\t", header = T, row.names = 1)
lm22_geneset <- lapply(lm22_genes, function(x) rownames(lm22_genes)[x == 1])

dir.create('results')

# Run ssgsea
ssgsea_result <- GSVA::gsva(expr = gene_exp, gset.idx.list = lm22_geneset,
                            method = "ssgsea", verbose = FALSE)
write.table(t(ssgsea_result), file.path("results", "ssGSEA_results.tsv"),
            sep = "\t", col.names = NA)
