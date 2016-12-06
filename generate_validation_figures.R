# Gregory Way 2016 - GBM Immune Validation
# generate_validation_figures.R
#
# Usage:
# Run in command line:
#
#       Rscript --vanilla generate_validation_figures.R
#
# Output:
# Compare ssGSEA results to %positivity IHC estimates

options(warn = -1)

suppressMessages(library(checkpoint))
suppressMessages(checkpoint("2016-08-16", checkpointLocation = "."))

library(ggplot2)
library(dplyr)
library(gridExtra)
source("util/base_theme.R")

# For random jitter control in plots
set.seed(123)

clinical <- readr::read_tsv(file.path("data", "clinical_dataframe.tsv"))
validation <- readr::read_tsv(file.path("data", "validation_gbm_data.tsv"))
ssgsea <-  readr::read_tsv(file.path("results", "ssGSEA_results.tsv"))

# Prepare validation data
kept_columns <- c("SUBTYPE", "CD4 T cells", "CD8 T", "Macrophages")
validation_sum <- validation[, kept_columns] %>%
  group_by(SUBTYPE) %>%
  summarize_each(funs(mean(.)))

validation_sum <- reshape2::melt(validation_sum, id = "SUBTYPE")
validation_sum[, "study"] <- rep("validation", nrow(validation_sum))
colnames(validation_sum) <- c("Mixture", "subtype", "proportion", "study")
colnames(ssgsea)[1] <- "sampleID"

subtype_colors <- c("blue3", "darkorchid4", "green4", "orange1")

ssgsea_clin <-  dplyr::full_join(clinical, ssgsea, by = "sampleID")
keep_cols <- c("sampleID", "GeneExp_Subtype", "CD4", "CD8", "Macrophages")
ssgsea_subtype <- ssgsea_clin[, colnames(ssgsea_clin) %in% keep_cols]

ssgsea_subtype <- reshape2::melt(ssgsea_subtype,
                                 id.vars = c("sampleID", "GeneExp_Subtype"),
                                 measure.vars = c("CD4", "CD8", "Macrophages"),
                                 value.name = "Enrichment",
                                 variable.name = "Cell Type")

# Process validation data to compare with ssGSEA results
colnames(validation)[2:3] <- c("CD4", "CD8")
val_plot <- reshape2::melt(validation, id.vars = c("TMA_ID", "SUBTYPE"),
                           measure.vars = c("CD4", "CD8", "Macrophages"),
                           variable.name = "Cell Type",
                           value.name = "Positivity")
val_plot$`Cell Type` <- factor(val_plot$`Cell Type`,
                               levels = c("CD4", "CD8", "Macrophages"))

# Plot ssGSEA results
ssGSEA_theme <-   theme(axis.text.x = element_blank(),
                        axis.text.y = element_text(size = rel(0.6)),
                        axis.title.x = element_blank(),
                        axis.title.y = element_text(size = rel(0.9)),
                        axis.ticks = element_line(color = "black"),
                        axis.ticks.x = element_blank(),
                        axis.ticks.margin = unit(8, "mm"),
                        legend.position = "right",
                        plot.margin = unit(rep(0.1, 4), "cm"),
                        legend.text = element_text(size = rel(0.5)),
                        legend.key = element_blank(),
                        legend.key.size = unit(3, "mm"), 
                        strip.text.x = element_text(size = rel(0.5)))

ssgsea_grob <- ggplot2::ggplot(ssgsea_subtype, aes(x = GeneExp_Subtype,
                                                   y = Enrichment,
                                                   fill = GeneExp_Subtype)) +
  facet_wrap(~`Cell Type`, scales = "free") +
  geom_jitter(aes(color = GeneExp_Subtype), width = 0.2, size = rel(0.1)) +
  geom_boxplot(outlier.size = -1, lwd = 0.1) + 
  scale_fill_manual(values = subtype_colors) +
  scale_color_manual(values = subtype_colors) +
  ylab("ssGSEA Enrichment Score") + 
  theme_gbm() + ssGSEA_theme

validation_grob <- ggplot2::ggplot(val_plot, aes(x = SUBTYPE,
                                                 y = Positivity,
                                                 fill = SUBTYPE)) +
  facet_wrap(~`Cell Type`, scales = "free") +
  geom_jitter(aes(color = SUBTYPE), width = 0.2, size = rel(0.1)) +
  geom_boxplot(outlier.size = -1, lwd = 0.1) + 
  scale_fill_manual(values = subtype_colors) +
  scale_color_manual(values = subtype_colors) +
  ylab("Percent Positivity") + 
  theme_gbm() + ssGSEA_theme

# Extract the legend from the ssGSEA plots
gtable <- ggplot_gtable(ggplot_build(validation_grob))
legend_grob <- which(sapply(gtable$grobs, function(x) x$name == "guide-box"))
legend_grob <- gtable$grobs[[legend_grob]]

# Save Multiple Facet Boxplot
layout <- matrix(c(rep(1, 85), rep(2, 15)), nrow = 5, ncol = 20)
pdf(file.path("figures", "boxplot_validation_TCGA_summary.pdf"), width = 5,
    height = 3.3)
grid.arrange(arrangeGrob(ssgsea_grob + theme(legend.position = "none"),
                         validation_grob + theme(legend.position = "none")), 
             layout_matrix = layout,
             legend_grob, nrow = 1)
dev.off()

# Also save png for viewing in github
png(file.path("figures", "boxplot_validation_TCGA_summary.png"), width = 750,
    height = 505)
grid.arrange(arrangeGrob(ssgsea_grob + theme(legend.position = "none"),
                         validation_grob + theme(legend.position = "none")), 
             layout_matrix = layout,
             legend_grob, nrow = 1)
dev.off()
