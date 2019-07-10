# Gregory Way 2016 - GBM Immune Validation
# immune_survival.R
#
# Usage:
# Run in command line:
#
#       Rscript --vanilla immune_survival.R
#
# Data Required:
#
#       clinical data: "data/clinical_dataframe.tsv"
#       ssGSEA results: "results/ssGSEA_results.tsv"
#
# Output:
# A table showing a survival analysis based on immune cell infiltration and
# several Kaplan Meier survival curves

options(warn = -1)

suppressMessages(library(checkpoint))
suppressMessages(checkpoint("2016-08-16", checkpointLocation = "."))

SummariseCoxModel <- function(coxData, cell_type, sig_digits = 2) {
  # Function modified from: https://github.com/greenelab/hgsc_subtypes
  # 
  # This function will write out the hazards ratios, confidence intervals,
  # pvalues, and Wald's P for each cox proportional hazards model
  #
  # Args: 
  #   coxData: a coxph object
  #   cell_type: a string indicating which cell type is being investigated
  #   sig_digits: how many significant digits to display
  #
  # Returns:
  #   a dataframe of pertinent information regarding the cox model summary

  cox_summary <- summary(coxData)
  confidence_intervals <- round(cox_summary$conf.int[ ,3:4], sig_digits)
  hazard <- round(exp(coxData$coefficients), sig_digits)
  Pvalues <- format(cox_summary$coefficients[ ,5], scientific = TRUE,
                    digits = sig_digits)
  WaldsP <- format(cox_summary$waldtest["pvalue"], scientific = TRUE,
                   digits = sig_digits)
  CoxPH_results <- cbind(confidence_intervals, hazard, Pvalues, WaldsP)
  rownames(CoxPH_results)[1] <- cell_type
  return(CoxPH_results)
}

PerformSurvivalAnalysis <- function(cell_type, survival_data, output = TRUE,
                                    data_origin = 'TCGA') {
  # Performs an adjusted and unadjusted survival analysis according to the
  # estimated enrichment of immune cell types into GBM tumors
  #
  # Args: 
  #   cell_type: a string indicating the cell type to analyze
  #   survival_data: full dataframe with tumor specific information
  #   output: boolean to write files to disk (defaults to TRUE)
  #   data_origin: for purposes of file name - origin of the data is important
  #
  # Returns:
  #   a list of length two with unadjusted and adjusted models, will also write
  #   summaries of each model to the harddrive if `output = TRUE`

  surv_obj <- survival::Surv(as.numeric(paste(survival_data$days_to_death)),
                             event = survival_data$`_EVENT`)

  # Unadjusted model
  unadj_model <- survival::coxph(surv_obj ~ survival_data[[cell_type]])

  # Adjusted model by age, gender, and subtype
  adj_model <- surv_obj ~ survival_data[[cell_type]] + 
    survival_data$age_at_initial_pathologic_diagnosis + survival_data$gender +
    survival_data$GeneExp_Subtype

  adj_model <- survival::coxph(adj_model)

  unadj_model_summary <- SummariseCoxModel(unadj_model, cell_type)
  adj_model_summary <- SummariseCoxModel(adj_model, cell_type)

  if (output) {
    unadj_file <- file.path("results", paste0("unadjusted_", data_origin,
                                              "_model_", cell_type, ".tsv"))
    write.table(unadj_model_summary, unadj_file, sep = "\t", col.names = NA)
    adj_file <- file.path("results", paste0("adjusted_",data_origin,
                                            "_model_", cell_type, ".tsv"))
    write.table(adj_model_summary, adj_file, sep = "\t", col.names = NA)
  }

  return(list("unadjusted" = unadj_model_summary,
              "adjusted" = adj_model_summary))
}

PlotKaplanMeierCellType <- function(cell_type, surv_data, tertiles = FALSE) {
  # Fit a survival model and plot using survminer package
  #
  # Args: 
  #   cell_type: a string indicating the cell type to analyze
  #   surv_data: full dataframe with tumor specific information
  #   tertiles: boolean of splitting data into thirds (high, low, mid)
  #
  # Returns:
  #   Plots a survminer object to the screen

  if (tertiles) {
    tert <- quantile(surv_data[[cell_type]], c(.33, .66))
    cell_type_content <- rep("", nrow(surv_data))
    cell_type_content[surv_data[[cell_type]] > tert[2]] <- "High"
    cell_type_content[surv_data[[cell_type]] < tert[1]] <- "Low"
    cell_type_content[!(cell_type_content %in% c("Low", "High"))] <- "Mid"
    legend_labels <- c("High", "Mid", "Low")

  } else {
    cell_type_content <- rep("Low", nrow(surv_data))
    cell_type_content[surv_data[[cell_type]] > 
                        median(surv_data[[cell_type]])] <- "High"
    legend_labels <- c("High", "Low")
  }

  cell_type_content <- factor(cell_type_content, levels = legend_labels)
  surv_data$cell_type_content <- cell_type_content
  celltype_km <- survival::survfit(
    survival::Surv(as.numeric(paste(days_to_death)), `_EVENT`) ~
      cell_type_content, data = surv_data)
  
  print(celltype_km)
  survminer::ggsurvplot(celltype_km, legend = "right",
                        legend.title = cell_type,
                        color = legend_labels,
                        legend.labs = legend_labels,
                        conf.int = FALSE,
                        main = '',
                        xlab = "Days to Death",
                        size = 0.8)
}

# Load Constants
sur_plot_height <- 3
sur_plot_width <- 5

# Load and process data
clinical_file <- file.path("data", "clinical_dataframe.tsv")
ssgsea_file <- file.path("results", "ssGSEA_results.tsv")
ihc_file <- file.path("data", "validation_gbm_data_version3_clean.tsv")

clinical <- suppressMessages(readr::read_tsv(clinical_file))
ssgsea <-  suppressMessages(readr::read_tsv(ssgsea_file))
ihc <- suppressMessages(readr::read_tsv(ihc_file))

colnames(ssgsea)[1] <- "sampleID"

survival_data <- dplyr::full_join(clinical, ssgsea, by = "sampleID")
gene_exp_levels <- c("Classical", "Mesenchymal", "Neural", "Proneural")
survival_data$GeneExp_Subtype <- factor(survival_data$GeneExp_Subtype,
                                        levels = gene_exp_levels)
ihc$SUBTYPE <- factor(ihc$SUBTYPE, levels = gene_exp_levels)

# Rename IHC columns to match references to TCGA data
colnames(ihc) <- c("TMA_ID", "CD4", "CD8", "CD68", "CD163", "GeneExp_Subtype",
                   "age_at_initial_pathologic_diagnosis", "days_to_death",
                   "_EVENT", "gender")

# Perform analysis and write results to file
immune_cell_types <- c("CD4", "CD8", "CD68", "CD163")
cell_surv <- list()
ihc_surv <- list()
for (cell in immune_cell_types) {
  
  if (cell %in% c("CD68", "CD163")) {
    tcga <- "Macrophages"
  } else {
    tcga <- cell
  }
  
  cell_surv[[tcga]] <- PerformSurvivalAnalysis(cell_type = tcga,
                                               survival_data = survival_data,
                                               output = TRUE,
                                               data_origin = "TCGA")
  
  ihc_surv[[cell]] <- PerformSurvivalAnalysis(cell_type = cell,
                                              survival_data = ihc,
                                              output = TRUE,
                                              data_origin = "IHC")

  # Plot Kaplan-Meier curves for each cell type
  PlotKaplanMeierCellType(cell_type = tcga, surv_data = survival_data)
  plot_file <- file.path("figures", paste0("TCGA_kaplanmeier_", tcga, ".pdf"))
  ggplot2::ggsave(plot_file, width = sur_plot_width, height = sur_plot_height)
  
  PlotKaplanMeierCellType(cell_type = cell, surv_data = ihc)
  plot_file <- file.path("figures", paste0("IHC_kaplanmeier_", cell, ".pdf"))
  ggplot2::ggsave(plot_file, width = sur_plot_width, height = sur_plot_height)
  
  # Plot Kaplan-Meier curves for each cell type in tertiles
  PlotKaplanMeierCellType(cell_type = tcga, surv_data = survival_data,
                          tertiles = TRUE)
  plot_file <- file.path("figures", paste0("TCGA_kaplanmeier_tertiles_", tcga,
                                           ".pdf"))
  ggplot2::ggsave(plot_file, width = sur_plot_width, height = sur_plot_height)
  
  PlotKaplanMeierCellType(cell_type = cell, surv_data = ihc, tertiles = TRUE)
  plot_file <- file.path("figures", paste0("IHC_kaplanmeier_tertiles_", cell,
                                           ".pdf"))
  ggplot2::ggsave(plot_file, width = sur_plot_width, height = sur_plot_height)
  
  # Also save png files of each figure
  PlotKaplanMeierCellType(cell_type = tcga, surv_data = survival_data)
  plot_file <- file.path("figures", paste0("TCGA_kaplanmeier_", tcga, ".png"))
  ggplot2::ggsave(plot_file, width = sur_plot_width, height = sur_plot_height)
  
  # Plot Kaplan-Meier curves for each cell type in tertiles
  PlotKaplanMeierCellType(cell_type = tcga, surv_data = survival_data,
                          tertiles = TRUE)
  plot_file <- file.path("figures", paste0("TCGA_kaplanmeier_tertiles_", tcga,
                                           ".png"))
  ggplot2::ggsave(plot_file, width = sur_plot_width, height = sur_plot_height)
}

# Subtype specific Kaplan-Meier Curve
subtype_km <- survival::survfit(survival::Surv(as.numeric(paste(days_to_death)),
                                               `_EVENT`) ~ GeneExp_Subtype,
                                data = survival_data)
subtype_colors <- c("blue3", "darkorchid4", "green4", "orange1")

survminer::ggsurvplot(subtype_km, palette = subtype_colors,
                      legend = "right", legend.labs = gene_exp_levels,
                      legend.title = "Subtypes", conf.int = TRUE,
                      main = "TCGA GBM Survival",
                      xlab = "Days to Death")
plot_file <- file.path("figures", "TCGA_kaplanmeier_subtypes.pdf")
ggplot2::ggsave(plot_file, width = sur_plot_width, height = sur_plot_height)

# Save .png
survminer::ggsurvplot(subtype_km, palette = subtype_colors,
                      legend = "right", legend.labs = gene_exp_levels,
                      legend.title = "Subtypes", conf.int = TRUE,
                      main = "TCGA GBM Survival",
                      xlab = "Days to Death")
plot_file <- file.path("figures", "TCGA_kaplanmeier_subtypes.png")
ggplot2::ggsave(plot_file, width = sur_plot_width, height = sur_plot_height)