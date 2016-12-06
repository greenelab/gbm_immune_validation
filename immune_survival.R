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
    survival_data$age_recode + survival_data$gender +
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

PlotKaplanMeierCellType <- function(cell_type, survival_data,
                                    tertiles = FALSE) {
  # Fit a survival model and plot using survminer package
  #
  # Args: 
  #   cell_type: a string indicating the cell type to analyze
  #   survival_data: full dataframe with tumor specific information
  #   tertiles: boolean of splitting data into thirds (high, low, mid)
  #
  # Returns:
  #   Plots a survminer object to the screen

  if (tertiles) {
    tert <- quantile(survival_data[[cell_type]], c(.33, .66))
    cell_type_content <- rep("", nrow(survival_data))
    cell_type_content[survival_data[[cell_type]] > tert[2]] <- "High"
    cell_type_content[survival_data[[cell_type]] < tert[1]] <- "Low"
    cell_type_content[!(cell_type_content %in% c("Low", "High"))] <- "Mid"
    legend_labels <- c("High", "Mid", "Low")
  } else {
    cell_type_content <- rep("Low", nrow(survival_data))
    cell_type_content[survival_data[[cell_type]] > 
                        median(survival_data[[cell_type]])] <- "High"
    legend_labels <- c("Low", "High")
  }

  cell_type_content <- factor(cell_type_content, levels = legend_labels)
  survival_data$cell_type_content <- cell_type_content
  celltype_km <- survival::survfit(survival::Surv(as.numeric(paste(days_to_death)),
                                                `_EVENT`) ~ cell_type_content,
                                 data = survival_data)
  
  print(celltype_km)
  survminer::ggsurvplot(celltype_km, legend = "right",
                        legend.title = cell_type,
                        legend.labs = legend_labels,
                        conf.int = TRUE,
                        main = '',
                        xlab = "Days to Death",
                        size = 0.3)
}

# Load Constants
surv_plot_height <- 3
surv_plot_width <- 5

# Load and process data
clinical_file <- file.path("data", "clinical_dataframe.tsv")
ssgsea_file <- file.path("results", "ssGSEA_results.tsv")
clinical <- suppressMessages(readr::read_tsv(clinical_file))
ssgsea <-  suppressMessages(readr::read_tsv(ssgsea_file))
colnames(ssgsea)[1] <- "sampleID"

survival_data <- dplyr::full_join(clinical, ssgsea, by = "sampleID")
gene_exp_levels <- c("Classical", "Mesenchymal", "Neural", "Proneural")
survival_data$GeneExp_Subtype <- factor(survival_data$GeneExp_Subtype,
                                        levels = gene_exp_levels)

# Recode age into discretized bins
survival_data$age_recode <- 
  car::recode(survival_data$age_at_initial_pathologic_diagnosis,
                "0:14 = 1; 15:39 = 2; 40:44 = 3; 45:49 = 4; 50:54 = 5;
                55:59 = 6; 60:64 = 7; 65:69 = 8; 70:74 = 9; 75:150 = 10")

# Perform analysis and write results to file
immune_cell_types <- c("CD4", "CD8", "Macrophages")
cell_surv <- list()
for (cell in immune_cell_types) {
  cell_surv[[cell]] <- PerformSurvivalAnalysis(cell_type = cell,
                                               survival_data = survival_data,
                                               output = TRUE,
                                               data_origin = "TCGA")

  # Plot Kaplan-Meier curves for each cell type
  PlotKaplanMeierCellType(cell, survival_data)
  plot_file <- file.path("figures", paste0("TCGA_kaplanmeier_", cell, ".pdf"))
  ggplot2::ggsave(plot_file, width = surv_plot_width, height = surv_plot_height)
  
  # Plot Kaplan-Meier curves for each cell type in tertiles
  PlotKaplanMeierCellType(cell, survival_data, tertiles = TRUE)
  plot_file <- file.path("figures", paste0("TCGA_kaplanmeier_tertiles_", cell,
                                           ".pdf"))
  ggplot2::ggsave(plot_file, width = surv_plot_width, height = surv_plot_height)
  
  # Also save png files of each figure
  PlotKaplanMeierCellType(cell, survival_data)
  plot_file <- file.path("figures", paste0("TCGA_kaplanmeier_", cell, ".png"))
  ggplot2::ggsave(plot_file, width = surv_plot_width, height = surv_plot_height)
  
  # Plot Kaplan-Meier curves for each cell type in tertiles
  PlotKaplanMeierCellType(cell, survival_data, tertiles = TRUE)
  plot_file <- file.path("figures", paste0("TCGA_kaplanmeier_tertiles_", cell,
                                           ".png"))
  ggplot2::ggsave(plot_file, width = surv_plot_width, height = surv_plot_height)
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
ggplot2::ggsave(plot_file, width = surv_plot_width, height = surv_plot_height)

# Save .png
survminer::ggsurvplot(subtype_km, palette = subtype_colors,
                      legend = "right", legend.labs = gene_exp_levels,
                      legend.title = "Subtypes", conf.int = TRUE,
                      main = "TCGA GBM Survival",
                      xlab = "Days to Death")
plot_file <- file.path("figures", "TCGA_kaplanmeier_subtypes.png")
ggplot2::ggsave(plot_file, width = surv_plot_width, height = surv_plot_height)