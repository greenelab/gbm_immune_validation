# Gregory Way 2016 - GBM Immune Profiles
# util/base_theme.R
#
# Usage:
# Sourced only
#
# Serves as the base themes used in ggplot2 outputs

require(ggplot2)

base_density <- geom_density(aes(group = GeneExp_Subtype,
                                 fill = GeneExp_Subtype),
                             alpha = 0.3, size = 0.1, na.rm = TRUE)

theme_gbm <- function(base_size = 12, base_family = "") {
  ggplot2::theme_bw(base_size = base_size, base_family = base_family) %+replace%
    ggplot2::theme(title = element_text(size = rel(0.8)),
                   axis.title = element_text(size = rel(0.9)),
                  axis.text.x = element_text(size = rel(0.7), angle = 45),
                  axis.text.y = element_text(size = rel(0.7)),
                  axis.line.x = element_line(color = "black", size = 0.2),
                  axis.line.y = element_line(color = "black", size = 0.2),
                  axis.ticks = element_blank(),
                  legend.title = element_blank(),
                  legend.text = element_text(size = rel(0.4)),
                  legend.key.size = unit(0.3, "cm"),
                  legend.position = c(0.85, 0.7),
                  panel.grid.major = element_line(color = "white", size = 0.3),
                  panel.grid.minor = element_line(color = "white", size = 0.3),
                  panel.background = element_rect(fill = "white"))
}