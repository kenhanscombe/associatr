
#' A "manhattan" plot for genome-wide genetic association results
#'
#' A scatter plot of negative log base 10 p-value against genomic position
#'
#' @param data A genetic association results dataframe.
#' @param trait_label A character vector to be used as label. Default = "phenotype".
#' @param trim Observations with p-values larger than the trim value are dropped. Default = 0.001.
#' @param chromosome An integer vector of chromosome codes.
#' @param snp SNP id.
#' @param bp Base-pair coordinate.
#' @param a1 Allele 1 (usually minor). See details.
#' @param effect Regression coefficient or odds ratio.
#' @param p Association p-value.
#' @param chr_color A character vector of 2 colours to be used on alternate chromosomes.
#' @param strip_color Color of the vertical strip describing the panel.
#' @param ... Additional arguments passed to \code{geom_point}.
#'
#' @details SNP id (\code{snp}), Allele 1 (\code{a1}), and regression coefficient or odds ratio (\code{effect}) columns are only used for labelling points. \code{a1} is described as it is in \href{https://www.cog-genomics.org/plink/1.9/formats#assoc_linear}{PLINK} documentation, but as this column is only used for labelling it can refer to any allele. PLINK output names are used as the default values for all columns.
#'
#' @import dplyr ggplot2
#' @importFrom ggrepel geom_text_repel
#' @importFrom magrittr "%>%"
#'
#' @export
gwa_manhattan <- function(data, trait_label = "phenotype", trim = 0.001,
                          chromosome = "#CHROM", snp = "SNP", bp = "POS",
                          a1 = "A1", effect = "BETA", p = "P",
                          chr_color = c("#878D92", "#49494D"),
                          strip_color = "grey25", ...) {
  df <- data[data[[p]] <= trim, ] %>%
    mutate(panel = trait_label)

  # Genomic position
  max_bp <- tapply(df[[bp]], df[[chromosome]], function(x) max(x))
  cum_max_bp <- cumsum(as.numeric(max_bp))
  df$POS <- df[[bp]] + c(0, cum_max_bp)[df[[chromosome]]]
  chr_midpoints <- tapply(df$POS, df[[chromosome]], function(x) mean(range(x)))

  # Colour scheme
  chr_col <- rep_len(chr_color, 22)
  names(chr_col) <- levels(factor(df[[chromosome]]))

  min_p_by_chr <- tapply(df[[p]], as.factor(df[[chromosome]]), min,
                         na.rm = TRUE)
  max_pos <- max(df$POS, na.rm = TRUE)

  # Scatterplot
  df %>%
    ggplot(aes(POS, -log10(df[[p]]))) +
    geom_hline(yintercept = -log10(5e-08), linetype = 3, col = "orange",
               size = 0.15) +
    annotate("text", x = max_pos, y = -log10(5e-08), vjust = 1,
             hjust = 1, label = "italic('p = 5e-08')", parse = TRUE,
             size = 2.5, col = "orange") +
    geom_point(aes(color = factor(df[[chromosome]])), na.rm = TRUE, size = 1,
               shape = 19, ...) +
    scale_colour_manual(values = chr_col) +
    guides(color = FALSE) +
    scale_x_continuous(name = "Chromosome", breaks = chr_midpoints,
                       labels = as.character(1:22)) +
    theme(panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_text(color = chr_color[1], size = 6, vjust = 2.5),
          strip.text.y = element_text(face = "plain", colour = "black",
                                      family = "Helvetica"),
          strip.background = element_rect(fill = strip_color)) +
    labs(y = expression(-log[10] (p))) +
    facet_grid(panel ~ ., scales = "free")
}
# element_rect(color = NULL, fill = grey(.975))

