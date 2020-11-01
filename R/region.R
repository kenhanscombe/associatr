globalVariables(c("#CHROM", ".", "P", "POS", "R2", "chr", "chromEnd",
                  "chromStart", "combined_rate", "end", "genetic_map_b37",
                  "glist_hg19", "grey", "name", "plot_loc", "position",
                  "start", "x_end", "x_start"))

#' Regional plot
#'
#' Plots is a negative log base 10 p-value against genomic position
#' scatterplot for a relatively small genomic "window". Below this
#' \bold{Association} layer, there are additional layers or "tracks"
#' providing biological context: hg19 build 37 combined
#' \bold{Recombination Rate}, known \bold{refSeq Genes}, and optional
#' tracks including genetic regulation information.
#'
#' @param data A genetic association results dataframe.
#' @param recomb_map A recombination map dataframe with header and
#' obligatory columns `position` (genomic position), `combined_rate`
#' (combined recombination rate cM/Mb), and `chr` (chromosome). Use
#' \code{\link{get_map}} to retrieve a build 37 map for any 1000
#' genomes phase 3 population.
#' @param target A character vector for the variant of interest.
#' @param chromosome An integer indicating the chromosome on which the
#' region sits.
#' @param target_bp A integer indicating base pair position in the
#' centre of the region.
#' @param pop A string argument passed to LDlinkR. A 1000 Genomes
#' population codes, or super population codes can be used. Multiple
#' codes allowed. Default is super population "EUR" (includes "CEU",
#' "TSI", "FIN", "GBR", "IBS").
#' @param window_kb A integer indicating window size in kb (default is
#' 150kb), to be included either side of the \code{target_bp} location.
#' Maximum is 500kb.
#' @param regulation A boolean (default \code{FALSE}) indicating
#' whether or not to include regulation tracks.
#' @param token An LDlinkR 'Personal Access Token'. Default is `NULL`.
#' See Details.
#'
#' @details You will need a genetic map. Use \code{\link{get_map}} to download 1000 Genomes phase 3 recombination map data. To retrieve LD data for the region of interest, you will need an LDlinkR 'Personal Access Token' described \href{https://cran.r-project.org/web/packages/LDlinkR/vignettes/LDlinkR_vignette_v8-2_TM.html}{here}. Apply for a token \href{https://ldlink.nci.nih.gov/?tab=apiaccess}{here}.
#'
#' @import dplyr tidyr purrr grid ggplot2 ggrepel RMySQL stringr DBI
#' @importFrom magrittr "%>%"
#' @importFrom LDlinkR LDproxy
#' @export
region <- function(data, recomb_map, target, chromosome, target_bp,
                   pop = "EUR", window_kb = 150, regulation = FALSE,
                   token = NULL) {
    if (is.null(token)) {
        stop("Argument `token` requires an LDlinkR 'Personal Access Token' described here https://cran.r-project.org/web/packages/LDlinkR/vignettes/LDlinkR_vignette_v8-2_TM.html. Apply for a token here https://ldlink.nci.nih.gov/?tab=apiaccess", call. = FALSE)
    }

    if (window_kb > 500) {
        stop("Maximum window is +/-500kb.")
    }

    window <- window_kb * 1e3
    window_min <- target_bp - window
    window_max <- target_bp + window
    plot_margin <- margin(0, 5.5, 0, 5.5)

    df <- dplyr::filter(data, ((`#CHROM` == chromosome) &
                                (POS %in% window_min:window_max)))
    df <- tidyr::unite(df, col = "variant", c("#CHROM", "POS"), sep = ":",
                       remove = FALSE)

    region_p_min <- df %>%
        dplyr::select(P) %>%
        min(.$P, na.rm = TRUE)

    track_list <- track_query(chromosome, window_max, window_min, regulation)
    ld <- ld_query(target, pop, token)
    df <- df %>% dplyr::left_join(ld, by = "variant")

    p_assoc <- plot_assoc(df, target_bp, chromosome, window_max, window_min,
                          window_kb, region_p_min)
    p_recom <- plot_recom(df, recomb_map, target_bp, chromosome, window_max,
                          window_min, window_kb, region_p_min, plot_margin)
    p_gene <- plot_gene(refGene = track_list[[1]], target_bp, chromosome,
                        window_max, window_min, window_kb, region_p_min,
                        plot_margin)

    p_list <- if (regulation) {
        p_cpg <- plot_cpg(cpgIslandExt = track_list[[2]], target_bp, chromosome,
                          window_max, window_min, window_kb, region_p_min,
                          plot_margin)

        p_tfbs <- plot_tfbs(wgEncodeRegTfbsClusteredV2 = track_list[[3]], target_bp,
                            chromosome, window_max, window_min, window_kb,
                            region_p_min, plot_margin)

        list(p_assoc, p_recom, p_gene, p_cpg, p_tfbs)
    } else {
        list(p_assoc, p_recom, p_gene)
    }

    g <- purrr::map(p_list, ggplotGrob) %>%
        purrr::reduce(rbind, size = "first")

    grid::grid.newpage()
    grid::grid.draw(g)
}




# Retrieves LD with target variant
ld_query <- function(target, pop, token) {
    LDlinkR::LDproxy(snp = target,
                     pop = pop,
                     r2d = "r2",
                     token = token) %>%
    dplyr::mutate(variant = stringr::str_replace(Coord, "chr", "")) %>%
    tidyr::separate(
        col = "variant",
        into = c("chr", "bp"),
        sep = ":",
        remove = FALSE,
        convert = TRUE
    )
}


# SQL query for requested tracks
track_query <- function(chromosome, window_max, window_min, regulation = FALSE) {
    if(regulation) {
        con_ucsc <- DBI::dbConnect(RMySQL::MySQL(), db = "hg19", user = "genome",
                            host = "genome-mysql.soe.ucsc.edu")

        cpgIslandExt <- suppressWarnings(DBI::dbGetQuery(
        con_ucsc,
        stringr::str_interp(
            "SELECT chrom, name, chromStart, chromEnd
            FROM cpgIslandExt
            WHERE chrom = 'chr${chromosome}' AND chromStart <= ${window_max} AND chromEnd >= ${window_min}")))

        wgEncodeRegTfbsClusteredV2 <- suppressWarnings(DBI::dbGetQuery(
        con_ucsc,
        stringr::str_interp(
            "SELECT chrom, chromStart, chromEnd, thickStart, thickEnd, name
            FROM wgEncodeRegTfbsClusteredV2
            WHERE chrom = 'chr${chromosome}' AND chromStart <= ${window_max} AND chromEnd >= ${window_min}")))

        DBI::dbDisconnect(con_ucsc)
        rm(con_ucsc)
        return(list(cpgIslandExt, wgEncodeRegTfbsClusteredV2))
    }
}


# Formats scale of BP axis to MB
pos_mb <- function(x) {format(round(x / 1e6, 1), nsmall = 1)}




# Subplots of the regional plot -------------------------------------------

# SNP association
plot_assoc <- function(data, target_bp, chromosome, window_max, window_min,
                       window_kb, region_p_min) {
    data <- data %>%
        dplyr::filter(`#CHROM` == chromosome & POS %in% window_min:window_max) %>%
        dplyr::mutate(panel = "Association")

    ggplot(data, aes(POS, -log10(P), color = R2)) +
        geom_point(data = filter(data, is.na(R2)), size = 1, color = grey(.9)) +
        geom_point(data = filter(data, !is.na(R2)), size = 1) +
        facet_grid(panel ~ ., scales = "free") +
        annotate("line", x = c(window_min, window_max), y = -log10(5e-8),
                 color = "grey", linetype = 3, size = 0.25) +
        annotate("text", x = window_max, y = -log10(5e-8) + 0.5, color = "grey",
                 label = "italic('p = 5e-08')", parse = TRUE, size = 2.5,
                 hjust = 1) +
        annotate("line", x = target_bp, y = c(0, -log10(region_p_min)),
                 color = "red", linetype = 3, size = 0.25) +
        annotate("point", x = target_bp, y = -log10(region_p_min), shape = 1,
                 size = 2, color = "red") +
        annotate("text", x = target_bp, y = -0.6, color = "grey35",
                 label = stringr::str_interp(
                     "+/- ${window_kb} Kb window around target SNP"), size = 3) +
        geom_rug(sides = "b", color = "grey", alpha = 0.25) +
        scale_color_gradient2(low = "grey", mid = "orange", high = "red",
                            midpoint = 0.5, na.value = "grey") +
        scale_x_continuous(labels = pos_mb, limits = c(window_min, window_max)) +
        scale_y_continuous(limits = c(-1.0, max(-log10(5e-08) + 2,
                                                -log10(region_p_min) + 2))) +
        labs(x = stringr::str_interp("Position on chromosome ${chromosome} (Mb)"),
            y = expression(-log[10] (p)), color = "LD") +
        guides(color = guide_colorbar(barwidth = 10, barheight = 0.5)) +
        theme(
            legend.position = "top",
            legend.title = element_text(face = "bold", vjust = 0.25),
            panel.background = element_blank(),
            panel.grid = element_blank(),
            axis.line = element_blank(),
            strip.text.y = element_text(face = "bold", colour = "white"),
            strip.background = element_rect(fill = "grey35"),
            plot.margin = margin(0, 5.5, 3, 5.5)
        )
}


# 1000G b37 combined recombination rate
plot_recom <- function(df, recomb_map, target_bp, chromosome, window_max,
                       window_min, window_kb, region_p_min, plot_margin) {
    p_recom <- recomb_map %>%
        dplyr::filter(chr == chromosome & position %in% window_min:window_max) %>%
        dplyr::mutate(panel = "Recombination") %>%
        ggplot(aes(position, combined_rate, group = chr)) +
        facet_grid(panel ~ ., scales = "free") +
        geom_vline(xintercept = target_bp, color = "red", linetype = 3, size = 0.25) +
        geom_path(color = "cornflowerblue") +
        ylim(0, 100) +
        scale_x_continuous(limits = c(window_min, window_max)) +
        labs(x = "", y = "Rate (cM/Mb)") +
        theme(
            panel.background = element_rect(color = NULL, fill = alpha("grey", 0.10)),
            panel.grid = element_blank(),
            axis.line = element_blank(),
            axis.ticks = element_blank(),
            axis.text.x = element_blank(),
            axis.title.x = element_blank(),
            strip.text.y = element_text(face = "bold", colour = "white"),
            strip.background = element_rect(fill = "cornflowerblue"),
            plot.margin = plot_margin,
            aspect.ratio = 0.55
        )

    mhc_min <- 29640147
    mhc_max <- 33115544
    class_1 <- c(29640147, 31478898)
    class_3 <- c(31478898, 32191844)
    class_2 <- c(32191844, 33115544)

    p_recom <- if(chromosome == 6 & (window_min < mhc_max | window_max > mhc_min)) {
        p_recom +
        annotate(
            "text",
            y = 85,
            x = (min(mhc_max, window_max) + max(mhc_min, window_min)) / 2,
            size = 3,
            label = 'italic("Extended MHC chr6: 29.64 - 33.12 Mb")',
            parse = TRUE
        ) +
        annotate("rect", xmin = max(mhc_min, window_min),
                xmax = min(mhc_max, window_max), ymin = -Inf, ymax = Inf,
                alpha = .2, fill = "grey")
    } else {
        p_recom
    }
}


# NCBI RefSeq genes
plot_gene <- function(refGene, target_bp, chromosome, window_max, window_min,
                      window_kb, region_p_min, plot_margin) {

    genes <- glist_hg19 %>%
        dplyr::filter(chr == chromosome & end > window_min & start < window_max) %>%
        dplyr::mutate(
        panel = "RefSeq Genes",
        plot_loc = seq_len(nrow(.)),
        x_start = ifelse(start < window_min, window_min, start),
        x_end = ifelse(end > window_max, window_max, end))

    ggplot(genes) +
        facet_grid(panel ~ ., scales = "free") +
        geom_segment(aes(x = x_start, xend = x_end, y = plot_loc, yend = plot_loc),
                    color = "blue", na.rm = TRUE, size = 2) +
        geom_vline(xintercept = target_bp, color = "red", linetype = 3, size = 0.25) +
        scale_y_discrete(limits = range(genes$plot_loc)) +
        scale_x_continuous(limits = c(window_min, window_max)) +
        ggrepel::geom_text_repel(aes(y = plot_loc, x = (x_start + x_end) / 2,
                                 label = name), color = "grey", size = 2.5,
                                 segment.size = 0.25, na.rm = TRUE) +
        theme(
            panel.background = element_rect(color = NULL, fill = alpha("grey", 0.10)),
            panel.spacing = unit(1, "lines"),
            panel.grid = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.title = element_blank(),
            strip.text.y = element_text(face = "bold", colour = "white"),
            strip.background = element_rect(fill = "blue"),
            plot.margin = plot_margin,
            aspect.ratio = 2
        )
}


# CpG islands
plot_cpg <- function(cpgIslandExt, target_bp, chromosome, window_max,
                     window_min, window_kb, region_p_min, plot_margin) {
    cpgIslandExt %>%
        dplyr::select(chromStart, chromEnd, name) %>%
        dplyr::mutate(panel = "CpG") %>%
        ggplot() +
        geom_hline(yintercept = target_bp, color = "red", linetype = 3, size = 0.25) +
        facet_grid(panel ~ ., scales = "free") +
        geom_linerange(aes(ymin = chromStart, ymax = chromEnd, x = 1), size = 4,
                       color = "#17807E") +
        ggrepel::geom_text_repel(aes(x = 1, y = chromStart, label = name),
                                 point.padding = 1, segment.size = 0.125,
                                 color = "grey", size = 2.5, direction = "y") +
        labs(x = "", y = "") +
        scale_y_continuous(limits = c(window_min, window_max)) +
        coord_flip() +
        theme(
            panel.background = element_rect(color = NULL, fill = alpha("grey", 0.10)),
            panel.spacing = unit(1, "lines"),
            panel.grid = element_blank(),
            panel.border = element_blank(),
            axis.line = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.title = element_blank(),
            strip.text.y = element_text(face = "bold", colour = "white"),
            strip.background = element_rect(fill = "#17807E"),
            aspect.ratio = 0.5,
            plot.margin = plot_margin
        )
}


# Transcription factor binding sites
plot_tfbs <- function(wgEncodeRegTfbsClusteredV2, target_bp, chromosome,
                      window_max, window_min, window_kb, region_p_min,
                      plot_margin) {
    wgEncodeRegTfbsClusteredV2 %>%
        dplyr::mutate(panel = "TFBS") %>%
        ggplot() +
        facet_grid(panel ~ ., scales = "free") +
        geom_hline(yintercept = target_bp, color = "red", linetype = 3, size = 0.25) +
        geom_linerange(aes(ymin = chromStart, ymax = chromEnd, x = 1), size = 2,
                    color = "forestgreen", na.rm = TRUE) +
        labs(x = "", y = "") +
        scale_y_continuous(limits = c(window_min, window_max)) +
        coord_flip() +
        theme(
            panel.background = element_rect(color = NULL, fill = alpha("grey", 0.10)),
            panel.spacing = unit(1, "lines"),
            panel.grid = element_blank(),
            panel.border = element_blank(),
            axis.line = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.title = element_blank(),
            strip.text.y = element_text(face = "bold", colour = "white"),
            strip.background = element_rect(fill = "forestgreen"),
            plot.margin = plot_margin,
            aspect.ratio = 0.25
        )
}
