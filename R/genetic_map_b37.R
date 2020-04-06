
#' Get genetic build 37 recombination map
#'
#' Retrieves \href{https://www.internationalgenome.org/category/population/}{1000 Genomes phase 3 recombination map} data for a specified population.
#'
#' @param pop A \href{https://www.internationalgenome.org/category/population/}{1000 Genomes population code}. Default "CEU".
#' @param write_map Write a serialized dataframe of the map data to disk. Default \code{FALSE}, does not write the map to disk but instead returns a dataframe.
#'
#' @return A dataframe. If \code{write_map} is \code{TRUE}, then the recombination map is written to the working directory as genetic_map_b37_<\code{pop}>.rds. Note. the map is ~22MB.
#'
#' @details Previously retrieved map data from https://mathgen.stats.ox.ac.uk/impute/1000GP_Phase3.tgz
#'
#' @importFrom stringr str_interp
#' @importFrom purrr map reduce
#' @importFrom data.table fread
#' @importFrom magrittr "%>%"
#'
#' @export
gwa_get_map <- function(pop = "CEU", write_map = FALSE) {

  population_code <- c("CHB", "JPT", "CHS", "CDX", "KHV", "CEU", "TSI", "FIN",
                       "GBR", "IBS", "YRI", "LWK", "GWD", "MSL", "ESN", "ASW",
                       "ACB", "MXL", "PUR", "CLM", "PEL", "GIH", "PJL", "BEB",
                       "STU", "ITU")

  if(!(pop %in% population_code)) {
    stop("Invalid population code. See https://www.internationalgenome.org/faq/which-populations-are-part-your-study/", call. = FALSE)
  }

  url <- stringr::str_interp("ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130507_omni_recombination_rates/${pop}_omni_recombination_20130507.tar")

  file <- basename(url)
  utils::download.file(url, file)

  tmp <- tempdir()
  utils::untar(file, exdir = tmp)
  list.files(tmp)

  map_files <- list.files(file.path(tmp, pop))

  genetic_map_b37 <- map_files %>%
    purrr::map(
      ~mutate(
        data.table::fread(file.path(tmp, pop, .), col.names = c(
          "position", "combined_rate", "genetic_map", "filtered")),
        chr = as.integer(
          gsub(stringr::str_interp("${pop}-|-final.txt.gz"), "", .)))) %>%
    purrr::reduce(rbind)

  file.remove(str_interp("${pop}_omni_recombination_20130507.tar"))
  unlink(tmp, recursive = TRUE)

  if (write_map) {
    saveRDS(genetic_map_b37, stringr::str_interp("genetic_map_b37_${pop}.rds"))
  } else {
    genetic_map_b37
  }
}
