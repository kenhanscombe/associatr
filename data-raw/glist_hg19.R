
library(tidyverse)
library(data.table)

# Gene range lists downloaded from PLINK 1.9 resources, generated from UCSC Table Browser RefSeq track https://www.cog-genomics.org/static/bin/plink/glist-hg19

glist_hg19 <- fread(
  "data-raw/glist-hg19",
  header = FALSE,
  col.names = c("chr", "start", "end", "name")
) %>% as_tibble()

usethis::use_data(glist_hg19, overwrite = TRUE)

