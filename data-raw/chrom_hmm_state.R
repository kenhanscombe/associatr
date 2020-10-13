
library(tidyverse)

# 15 hmm states
state <- read_csv(
  "data-raw/broad_chromhmm_state.csv",
  col_names = c("state", "color", "annotation")
)

state_col <- c("red1", alpha("red1", 0.65), "magenta3", rep("orange1", 2),
               rep("yellow", 2), "dodgerblue", rep("green4", 2), "lawngreen", "grey35", rep("grey85", 3))


# 9 cell types
html_doc <- XML::htmlParse(
  "http://rohsdb.cmb.usc.edu/GBshape/cgi-bin/hgEncodeVocab?term=GM12878,H1-hESC,HepG2,HUVEC,HMEC,HSMM,K562,NHEK,NHLF")

html_table_nodes <- XML::getNodeSet(html_doc, "//table")

cell <- XML::readHTMLTable(
  html_table_nodes[[1]],
  as.data.frame = TRUE,
  stringsAsFactors = FALSE
) %>%
  as_data_frame() %>%
  rename_all(tolower) %>%
  select(cell, lineage, tissue, sex, everything())