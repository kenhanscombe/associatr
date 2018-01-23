
# associatr: Visualization of Human Genetic Association Evidence

[![Travis-CI Build Status](https://travis-ci.org/kenhanscombe/associatr.svg?branch=master)](https://travis-ci.org/<USERNAME>/<REPO>)

<br>

### Overview

<br>

`associatr` is an R package for visualization of human genetic association results. Manhattan, qq, and regional plots, as well as some less commonly seen visualizations that enable the reader to evaluate association evidence. Emphasis is placed on a fast uninterupted workflow: biological context is retrieved with an `SQL` query and by default large p-values are trimmed. Individual plotting functions can be inserted into an analysis pipeline, and because the plots are all generated with `ggplot2`, they can be saved to any resolution, size, scale, and format with `ggsave`.

<br>
<br>




### Plots

<br>

_Genome-wide_

1. `gwa_manhattan`

2. `gwa_qq`


_Regional_

3. `gwa_region`

4. `gwa_signal` chisq / ld (rsq) signal support


_To do_

- `gwa_freq_effect` effect size, MAF, n at each frequency 

- `gwa_condition` conditional analysis series

- `gwa_hla` overlay hla alleles on SNP data

- `gwa_risk_score` PRS/ LD-score/ heritability/ genetic correlation

- `gwa_haplotype` a haplotype bifurcation plot

- `gwa_ld` ld with `target` SNP, in a `window`, retrieved or calculated from a `panel = c("reference", "local")`. Implement call to PLINK as R subprocess. 

- `gwa_ideogram`

<br>
<br>



### Input data

Update to accept standard formats:

* PLINK, BGENIE, SNPTEST

* PRSICE, LD-SCORE, GCTA?
