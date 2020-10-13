
# associatr: Genetic Association Plots

<!-- badges: start -->
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)
[![Travis build status](https://travis-ci.com/kenhanscombe/associatr.svg?branch=master)](https://travis-ci.com/kenhanscombe/associatr)
[![Codecov test coverage](https://codecov.io/gh/kenhanscombe/associatr/branch/master/graph/badge.svg)](https://codecov.io/gh/kenhanscombe/associatr?branch=master)
<!-- badges: end -->

## Overview

associatr is an R package for manhattan and regional plots. By default only variants with p < 0.001 are plotted in the manhattan plot. Region plots use LDlinkR to retrieve LD and SQL queries to UCSC genome browser to get biological context.

<br>

To install

```{R}
devtools::install_github("kenhanscombe/associatr", dependencies = TRUE, force = TRUE)
```

**Note**. `region` requires an [LDlinkR personal access token](https://cran.r-project.org/web/packages/LDlinkR/vignettes/LDlinkR_vignette_v8-2_TM.html). [Request token](https://ldlink.nci.nih.gov/?tab=apiaccess.).
