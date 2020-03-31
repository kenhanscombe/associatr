
# associatr: Genetic Association Plots

<!-- badges: start -->
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)
[![Travis build status](https://travis-ci.org/kenhanscombe/associatr.svg?branch=master)](https://travis-ci.org/kenhanscombe/associatr)
<!-- badges: end -->

## Overview

`associatr` is an R package for manhattan and regional plots. By default only variants with p < 0.001 are plotted in the manhattan plot. Region plots use LDlinkR to retrieve LD and SQL queries to UCSC genome browser to get biological context.
