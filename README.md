# Analysis for NIST Integration pipeline V3.3.2 manuscript.

This repository contains the code used to generate Supplemental Table XYZ, Supplemental Figure XYZ, and Supplemental Table 2.  

The R package [drake](https://github.com/ropensci/drake) was used to organize the run the analysis. 
The `R` directory contains scripts with the functions used in the analysis pipeline (`functions.R`) and the drake plan (`plan.R`) which defines the analysis workflow. The analysis is in the `analysis_results.Rmd` Rmarkdown document. 