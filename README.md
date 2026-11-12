# 100days_group3

> **Archived.** No longer maintained; last updated December 2022.

Analysis code from group 3 of a "100 days" outbreak analysis exercise in December 2022, working on a line list, contact tracing data and pathogen sequences from an outbreak.

- `R/finlay_script.R` and `R/finlay_functions.R`: data cleaning, transmission chains with `epicontacts`, a phylogeny and transmission tree, serial interval and offspring distribution estimates.
- `epinow2.r` (and `R/epinow2.r`): estimation of the reproduction number with `EpiNow2`, using the serial interval from the script above.
- `R/fix_onset_date.R`: corrects onset dates in the line list.
- `data-raw/`: the line list, contact tracing data, sequences and tree file, and saved delay distributions.
- `figures/`: resulting plots, including the phylogenetic and transmission trees and offspring distributions.
