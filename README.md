# Wildlife trade drives animal-to-human pathogen transmission over 40 years

This repository contains the **data and R code** used in:

> **Gippet, J.M.W. et al. (2026)**. *Wildlife trade drives animal-to-human pathogen transmission over 40 years*.

---

## Data description

### Main dataset

**`Dataset_complete.RDS`**

A species-level dataset covering extant wild mammals, integrating information on:

* Occurrence in the global wildlife trade (legal, illegal, live animals, products)
* Zoonotic pathogen sharing with humans (from CLOVER and related sources)
* Research effort (bibliometric indices)
* Synanthropy and wild meat use
* Phylogenetic information (MDD identifiers and code for calculating phylogenetic eigenvectors from VertLife phylogenetic trees)
* Biogeographic realm and taxonomic classification

### Data table

**`DataTableS1_Gippet2025_WildlifeTradeZoonoses.xlsx`**

Data Table S1 referenced in the manuscript, providing detailed variable definitions and sources.

---

## R script

**`script_WTZ_dec2025.R`**

This script reproduces all analyses and raw versions of the figures presented in the paper, including:

1. Data cleaning and preprocessing
2. Phylogenetic tree processing and grafting of missing species
3. Computation of phylogenetic eigenvectors (PVR)
4. Statistical analyses:
   * Binomial and negative binomial GLMMs
   * Structural equation models (piecewise SEM)
5. Generation of main text and supplementary figures

The script is structured by figure and follows the order of the manuscript. 
The software Inkscape was used to create the final rendering of the figures.

---

## Requirements

### R version

* R ≥ 4.2.0 (recommended)

### Main R packages

Key packages used include:

* `glmmTMB`
* `piecewiseSEM`
* `ggplot2`, `ggtree`, `ComplexHeatmap`
* `phytools`, `ape`, `PVR`
* `dplyr`, `tidyr`, `purrr`

A full list of required packages is loaded at the beginning of the script.

---

## Reproducibility notes

* Random processes (e.g. phylogenetic grafting) use fixed seeds where relevant.
* Phylogenetic trees are sourced from **VertLife** (see script comments for download instructions).
* Results may vary slightly depending on platform and package versions.

---

## Citation

If you use this data or code, please cite:

> Gippet, J.M.W. et al. (2026). *Wildlife trade drives animal-to-human pathogen transmission over 40 years*.

---

## License

* **Code**: MIT License
* **Data**: CC-BY 4.0, unless otherwise specified by original data providers

Please respect the terms of use of the original data sources (CITES, LEMIS, CLOVER, VertLife).

---

## Contact

For questions or issues, please contact:

**Jérôme M. W. Gippet**
[jerome.gippet.com]
