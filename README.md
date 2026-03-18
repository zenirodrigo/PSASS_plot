# PSASS_plot

![R](https://img.shields.io/badge/R-%3E%3D%204.1-276DC3?logo=r\&logoColor=white)
![Platform](https://img.shields.io/badge/platform-Linux%20%7C%20macOS-lightgrey)
![Output](https://img.shields.io/badge/output-PNG-success)
![License](https://img.shields.io/badge/license-MIT-green)

A script  for **PSASS** visualization (https://github.com/SexGenomicsToolkit/PSASS), generating:

## Overview

`PSASS_plot` takes the standard PSASS output file:

`psass_window.tsv`

and produces:

1. Multi-panel Manhattan plot
2. Circos plot

This script improves visualization of genomic regions with:

* copy number variation signals
* depth ratio shifts
* SNP-density variation

---


## Installation

This script requires **R** and the following packages:

- `sgtr`
- `ggplot2`
- `dplyr`
- `patchwork`
- `tidyr`
- `scales`
- `png`
- `grid` *(base R package, no separate installation required)*

---

## Conda Environment

A simple way to run the script is to create a dedicated Conda environment.

### 1. Create the environment

```bash
conda create -n psass_plot_env -c conda-forge r-base=4.3 r-ggplot2 r-dplyr r-tidyr r-scales r-png r-patchwork r-sgtr
conda activate psass_plot_env
```


## Input Requirements

Required file:

`psass_window.tsv`

The file must be a tab-separated table generated from PSASS window analyses.

---

## Usage

This script must be executed inside the PSASS analysis directory output containing the  file:

`psass_window.tsv`

Examples:

`Rscript run_psass_plot_auto.R --n 25`

`Rscript run_psass_plot_auto.R`

Behavior:

- if `--n` is provided, the script keeps only the first `n` ordered contigs
- if `--n` is omitted, the script uses all available contigs



## Outputs

The script generates:

* circos_depth_ratio.png
<img width="2400" height="2400" alt="HAP2_PSASS_circos_FST_SNPf_SNPm_pastelAlt" src="https://github.com/user-attachments/assets/e4dc57e6-aeb8-4471-a615-86edfec80462" />

* manhattan_depth_ratio.png
<img width="4800" height="5100" alt="HAP2_PSASS_manhattan_FST_SNPf_SNPm_pastelAlt" src="https://github.com/user-attachments/assets/dcec5159-666e-410c-bc72-6696b2b985c7" />

---


The script generates both intermediate files and final figures in the current working directory.

Intermediate output files

psass_window.clean.tsv
Cleaned version of the input file after removing a duplicated header line, if present.

selected_contigs.txt
List of contigs selected for plotting after ordering and optional filtering with --n.

psass_window.filtered.tsv
Filtered table containing only the selected contigs.

chromosomes.tsv
Contig-to-label mapping file used for Circos plotting.

Final output figures

The final figure names depend on the detected folder name:

HAP1_PSASS_manhattan_FST_SNPf_SNPm_pastelAlt.png

HAP1_PSASS_circos_FST_SNPf_SNPm_pastelAlt.png

or

HAP2_PSASS_manhattan_FST_SNPf_SNPm_pastelAlt.png

HAP2_PSASS_circos_FST_SNPf_SNPm_pastelAlt.png

or, if no haplotype-specific pattern is detected in the folder path:

PSASS_manhattan_FST_SNPf_SNPm_pastelAlt.png

PSASS_circos_FST_SNPf_SNPm_pastelAlt.png


Both figures are exported as PNG files.

## License

MIT License
