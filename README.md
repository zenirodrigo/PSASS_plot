# PSASS_plot

![R](https://img.shields.io/badge/R-%3E%3D%204.1-276DC3?logo=r\&logoColor=white)
![Platform](https://img.shields.io/badge/platform-Linux%20%7C%20macOS-lightgrey)
![Output](https://img.shields.io/badge/output-PNG-success)
![License](https://img.shields.io/badge/license-MIT-green)

A script for visualization of **PSASS** results (https://github.com/SexGenomicsToolkit/PSASS)

## Overview

`PSASS_plot` uses the PSASS output file:

`psass_window.tsv`

and produces:

1. Multi-panel Manhattan plot
2. Circos plot

This script improves the visualization of genomic regions with:

* fixation index (Fst)
* depth ratio shifts
* SNP-density variation

---


## Installation

Clone the repository
```bash
git clone https://github.com/zenirodrigo/PSASS_plot.git
cd PSASS_plot
```

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

A simple way to run the script is by creating a dedicated Conda environment.

### 1. Create the environment

```bash
conda create -n psass_plot_env -c conda-forge -c bioconda -c defaults r-base=4.3 r-ggplot2 r-dplyr r-tidyr r-scales r-png r-patchwork r-sgtr
conda activate psass_plot_env
```


## Input Requirements

Required file:

`psass_window.tsv`

The file must be a tab-separated table generated from PSASS window analyses. No further formatting is needed if it was created from the PSASS original pipeline.

---

## Usage

This script must be executed inside the PSASS analysis output directory containing the file:
Parameters:

--n

Plots only the first N sequences/contigs from the ordered dataset.

Rscript run_psass_plot_auto.R --n 25

Example: use only the first 25 sequences for plotting.

--chr

Plots only the sequence/contig at index X in the ordered TSV.

Rscript run_psass_plot_auto.R --chr 11

Example: if --chr 11 is used, the script will plot only the 11th sequence in the dataset order.

--region

Must be used together with --chr.
Plots only a specific genomic region of the selected sequence.

Rscript run_psass_plot_auto.R --chr 11 --region 1:100000

Example: plot only positions 1 to 100000 from the 11th sequence.

Behavior summary
No arguments → plot all sequences normally
--n X → plot only the first X sequences
--chr X → plot only sequence X
--chr X --region A:B → plot only region A:B from sequence X
Important notes
--chr uses the numeric order of sequences in the TSV, not the chromosome name.
--region only works together with --chr.
--n and --chr should not be used together

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




Both figures are exported as PNG files.

## License

MIT License
