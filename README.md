# Replication Files

This directory reproduces the main tables and figures for *Rules Under Pressure: Euroscepticism and the Enforcement of EU Fiscal Discipline* by Xiao LU.

## Contents

- `replication_data.RData`: the main dataset.
- `main_analysis.R`: generates the main figures and tables.

## Running the code

Use this directory as the working directory, then run:

```sh
Rscript main_analysis.R
```
The required R packages are `ggplot2`, `lmtest`, `sandwich`, and `survival`. No external data or internet connection is required. Each analysis script writes its listed files to `outputs/`.
