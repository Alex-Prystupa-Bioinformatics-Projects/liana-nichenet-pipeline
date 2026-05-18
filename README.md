# LRP — Ligand-Receptor Pathway Analysis Pipeline

A single-cell RNA-seq cell-cell communication pipeline integrating **LIANA** and **NicheNet** to identify and prioritize intercellular signaling across user-defined cellular niches.

## Overview

This pipeline is designed around a two-tool complementary strategy:

| Tool | Question | Output |
|---|---|---|
| **LIANA** | Which ligand-receptor pairs are actively expressed between these cell types? | `aggregate_rank` (lower = more consensus across tools) |
| **NicheNet** | Which ligands are most likely driving transcriptional changes in the receiver? | `prioritization_score` (combines expression + functional activity) |

Overlapping hits between LIANA and NicheNet are the strongest candidates for mechanistic claims.

## Directory Structure

```
config.yml                      # User-facing configuration (niches, thresholds, file paths)
pipeline-config.yml             # Internal paths (NicheNet reference matrices)
scripts/
  main.R                        # Entry point: loads Seurat object, dispatches analyses
  run-liana.R                   # LIANA LR interaction analysis → CSV + bar plots
  run-nichenet.R                # NicheNet ligand activity → circos plots
  run-liana-nn-combined.R       # Combined LIANA + NicheNet prioritization
  run-differential-nichenet.R   # Differential NicheNet (cross-niche comparison mode)
  helper-functions/
    helper-functions.R          # Shared utilities: dir creation, plotting helpers
routes/
  run-liana.sh                  # LSF job submission script
notes/                          # Interpretation guides and analysis notes
init                            # One-time renv environment setup
run                             # Job submission entry point
```

## Configuration

All user-facing settings live in `config.yml`. The key sections are:

```yaml
paths:
  RDS_File: /path/to/seurat_object.RDS
  DefaultAssay: RNA
  CellTypeCol: Final_Annotations   # metadata column holding cell type labels

nichenet:
  mode: standard        # "standard" (per-niche) or "differential" (cross-niche)
  organism: human       # or "mouse"
  lfc_cutoff: 1.5
  expression_pct: 0.25
  top_n_ligands: 40
  top_n_receptors_per_ligand: 3
  top_n_targets: 250

  niches:
    Squamous_niche:
      sender: [Tc1, Tc1-GZMA, Th17-MAIT]
      receiver: Squamous_Epithelia
    Columnar_niche:
      sender: [Tc1, Tc1-GZMA, Th17-MAIT]
      receiver: Columnar_Epithelia

combined:
  top_n_liana_interactions: 250
  aggregate_rank_threshold: 0.10
  top_n_nn_ligands: 25
```

Each niche has one receiver (unique across niches) and one or more senders. Cell type names must exactly match labels in the `CellTypeCol` metadata column.

## Usage

**1. Set up the environment (first time only)**
```bash
bash init
```

**2. Edit `config.yml`** to point to your Seurat RDS file and define your niches.

**3. Submit the job**
```bash
bash run
```

This submits `routes/run-liana.sh` to the LSF scheduler, which calls `Rscript scripts/main.R`.

**4. Outputs** are written to a directory named by your cell type column and niche names:
```
lrp-outs-{CellTypeCol}-{niche1}-{niche2}/
  liana/
  niche-net/
  liana-nn-combined/
```

## Dependencies

- R packages: `Seurat`, `liana`, `nichenetr`, `presto`, `clusterProfiler`, `org.Hs.eg.db`, `ggplot2`, `dplyr`, `yaml`, `glue`, `openxlsx`, `cowplot`, `RColorBrewer`
- NicheNet reference matrices (human/mouse) at paths defined in `pipeline-config.yml`
- HPC environment with LSF scheduler

## Biological Context

This pipeline was developed to study **columnar-to-squamous epithelial metaplasia in perianal fistula** (Simmons/Cho dataset). The niche design captures T cell → epithelial signaling across squamous and columnar niches, with the central question being what drives and maintains the metaplastic transition inside the fistula tract.

## References

- LIANA: [Tutorial](https://saezlab.github.io/liana/articles/liana_tutorial.html) | [GitHub](https://github.com/saezlab/liana)
- NicheNet: [Vignette](https://github.com/saeyslab/nichenetr/blob/master/vignettes/differential_nichenet.md) | [GitHub](https://github.com/saeyslab/nichenetr)
- Combined LIANA + NicheNet: [Vignette](https://saezlab.github.io/liana/articles/liana_nichenet.html)
- NicheNet methods paper: [Nature Protocols 2024](https://www.nature.com/articles/s41596-024-01121-9)
