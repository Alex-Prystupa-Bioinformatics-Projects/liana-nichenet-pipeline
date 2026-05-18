# CLAUDE.md

## Overview

This is a **Ligand-Receptor Pathway (LRP) analysis pipeline** for cell-cell communication analysis in single-cell RNA-seq data. It integrates two complementary tools — **LIANA** and **NicheNet** — to identify and prioritize intercellular signaling in single cell datasets.

The main paper to follow for ideation with user on what analysis steps for NICHENET that is best is at this link: https://www.nature.com/articles/s41596-024-01121-9

The goals are:
    1. Run LIANA which is a package that aggreggates results across multiple cell-cell communication packages
      a. Liana REFERENCE VIGNETTE: https://saezlab.github.io/liana/articles/liana_tutorial.html
      b. Liana GITHUB PAGE: https://github.com/saezlab/liana

    2. Run NicheNet
      a. NicheNet REFERENCE VIGNETTE: https://github.com/saeyslab/nichenetr/blob/master/vignettes/differential_nichenet.md
      b. NicheNet GITHUB PAGE: https://github.com/saeyslab/nichenetr

    3. Combine LIANA & NicheNet
      a. Combine LIANA & NicheNet REFERENCE VIGNETTE: https://saezlab.github.io/liana/articles/liana_nichenet.html

VERY CRITICAL: BEFORE implementing ANY analysis step involving these tools, you MUST fetch and read
the relevant vignette link above. Do not assume you know the correct approach from memory.
Do not implement based on prior knowledge alone. READ THE LINK FIRST, every single time.
If the user provides a vignette or documentation URL during the conversation, fetch and read it
immediately before writing any code. Failing to do this has caused major implementation errors.


## Architecture
```
scripts/main.R                  # Entry point: loads Seurat object, dispatches to analysis scripts
scripts/run-liana.R             # LIANA LR interaction analysis → CSV + PDF bar plots (DONE)
scripts/run-nichenet.R          # NicheNet ligand activity analysis → circos plots (WORK IN PROGRESS)
scripts/run-liana-nn-combined.R # Placeholder (not yet implemented)
scripts/helper-functions/
  helper-functions.R            # Shared utilities: dir creation, sender-receiver combos, plotting
routes/run-liana.sh             # LSF job script (resource specs + Rscript call)
config.yml                      # Pipeline configuration
init                            # One-time environment setup script
run                             # Job submission entry point
```

## Personal Coding Preferences
I like to have comments on top of big chunks of code. For example comments infront of large lapplys, and all I like to section out comments 1. 2. 3. Sometimes if they have steps

Always have seurat objects be named seu_obj . Any deviation from that should start with seu_obj ex. seu_obj_list, seu_obj_sub etc.
