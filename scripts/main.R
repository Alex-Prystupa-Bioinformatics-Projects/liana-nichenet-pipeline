# main.R
library(yaml)
configs <- yaml::read_yaml("config.yml")
pipeline_configs <- yaml::read_yaml("pipeline-config.yml")
library(Seurat)
library(dplyr)
library(ggplot2)
library(glue)
library(openxlsx)
library(tibble)
library(liana)
library(nichenetr)
library(presto)
library(forcats)
library(cowplot)
library(RColorBrewer)
library(clusterProfiler)
library(org.Hs.eg.db)

source("scripts/helper-functions/helper-functions.R")
source("scripts/run-liana.R")
source("scripts/run-nichenet.R")
source("scripts/run-differential-nichenet.R")
source("scripts/run-liana-nn-combined.R")
options(timeout=600)

# Output base directory — named by cluster column and niche names
out_base <- glue("lrp-outs-{configs$paths$CellTypeCol}-{paste(names(configs$nichenet$niches), collapse='-')}")

# Create Output Directories (in helper-functions)
create.all.dirs(out_base)

#### Part 1: Load Data
rds_path <- configs$paths$RDS_File
if (!file.exists(rds_path)) {
  message(glue("ERROR: RDS file not found: {rds_path}"))
  quit(status = 1)
}
seu_obj <- readRDS(rds_path)

# Join Layers if not joined already
DefaultAssay(seu_obj) <- configs$paths$DefaultAssay
seu_obj <- JoinLayers(seu_obj, assay = configs$paths$DefaultAssay)

cell_type_col <- configs$paths$CellTypeCol
if (!cell_type_col %in% colnames(seu_obj@meta.data)) {
  message(glue("ERROR: Column '{cell_type_col}' not found in Seurat object metadata. Available columns: {paste(colnames(seu_obj@meta.data), collapse=', ')}"))
  quit(status = 1)
}
Idents(seu_obj) <- seu_obj@meta.data[[cell_type_col]]

# Subset to only cell types referenced in niches — LIANA/NicheNet don't need the full atlas
niche_celltypes <- unique(c(
    unlist(lapply(configs$nichenet$niches, function(n) c(n$sender, n$receiver)))
))
seu_obj_niche <- subset(seu_obj, idents = niche_celltypes)
message(glue("Subsetting to {length(niche_celltypes)} niche cell types: {paste(niche_celltypes, collapse=', ')}"))

#### Part 2: Run LIANA
run_liana(seu_obj_niche)

#### Part 3: Run NicheNet
if (configs$nichenet$mode == "differential") {
  run_differential_nichenet(seu_obj_niche)
} else {
  run_nichenet(seu_obj_niche)
}

#### Part 4: Combined LIANA + NicheNet Analysis
run_liana_nn_combined(seu_obj_niche)
