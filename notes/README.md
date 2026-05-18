# LRP Pipeline

## NicheNet Parameter Reference

| Parameter | What it does | Importance | Typical values |
|---|---|---|---|
| `lfc_cutoff` | Minimum log fold-change for a gene to enter the geneset of interest (receiver upregulated genes used to score ligand activity). The resulting geneset should be 20–1000 genes. | High — directly controls which genes NicheNet tries to explain with ligand activity | `0.15` for 10x with 3+ niches; `0.25` for only 2 niches; higher for Smart-seq2 |
| `expression_pct` | Minimum fraction of cells in a cell type that must express a gene for it to be included. Genes below this threshold receive reduced prioritization scores. | Medium — prevents sporadic/lowly expressed genes from inflating scores | `0.10` (default); raise to `0.20` for stricter filtering |
| `top_n_ligands` | Number of top-ranked ligands shown per niche in the circos plot. | Low — visualization only, does not affect underlying scores | `15` (vignette default); increase if you want a denser plot |
| `top_n_receptors_per_ligand` | Number of top receptors shown per ligand in the circos plot. | Low — visualization only | `2` (vignette default) |
| `top_n_targets` | Number of top predicted target genes per ligand used in activity scoring. | Low-Medium — controls compute and focus, does not change biology much | `250` (vignette default) |

## NicheNet Prioritization Weights Reference

These are hardcoded in `scripts/run-nichenet.R` and should not need to be changed for standard non-spatial 10x analyses.

| Weight | What it scores | Value |
|---|---|---|
| `scaled_ligand_score` | Niche-specific DE of the ligand (min LFC vs all other niches) — most important signal | 5 |
| `scaled_ligand_expression_scaled` | Scaled expression level of the ligand across all cell types | 1 |
| `ligand_fraction` | Fraction of sender cells actually expressing the ligand | 1 |
| `scaled_ligand_score_spatial` | Spatial DE specificity of the ligand — zeroed (no spatial data) | 0 |
| `scaled_receptor_score` | Niche-specific DE of the receptor | 0.5 |
| `scaled_receptor_expression_scaled` | Scaled expression level of the receptor across all cell types | 0.5 |
| `receptor_fraction` | Fraction of receiver cells expressing the receptor | 1 |
| `ligand_scaled_receptor_expression_fraction` | Receptor expression strength within the ligand's binding family | 1 |
| `scaled_receptor_score_spatial` | Spatial DE specificity of the receptor — zeroed (no spatial data) | 0 |
| `scaled_activity` | Absolute ligand activity score — not recommended raw | 0 |
| `scaled_activity_normalized` | Normalized ligand activity relative to other ligands | 1 |
