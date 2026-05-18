# NicheNet — Differential vs Standard: What Changes and What You Lose

## Overview

The current `scripts/run-nichenet.R` is implemented as **differential NicheNet** — it compares signaling across two niches simultaneously using cross-niche DE. The alternative is **standard NicheNet**, where you run one independent analysis per niche and compare results manually afterward.

The choice matters because the two approaches ask different questions, use different code paths, and produce different outputs.

---

## Which Approach to Use

| | **Standard NicheNet (per niche)** | **Differential NicheNet (current code)** |
|---|---|---|
| **Question asked** | "What ligands best explain this receiver's DE genes?" | "What ligands are specifically enriched in niche A vs niche B?" |
| **Senders can be identical across niches** | Yes | No — breaks `calculate_niche_de` |
| **Run structure** | One independent run per niche | Single run comparing both niches simultaneously |
| **Ligand prioritization basis** | Ligand activity against receiver DE genes | Ligand activity + sender-side LFC + receptor-side LFC |
| **Circos plots** | Yes — per niche | Yes — but highlights niche-specific ligands via `top_niche` coloring |
| **Identifies niche-specific signaling** | Only by manual comparison after the fact | Built into the prioritization score |
| **Receiver DE genes** | Computed by comparing receiver vs all other cells (per niche) | Computed by comparing receiver across niches |
| **Cross-niche contrast in one ranked list** | No — two separate ranked lists | Yes — single ranked list of what is differential |
| **Strongest paper claim** | "These are the top interactions in each niche" | "These interactions are specifically elevated in niche A over niche B" |
| **Risk of code breakage with identical senders** | None | Yes — identical senders cause self-vs-self DE in `calculate_niche_de` |

---

## Why Identical Senders Break Differential NicheNet

`calculate_niche_de` computes pairwise DE by comparing each sender in niche 1 against each sender in niche 2. If the same cell type appears as a sender in both niches, the function constructs a self-vs-self comparison (e.g. `Monocytes vs Monocytes`). This causes the DE function to either throw an error or return LFC = 0 across the board, which collapses `min_lfc` — the core specificity metric — for every ligand in that sender.

---

## Code Flow Differences

The differential-specific steps in `run-nichenet.R` that get replaced or removed in the standard version:

| Step | Differential (current code, lines) | Standard replacement |
|---|---|---|
| Sender DE | `calculate_niche_de` + `process_niche_de` + `combine_sender_receiver_de` (lines 38–44) | Not needed — senders only used to get list of expressed ligands |
| Receiver DE / geneset | `calculate_niche_de_targets` + `process_receiver_target_de` + `niche_geneset_list` (lines 57–78) | `FindMarkers` on receiver vs all other cells — one call per niche |
| Ligand activity | `get_ligand_activities_targets` | `predict_ligand_activities` — same logic |
| Expression tables | DotPlot — same | DotPlot — same |
| Prioritization | `get_prioritization_tables` with full weight vector including `scaled_ligand_score` and `scaled_receptor_score` | Rank directly by ligand activity score — no unified prioritization function |
| Output structure | Single table across both niches with `niche` column | Separate ranked list per niche |

Lines that disappear entirely in the standard version: **38–44** (sender DE block) and **57–78** (target DE + geneset block).

---

## Output File Differences

| File | Differential | Standard |
|---|---|---|
| `nichenet-prioritization-lr.csv` | Has `niche`, `scaled_ligand_score` (from cross-niche sender LFC), `scaled_receptor_score` (from cross-niche receiver LFC), and a unified `prioritization_score` combining all weights | Replaced by a simpler per-niche ligand activity table — no sender LFC columns, no cross-niche `prioritization_score` |
| `nichenet-prioritization-lt.csv` | Ligand → target genes, niche-aware in a single table | Still exists, but produced independently per niche |
| Circos plots (`{niche}-circos.pdf`) | LR pairs are color-coded by `top_niche` — which niche each pair is most specific to, determined by comparing `prioritization_score` across niches | Still exists, but `top_niche` coloring is absent since there is no cross-niche comparison |

The `top_niche` assignment in the circos loop (lines 170–171 of `run-nichenet.R`) specifically depends on comparing `prioritization_score` across both niches to assign each LR pair to its dominant niche. This logic has no equivalent in the standard version.

---

## Summary

- If your two niches have **distinct senders** and you want a single ranked list of what is differentially active between them → keep differential NicheNet as-is
- If your senders overlap or are identical, or you just want per-niche top interactions for circos plots → switch to standard NicheNet (run independently per niche)
- In either case, circos plots are achievable; what you lose in the standard version is the cross-niche specificity scoring and the `top_niche` coloring
