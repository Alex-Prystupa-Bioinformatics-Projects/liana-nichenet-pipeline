# LIANA vs NicheNet — Interpretation Framework

## The Core Distinction

**LIANA and NicheNet are asking fundamentally different questions.**

### LIANA
**Question:** *"Which ligand-receptor pairs are actively expressed between these cell types?"*

- The output metric is `aggregate_rank` — this is a **p-value** from Robust Rank Aggregation (RRA)
- **Lower is better** — a pair with a low aggregate_rank was ranked consistently high across multiple tools (NATMI, CellChat, SingleCellSignalR, etc.)
- It is **expression-driven** — it tells you what interactions are happening, but not what they are doing downstream
- Use `aggregate_rank < 0.05` as your significance threshold

### NicheNet
**Question:** *"Of all the ligands a sender could be producing, which ones are most likely driving transcriptional changes in the receiver?"*

- The output metric is `prioritization_score` — combines multiple signals:
  - Ligand differential expression in the sender niche
  - Receptor differential expression in the receiver
  - Expression levels (ligand and receptor fraction)
  - **Ligand activity** — how well that ligand's known target genes predict the receiver's actual DE gene set
- It is **functionally-driven** — it adds the downstream biology that LIANA cannot see

---

## How to Interpret Them in Parallel (Without the Combined Method)

### Step 1 — Start with NicheNet
NicheNet gives you a ranked list of ligands most likely to be biologically meaningful in your niches (e.g. Squamous vs Columnar). This is your hypothesis list.

### Step 2 — Cross-validate with LIANA
For your top NicheNet ligands, check if they also have a low `aggregate_rank` in LIANA for the same sender-receiver pair. If a ligand is both:
- Highly prioritized by NicheNet (high `prioritization_score`) AND
- Has strong LIANA consensus (low `aggregate_rank`)

That is your most confident hit — strongest candidate for a paper claim.

### Step 3 — Interpret the gaps

| Scenario | Interpretation |
|---|---|
| Low LIANA rank + High NicheNet score | **Strongest hits** — expressed AND functionally active |
| Low LIANA rank only | Interaction is expressed but may not be driving receiver transcription — interesting but lower confidence |
| High NicheNet score only | Predicted functional activity but expression evidence is weaker — flag but treat cautiously |

> **Key caveat from the NicheNet docs:** *"A ligand can have a great activity score but we don't actually know if it's mediated by the receptor that LIANA predicted."* The overlap between the two is where you make your strongest mechanistic claims.

---

## Finding Downstream Target Genes

Once you have a ligand that appears in both LIANA and NicheNet, the next question is: **which genes is it actually regulating downstream in the receiver?**

This is answered by the **ligand-target matrix** in NicheNet. The pipeline already computes and saves this:

```
lrp-outs/niche-net/results/nichenet-prioritization-lt.csv
```

This file links each ligand to its predicted target genes in the receiver cell type. NicheNet scores ligand activity by asking how well a ligand's known downstream targets overlap with the receiver's actual DE gene set — so the targets in this table are the genes that drove the activity score.

### Full Workflow

1. Pull top LIANA hits (`aggregate_rank < 0.05`) filtered to your niche sender-receiver pairs
2. Cross with top NicheNet prioritized ligands (`prioritization_score`)
3. For overlapping ligands → filter `nichenet-prioritization-lt.csv` by ligand to get predicted target genes in the receiver
4. Run pathway enrichment (GO, KEGG, Reactome) on those target genes to get the biological story

**This is exactly what the combined LIANA+NicheNet method (`run-liana-nn-combined.R`) automates** — but steps 1-3 can be done manually with the CSVs right now as a sanity check before committing to that analysis.

---

## In the Context of This Project

- **Squamous niche (hub_2):** Monocyte-dominated, innate/myeloid skewed — expect NicheNet to prioritize myeloid-derived ligands acting on squamous epithelium
- **Columnar niche (hub_4):** Adaptive mucosal immune environment (Th17-MAIT, CD8aa IEL, ILC3, cDC2) — expect mucosal cytokine signaling (IL-17, IL-22 axis) to dominate

Cross-referencing LIANA and NicheNet hits between these two niches will let you make a direct comparison of the signaling logic governing squamous vs columnar epithelial maintenance in the fistula context.
