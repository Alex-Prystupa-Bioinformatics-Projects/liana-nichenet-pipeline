# NicheNet — Feature Extraction Guide (from NicheNet v2 Protocol Paper)

## Flowchart: How to Define Ligands and Gene Set of Interest

```
START
└─ For each receiver cell population
   ├─ Ask: Are the sender cell population(s) and their secreted ligands captured in the dataset?
   │  ├─ Yes
   │  │  └─ Potential ligands:
   │  │     Sender-focused approach:
   │  │     ligands expressed in sender(s) whose cognate receptors are expressed in the receiver
   │  └─ No / I don't know
   │     └─ Potential ligands:
   │        Sender-agnostic approach:
   │        all ligands in LR database whose cognate receptors are expressed in the receiver
   │
   └─ Then ask: Does the dataset have two or more conditions?
      ├─ Yes
      │  ├─ Gene set of interest:
      │  │  DE genes in condition of interest within the receiver
      │  └─ Background gene set:
      │     All expressed genes in the receiver
      │
      └─ No
         └─ Ask: Is the CCC event of interest one of the following?
            a) Cell differentiation:
               differences between a progenitor and differentiated cell type
            b) Cell localization:
               differences between cell subtypes located in different niches or areas
            ├─ Yes
            │  ├─ Gene set of interest:
            │  │  DE genes between the receiver and another cell (sub)type
            │  └─ Background gene set:
            │     The union of all expressed genes in both cell types
            │     OR
            │     All genes in the transcriptome
            │     (if gene set size is not ≥10× the gene set of interest)
            │
            └─ No
               └─ NicheNet analysis not suitable
```

---

## How This Maps to This Project

### Decision 1 — How to define potential ligands

> "Are the sender cell populations captured in the dataset?"

**Answer: Yes.** Senders are known — monocytes dominate hub_2 (squamous niche), Th17/ILC3/cDC2 dominate hub_4 (columnar niche).

**→ Sender-focused approach.** Only consider ligands expressed in those senders whose cognate receptors are expressed in the receiver. The current code already does this correctly.

---

### Decision 2 — How to define the gene set of interest

> "Does the dataset have two or more conditions?"

**Answer: No.** This is one integrated dataset — there is no treated vs untreated, no timepoint comparison.

> "Is the CCC event cell localization — differences between cell subtypes located in different niches or areas?"

**Answer: Yes.** The squamous niche (hub_2) and columnar niche (hub_4) are anatomically distinct locations within the perianal fistula. This is a localization comparison.

**→ The paper prescribes:**
- Gene set of interest = DE genes between the receiver cell type and the other receiver cell type (`FindMarkers` directly between the two receivers)
- Background = union of all expressed genes in both cell types (if that union is ≥10× the size of the gene set; otherwise use the full transcriptome)

---

## Implication for the Current Code

The current `run-nichenet.R` uses the **differential NicheNet** pipeline, which constructs the receiver gene set via `calculate_niche_de_targets` — a cross-niche DE approach. According to the protocol paper, this is not the prescribed path for a localization comparison in a single dataset.

The paper points toward **standard NicheNet run independently per niche**, where:

1. `FindMarkers` is called directly between the two receiver cell types to get the gene set of interest
2. The background is the union of expressed genes in both receivers
3. Ligand activity is scored against that gene set
4. No cross-niche sender DE is needed — which also resolves the identical sender problem entirely

This means the geneset construction no longer depends on `calculate_niche_de`, `calculate_niche_de_targets`, or `min_lfc`. The analysis becomes simpler, more aligned with the paper's recommendations, and avoids the code breakage risk from identical senders.

See `NICHENET-DIFFERENTIAL-VS-STANDARD.md` for a full breakdown of what changes in the code and what output files differ between the two approaches.
