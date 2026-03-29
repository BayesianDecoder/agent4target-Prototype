# Agent4Target — Prototype

**Agent-based Evidence Aggregation for Therapeutic Target Identification**

GSoC 2026 Proposal Prototype · UC OSPO / UCI · Mentor: [Ziheng Duan](https://ucsc-ospo.github.io/author/ziheng-duan/)

[![Python](https://img.shields.io/badge/Python-3.10%2B-blue)](https://python.org)
[![License: MIT](https://img.shields.io/badge/License-MIT-green)](LICENSE)

---

## What This Is

This prototype demonstrates a modular, agent-based pipeline for therapeutic target identification. Given any protein-coding gene as input, it:

1. Collects evidence from three independent biomedical sources in parallel
2. Harmonises heterogeneous disease identifiers to a unified MONDO namespace
3. Scores gene-disease associations using two complementary methods
4. Produces structured, source-traceable explanations for every score

Validated on **BRAF (ENSG00000157764)** against five FDA-approved drug indications.

---

## Pipeline Architecture

```
Gene Input (e.g. BRAF)
        │
        ├──► Agent 1 — OpenTargetsAgent     (GraphQL · per-datasource scores)
        ├──► Agent 2 — PharosAgent          (TDL druggability tier · MONDO IDs)
        └──► Agent 3 — DepMapAgent          (CRISPR essentiality · cancer lineages)
                │
                ▼
     DiseaseAlignmentAgent
     EFO → MONDO via static map + OLS4 API
                │
                ▼
     ┌─────────────────────┐  ┌──────────────────────────┐
     │  ScoringAgent        │  │  SimpleScoringAgent       │
     │  Noisy-OR + Shrinkage│  │  Confidence-weighted avg  │
     │  → Tier A / B / C    │  │  (reference baseline)     │
     └─────────────────────┘  └──────────────────────────┘
                │
                ▼
     ExplanationAgent
     Structured provenance — every score linked to source + evidence type + confidence
                │
                ▼
     Ranked gene-disease association list
     score · tier · n_sources · jackknife CI · novel flag · JSON export
```

---

## Key Design Decisions

| Decision | Rationale |
|---|---|
| **Per-datasource scores** (not aggregated) | Preserves evidence provenance — the exact limitation that motivated this project over DrugnomeAI |
| **Source-specific confidence weights** | ChEMBL (0.95) ≠ EuropePMC (0.55) — not all datasets are equally reliable |
| **Noisy-OR fusion** | Models each source as independent noisy channel — multiple weak signals jointly produce strong evidence |
| **Shrinkage toward prior 0.22** | Regularises sparse evidence using empirical druggable genome base rate (Finan et al. 2017) |
| **Tier system** | Tier A (n_eff ≥ 2.0) / Tier B (0.8–2.0) / Tier C (<0.8, suppressed) — prevents overconfident single-source rankings |
| **Structured explanations** | No free-text LLM generation — every score component traces to a specific EvidenceRecord |
| **EFO→MONDO alignment** | Static map (38 curated pairs, ~95% coverage) + OLS4 xref fallback — prerequisite for multi-source merging |

---

## Benchmark Results — BRAF

Ground truth: 4 FDA-approved BRAF drug indications (melanoma, non-small cell lung carcinoma, colorectal cancer, thyroid carcinoma)

| Strategy | Hits | Recall@10 | Role |
|---|---|---|---|
| Uniform average | 3 | 75% | naive baseline |
| **Confidence-weighted average** | **3** | **75%** | **reference — current default** |
| Noisy-OR post-alignment | 2 | 50% | primary — underperforms due to incomplete alignment |
| Max score | 3 | 75% | upper bound |

**Key finding:** Noisy-OR scores 50% because disease alignment is still incomplete (EFO subterm merging not fully covered). The confidence-weighted method is the reliable current ranker. The Noisy-OR architecture is correct — it requires ≥3 aligned sources to outperform the simple baseline. This will be the focus of the GSoC work.

---

## Unified Evidence Schema

Every agent outputs records conforming to a single Pydantic schema:

```python
class EvidenceRecord(BaseModel):
    target_id:        str          # Ensembl gene ID
    target_symbol:    str          # Gene symbol
    disease_id:       str          # Disease ID (EFO/MONDO)
    disease_name:     str
    evidence_type:    EvidenceType # genetic_association | somatic_mutation | known_drug | ...
    source:           EvidenceSource  # open_targets | depmap | pharos
    raw_score:        float        # 0.0 – 1.0
    confidence:       float        # source reliability weight
    calibrated_score: Optional[float]
    metadata:         Dict         # source-specific extras (TDL, CERES, datasource)
```

---

## Installation

```bash
git clone https://github.com/BayesianDecoder/agent4target-Prototype
cd agent4target-Prototype

pip install pydantic requests numpy pandas matplotlib seaborn
```

**Optional — DepMap agent:**

Download `CRISPRGeneEffect.csv` and `sample_info.csv` from [DepMap portal](https://depmap.org/portal/download/) and place them in the project root.
Data: DepMap CRISPRGeneEffect.csv (download from official source)

---

## Running the Notebook

```bash
jupyter notebook Agent4Target_Prototype.ipynb
```

Run all cells top to bottom. The notebook will:
- Fetch Open Targets evidence for BRAF via GraphQL
- Run PHAROS druggability query
- (If DepMap CSV is present) aggregate CRISPR scores by cancer lineage
- Align all disease IDs to MONDO
- Score with both methods
- Print the Recall@10 benchmark
- Export `BRAF_agent4target_report.json`

Without DepMap CSV, the notebook runs on Open Targets + PHAROS only and still produces meaningful rankings.



---

## Project Structure

```
agent4target-Prototype/
├── Agent4Target_Prototype.ipynb   ← main notebook (all agents + scoring + benchmark)
├── requirements.txt
├── README.md
└── results/                       ← created on run
    ├── BRAF_agent4target_report.json
    ├── braf_evidence_scores.png
    ├── braf_confidence_weighted_scores.png
    └── braf_evidence_heatmap.png
```

---

## Scoring — How It Works

### Source calibration (Layer 1)
Each source's raw score is transformed into a calibrated probability using source-specific functions:
- **Open Targets** — piecewise linear, anchored at raw=0.4 → P=0.22 (empirical prior)
- **DepMap** — sigmoid on CERES score, anchored at -1.0
- **PHAROS** — ordinal tier mapping: Tclin→1.00, Tchem→0.62, Tbio→0.28, Tdark→0.05

### Noisy-OR fusion (Layer 2)
```
P = 1 − Π(1 − wᵢ × pᵢ)
```
Source weights: PHAROS=0.40, DepMap=0.35, Open Targets=0.25

### Shrinkage toward prior (Layer 3)
```
α = K / (n_eff + K),   K = 1.5
final_score = α × 0.22 + (1 − α) × noisy_or_score
```
Where `n_eff = Σ confidence_i` across best record per source.

### Tier allocation
| Tier | n_eff | Meaning |
|---|---|---|
| A | ≥ 2.0 | Strong multi-source — evidence dominates |
| B | 0.8–2.0 | Moderate — prior dominates, flagged |
| C | < 0.8 | Suppressed — score set to None |

---

## Disease Alignment — Five-Step Resolution

| Step | Method | Coverage |
|---|---|---|
| 1 | MONDO passthrough | IDs already canonical |
| 2 | Static EFO→MONDO map | ~95% of common OT disease IDs |
| 3 | OLS4 xref lookup (`annotation.oboInOwl:hasDbXref`) | Long-tail EFO IDs |
| 4 | Label matching | Cases where xref index is incomplete |
| 5 | Graceful fallback + warning flag | Never breaks the pipeline |

> **Fixed bug:** Previous OLS4 implementation used `queryFields=obo_xref` (invalid field, silently returns empty). Correct field: `annotation.oboInOwl:hasDbXref`.

---


## References

1. Raies et al. (2022). DrugnomeAI. *Communications Biology* — https://www.nature.com/articles/s42003-022-04245-4
2. Finan et al. (2017). The druggable genome. *Science Translational Medicine* — empirical prior 0.22
3. Open Targets Platform — https://platform.opentargets.org
4. DepMap — https://depmap.org/portal
5. PHAROS / TCRD — https://pharos.nih.gov
6. MONDO Ontology — https://mondo.monarchinitiative.org

---

## Author

**A Vijay Aditya** · IIT (BHU) Varanasi  
GitHub: [BayesianDecoder](https://github.com/BayesianDecoder)  
Email: vijayaditya22012004@gmail.com
