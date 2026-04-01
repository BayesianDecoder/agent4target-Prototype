# Agent4Target — Prototype

**Agent-based Evidence Aggregation for Therapeutic Target Identification**

GSoC 2026 Proposal Prototype · UC OSPO / UCI · Mentor: [Ziheng Duan](https://ucsc-ospo.github.io/author/ziheng-duan/)

[![Python](https://img.shields.io/badge/Python-3.10%2B-blue)](https://python.org)
[![License: MIT](https://img.shields.io/badge/License-MIT-green)](LICENSE)

---

## What this is

Agent4Target is a prototype pipeline that automatically gathers evidence from three biomedical databases — Open Targets, DepMap, and PHAROS — to rank which diseases a gene is most likely to be a therapeutic target for.

The core idea: instead of flattening everything into a feature matrix (as DrugnomeAI does), each data source gets its own dedicated agent that preserves source identity, evidence type, and confidence. Scores are then fused using a probabilistic model that rewards multi-source convergence and penalizes single-source inflation.

The validation case is **BRAF** — a well-characterised oncology gene with 4 FDA-approved drug indications (melanoma, non-small cell lung carcinoma, colorectal cancer, thyroid carcinoma). Every scoring method is benchmarked against these ground truth diseases using Recall@20.

---

## Pipeline at a glance

```
Gene Symbol + Ensembl ID
        │
        ├─► OpenTargetsAgent     →  137 records (30 diseases, GraphQL API)
        ├─► DepMapAgent          →   28 records (28 cancer lineages, local CSV)
        └─► PharosAgent         →   20 records (20 diseases, GraphQL API)
                                         │
                              MondoNormalizationAgent
                         (68 raw disease IDs → 35 MONDO IDs)
                                         │
                    ┌────────────────────┴────────────────────┐
             ScoringAgent                             SimpleScoringAgent
        (Noisy-OR + shrinkage                    (confidence-weighted avg)
           + Tier A/B/C)                              [reference]
                    │
             ExplanationAgent
        (structured provenance JSON)
                    │
          Ranked disease list + JSON export
```

---

## Quickstart

### 1. Install dependencies

```bash
pip install pydantic requests numpy pandas matplotlib seaborn
```

### 2. Get the data files

**DepMap CRISPR Gene Effect** (required for DepMap agent):

```bash
# Download from https://depmap.org/portal/download/all/
# File: CRISPRGeneEffect.csv  (~200 MB)
# Place it in the same directory as the notebook
```

**DepMap sample metadata** (optional — needed for cancer lineage breakdown):

```bash
# Same download page: sample_info.csv
# Without it, lineage labels fall back to cell line IDs
```

Open Targets and PHAROS are accessed via live API — no download needed.

### 3. Run the notebook

```bash
jupyter notebook Agent4Target_Prototype.ipynb
```

Run all cells in order. The full pipeline (including API calls) takes about 30–60 seconds depending on network speed.

### 4. Try a different gene

```python
# Cell 8 — change these two lines
BRAF_ENSEMBL_ID = "ENSG00000157764"   # ← replace with any Ensembl ID
GENE_SYMBOL     = "BRAF"              # ← replace with the gene symbol
```

Then re-run the pipeline from cell 8 onwards.

---

## Notebook structure

The notebook has evolved during iteration, so section-level structure is more reliable than fixed cell ranges.

| Section | What it does |
|------|--------------|
| Setup & Imports | Loads dependencies and shared utilities |
| Evidence Schema | Defines `EvidenceRecord`, explanation schema, and final report schema |
| Agent 1 — Open Targets | Queries GraphQL associations and converts datasource scores to evidence records |
| Open Targets Tables | Full aligned table + de-duplicated unique views |
| Agent 2 — DepMap | Aggregates CRISPR dependency by lineage and derives normalized score/confidence |
| Agent 3 — PHAROS | Fetches TDL and disease associations via GraphQL |
| Disease Normalization | Route-based MONDO harmonization + second-pass fallback + targeted overrides |
| Scoring Agents | Advanced Noisy-OR + shrinkage + tiering and reference confidence-weighted scoring |
| Explanation Agent | Human-readable and structured provenance breakdown |
| Pipeline Orchestrator | Executes agents end-to-end and tracks stage status |
| Results & Visualizations | Ranked outputs, bar charts, and evidence heatmap |
| Benchmark | Expanded MONDO candidate-set `Recall@20` validation |
| Export & Summary | JSON report generation and execution status summary |

---

## Evidence sources

### Open Targets Platform
- **API:** `https://api.platform.opentargets.org/api/v4/graphql`
- **What it provides:** per-datasource evidence scores linking genes to diseases across 17 evidence types (genetics, somatic mutations, known drugs, literature, pathways)
- **Output for BRAF:** 137 records across 30 diseases
- **Confidence weights:** each datasource gets its own reliability weight — `chembl` (0.95) and `eva_somatic` (0.90) rank highest; `europepmc` (0.55) lowest

### DepMap Cancer Dependency Map
- **Source:** locally downloaded `CRISPRGeneEffect.csv`
- **What it provides:** CRISPR knockout essentiality scores across 1,186 cancer cell lines, grouped by cancer lineage
- **Key metric:** CERES score (negative = essential = good drug target candidate)
- **Output for BRAF:** 28 cancer type records; skin neoplasm scores highest (mean CERES −0.34)
- **Requires:** `CRISPRGeneEffect.csv` + optionally `sample_info.csv` from [depmap.org/portal/download](https://depmap.org/portal/download/all/)

### PHAROS / TCRD
- **API:** `https://pharos-api.ncats.io/graphql`
- **What it provides:** Target Development Level (TDL) — a classification of how druggable a gene is
- **TDL tiers and what they mean:**

| Tier | Meaning | Confidence |
|------|---------|-----------|
| Tclin | Approved drug exists with known mechanism | 0.95 |
| Tchem | Active compound in ChEMBL or DrugCentral | 0.75 |
| Tbio | Biological annotation only, no drug | 0.50 |
| Tdark | Poorly characterised, essentially unknown | 0.20 |

- **BRAF result:** Tclin — maximum confidence, approved drugs exist

---

## Scoring methods

### Primary — Noisy-OR + Shrinkage + Tiers

The advanced scoring method models each source as an independent probabilistic signal.

**Layer 1 — Source-specific calibration**

Each source gets its own transform appropriate to its data type:
- Open Targets: piecewise linear calibration mapping raw score → P(Tclin | score)
- DepMap: sigmoid transform within the essential gene subpopulation (µ=−1.0, σ=0.30)
- PHAROS: direct tier lookup (Tclin → p=1.00, Tchem → p=0.62, Tbio → p=0.28, Tdark → p=0.05)

**Layer 2 — Noisy-OR fusion**

```
P(druggable) = 1 − ∏ᵢ (1 − wᵢ × pᵢ)
```

Source weights (sum to 1.0):

| Source | Weight | Rationale |
|--------|--------|-----------|
| PHAROS | 0.40 | Direct druggability tier — strongest single signal |
| DepMap | 0.35 | Functional CRISPR evidence — orthogonal to annotations |
| Open Targets | 0.25 | Pre-aggregated, strong disease context |

Two independent moderate signals compound (e.g. 0.6 + 0.6 → 0.84) rather than average to 0.6. Absent sources are omitted from the product — not treated as zero.

**Layer 3 — Shrinkage toward empirical prior**

```
final = α × P(druggable) + (1 − α) × 0.22

α = n_eff / (n_eff + 1.5)
```

The prior (0.22) is the empirical fraction of the human genome that is druggable, from Finan et al. 2017. With n_eff < 0.8, the prior dominates and the disease is flagged as Tier C (suppressed). This prevents single weak signals from producing misleadingly confident scores.

**Tier allocation**

| Tier | n_eff threshold | Action |
|------|----------------|--------|
| A | ≥ 2.0 | Report score — multi-source evidence confirmed |
| B | 0.8 – 2.0 | Report with warning — prior-dominated |
| C | < 0.8 | Suppress — insufficient evidence |

### Reference — Confidence-Weighted Average

```
score(d) = Σ(raw_score × confidence) / Σ(confidence)
```

Simpler and more robust when disease ID alignment is imperfect — it doesn't require records to share the same disease ID to contribute to a score. Used as a baseline comparison and achieves 100% Recall@20 on BRAF.

---

## Disease normalization

The three sources use different disease ontologies:

- Open Targets → EFO IDs (`EFO_0000756`)
- PHAROS → MONDO IDs (`MONDO:0005105`)
- DepMap → free-text lineage names (`"skin neoplasm"`)

`MondoNormalizationAgent` converts all disease inputs to canonical MONDO IDs using route-based OLS4 resolution plus fallback:

1. MONDO passthrough: if the input is already MONDO, keep it (validated when possible).
2. External ID route: for IDs such as EFO/DOID/ORDO-style values, query the source ontology in OLS4 (`queryFields=obo_id`), retrieve the label, then map that label to MONDO.
3. Text route: for free text or unresolved IDs, search MONDO by label/synonym (`queryFields=label,synonym`).
4. Second-pass fallback: if an ID-based attempt still returns `MONDO_UNKNOWN`, retry using `disease_name` text.
5. Targeted overrides: curated mappings for known unresolved Open Targets terms in this dataset (for example, rasopathy and vascular malformation entries).

Every normalization decision is written into metadata (`original_disease_id`, `original_disease_name`, `normalization_match_status`, `normalization_matched_from`) to keep the process auditable.

Normalization APIs used in this prototype:

- OLS4 API docs: https://www.ebi.ac.uk/ols4/api-docs
- OLS4 MCP endpoint: https://www.ebi.ac.uk/ols4/mcp

**Current behavior in this notebook:** Open Targets unresolved terms are reduced to zero after fallback + targeted overrides, while most remaining `MONDO_UNKNOWN` rows come from DepMap lineage/anatomy-style terms that are not direct disease entities.

---

## Benchmark results — BRAF

Ground truth: 4 FDA-approved BRAF drug indications — melanoma, non-small cell lung carcinoma, colorectal cancer, thyroid carcinoma.

| Strategy | Hits | Recall@20 | Role |
|----------|------|-----------|------|
| Uniform average | 3 | 75% | Baseline |
| Confidence-weighted average | 4 | **100%** | Reference |
| Noisy-OR + shrinkage (PRIMARY) | 2 | 50% | Target method |
| Max score | 3 | 75% | Upper bound |

The 50% Noisy-OR result is a known upstream issue, not a model failure. Because disease alignment between Open Targets (EFO IDs) and DepMap (tissue labels) is incomplete, most diseases remain in a single source, keeping n_eff below 0.8 and triggering Tier C suppression. The Noisy-OR model is scientifically correct — it just requires complete alignment to perform. This is the primary motivation for prioritising alignment work in the GSoC plan.

---

## Output

### Console output (top 10 PRIMARY)

```
Disease                                        Score  Tier  N Sources
--------------------------------------------- ------  ----  ---------
melanoma                                      0.4000     B          2
cardiofaciocutaneous syndrome                 0.3988     B          2
colorectal cancer                             0.3972     B          2
cancer                                        0.3970     B          2
LEOPARD syndrome 3                            0.3962     B          2
Noonan syndrome 7                             0.3959     B          2
...
```

### JSON export

Each pipeline run saves `BRAF_agent4target_report.json`:

```json
{
  "target_symbol": "BRAF",
  "n_diseases": 67,
  "n_evidence_records": 185,
  "overall_score": 0.1644,
  "sources": ["open_targets", "depmap", "pharos"],
  "top_evidence_types": ["functional_dependency", "known_drug", "somatic_mutation", ...],
  "per_disease": [
    {
      "disease_name": "melanoma",
      "final_score": 0.4000,
      "tier": "B",
      "n_sources": 2,
      "sources_present": ["open_targets", "pharos"],
      "sources_missing": ["depmap"],
      "jackknife_ci": [0.28, 0.44],
      "breakdown": [
        { "source": "open_targets", "evidence_type": "somatic_mutation",
          "raw_score": 0.9118, "confidence": 0.88, "weighted_contribution": 0.8024 },
        ...
      ]
    }
  ]
}
```

---

## Known limitations

**PHAROS raw_score is flat (0.5) for all tiers.** The calibration step correctly adjusts this via `calibrate_pharos()`, but the raw value in the EvidenceRecord itself is not tier-specific. This is fixed in the improved version by mapping raw_score directly to `TDL_RAW_SCORE_MAP`.

**Residual `MONDO_UNKNOWN` rows are now mostly DepMap-derived lineage/anatomy terms.** Terms like tissue or cell-context labels are often not disease ontology entities, so they do not map cleanly to MONDO without an additional lineage-to-disease translation step. This is now the main normalization bottleneck for multi-source merging.

**Jackknife CI with 2–3 sources has high variance.** The confidence interval is computed by leave-one-source-out — with only 3 sources this gives 3 replicates, producing a wide interval. The CI becomes more meaningful as the pipeline grows to 5 agents.

**DepMap requires a local CSV download.** The file is ~200 MB and is not included in the repository. The pipeline gracefully skips DepMap if the file is absent.

---

## What's coming in GSoC

The prototype uses 3 sources and a simple orchestrator. The full GSoC plan extends this into a 5-agent, evaluation-driven pipeline with stronger normalization and calibration:

- Agent 4 — Literature evidence (PubMed retrieval + LLM-assisted extraction/RAG).
- Agent 5 — Network/pathway evidence (PPI + pathway topology support).
- LangGraph orchestration replacing manual sequential execution (stateful retries, branch control, better observability).
- Normalization expansion:
  - broader external-ID coverage (EFO/DOID/OMIM/ORDO and related xrefs),
  - ontology-aware lineage/tissue-to-disease translation for DepMap terms,
  - improved MONDO merge coverage to reduce Tier C suppression caused by unresolved alignment.
- Score calibration upgrade:
  - empirically fit calibration curves and source weights on held-out Tclin genes,
  - replace heuristic breakpoints with data-driven calibration,
  - retune shrinkage and tier thresholds after calibration.
- Evaluation expansion:
  - beyond Recall@20 to include precision-oriented metrics and ranking diagnostics,
  - ablations (with/without each agent, with/without normalization layers),
  - per-source and per-disease-family error analysis.
- Aggregation outputs:
  - disease-level ranking (current),
  - gene-level druggability roll-up across diseases,
  - richer explanation artifacts for each prediction.
- Productionization:
  - CLI workflow + reproducible config,
  - persistent evidence store/knowledge cache,
  - automated tests and validation checks for each pipeline stage.

### Planned Agent 4 — Literature Evidence Agent (Detailed)

Agent 4 will convert biomedical literature into structured, disease-linked evidence records that can be fused with Open Targets, DepMap, and PHAROS.

What Agent 4 will do:

- Retrieve gene-relevant papers from PubMed and related metadata endpoints.
- Build a retrieval index (RAG) over abstracts/full text where available.
- Extract explicit target-disease statements and classify evidence strength.
- Attach provenance for each extracted claim (PMID, snippet, extraction confidence).
- Normalize disease mentions through the same MONDO normalization pipeline used by other agents.

Planned input/output contract:

- Input:
  - gene symbol and Ensembl ID,
  - optional disease filters and date window.
- Output (`EvidenceRecord`-compatible):
  - `source = literature`,
  - disease ID/name (MONDO-normalized downstream),
  - evidence type (literature-driven association classes),
  - raw score (claim strength / retrieval relevance),
  - confidence (model confidence + evidence redundancy),
  - metadata (PMID list, sentence spans, extraction method).

Scoring integration plan:

- Agent 4 will contribute as an additional calibrated source in Noisy-OR fusion.
- During calibration phase, its weight and calibration curve will be fitted against held-out Tclin targets.
- In low-evidence settings, repeated independent literature support will increase effective support instead of acting as a single anecdotal signal.

Proposal links for Agent 4 stack:

- PubMed: https://pubmed.ncbi.nlm.nih.gov/
- BioMistral: https://huggingface.co/BioMistral

### Planned Agent 5 — Network/Pathway Evidence Agent (Detailed)

Agent 5 will add systems-level support by quantifying how strongly a target is connected to disease-relevant biology in interaction and pathway graphs.

What Agent 5 will do:

- Build a gene-centric neighborhood from curated PPI resources.
- Map neighborhood genes to disease pathways and disease-associated modules.
- Compute network support features (connectivity, pathway overlap, neighborhood enrichment).
- Convert graph features into calibrated evidence records aligned to MONDO diseases.

Planned input/output contract:

- Input:
  - target gene,
  - interaction graph and pathway catalog,
  - optional disease seed gene sets.
- Output (`EvidenceRecord`-compatible):
  - `source = network`,
  - disease ID/name (MONDO-normalized downstream),
  - evidence type (network/pathway support),
  - raw score (network support magnitude),
  - confidence (graph quality + consistency across resources),
  - metadata (top supporting neighbors/pathways, enrichment stats).

Scoring integration plan:

- Agent 5 will be fused as an orthogonal signal to literature and clinical/genetic sources.
- Calibration and weighting will be fitted jointly with other agents to avoid over-crediting hub genes.
- Explanations will include top pathways and interaction-derived justifications for auditability.

Proposal links for Agent 5 stack:

- InWeb (network resource): https://www.intomics.com/inbio/map/
- Reactome: https://reactome.org/

### Proposed 5-Agent Fusion Weights (Initial, To Be Fitted)

Initial allocation discussed for the full pipeline (subject to empirical fitting):

- PHAROS: 0.30
- DepMap: 0.28
- Open Targets + OMIM: 0.22
- Agent 5 (Network): 0.10
- Agent 4 (Literature): 0.10

Final weights will be learned during calibration/evaluation milestones rather than being fixed.

---

## References

1. Raies et al. (2022). DrugnomeAI — the foundational framework this project improves upon. *Communications Biology*. https://doi.org/10.1038/s42003-022-04245-4
2. Open Targets Platform. https://platform.opentargets.org
3. DepMap Cancer Dependency Map. https://depmap.org/portal
4. PHAROS / TCRD. https://pharos.nih.gov
5. MONDO Disease Ontology. https://mondo.monarchinitiative.org
6. Finan et al. (2017). The druggable genome. *Science Translational Medicine*. https://doi.org/10.1126/scitranslmed.aag1166 — source of the 0.22 empirical prior

---

## Repository

**Notebook:** `Agent4Target_Prototype.ipynb`  
**Output:** `BRAF_agent4target_report.json` (generated on run)  
**GitHub:** https://github.com/BayesianDecoder/agent4target-Prototype
