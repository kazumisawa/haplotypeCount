
# haplo_phylo_v11 — Pivot-based haplotype phylogeny builder

`haplo_phylo_v11.py` is a command‑line tool to **infer derived alleles and build a haplotype phylogeny** from a **bgzip+tabix-indexed VCF (`.vcf.gz`)**.  
It implements a **pivot-based ancestor/descendant inference** (two-pass) and outputs:

- **Haplotype frequency table** (with first-owner naming, variant set, window, and missing-policy label)
- **Per-sample haplotype assignment** (phase order preserved)
- **Filtered sites report** (optional)
- **Newick tree** over **unique haplotypes only** (wavelet-tree style)
- **Step-wise timing report** (optional)

> **Key changes in v11**  
> - Haplotype IDs use the **first owner**: `<sample_id>_0` and `<sample_id>_1` (phase order), reused for identical patterns found later.  
> - `hapcount_by_sample.tsv` uses **phase order** (not lexicographic).  
> - Phylogeny is built over **unique haplotypes** (no duplicates in the tree).  
> - Optional **`--time-report`** generates `{out}.timings.tsv`.

---

## Table of contents

1. Requirements  
2. [nstallation  
3. Input assumptions  
4. Quick start  
5. [ommand-line options  
6. Outputs  
7. Algorithm (pivot-based derived inference)  
8. [avelet-tree phylogeny  
9. Performance & timing  
10. Reproducibility & determinism  
11. [imitations & edge cases  
12. Troubleshooting  
13. License & citation

---

## Requirements

- **Python**: 3.8+ recommended
- **Packages**:
  - `cyvcf2` — VCF reader (must be pre-installed)
- **Input data**:
  - **BGZF-compressed VCF**: `*.vcf.gz`
  - **tabix index**: `*.vcf.gz.tbi` present

> The tool does **not** install packages for you. Ensure `cyvcf2` is available in your environment before running.

---

## Installation

Clone your repository or download `haplo_phylo_v11.py`, then (if needed) install `cyvcf2`:

```bash
pip install cyvcf2
