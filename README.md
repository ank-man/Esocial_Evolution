# 🐝 Scripts for “Eusociality and Genome Architecture in Hymenoptera”

This repository contains the **reproducible pipeline** and all helper utilities used in the Materials & Methods section of our study on how eusocial behaviour shapes genome organisation across Hymenoptera.

---

## 📁 Directory layout

| Step                          | Folder / key script    | Purpose                                         | Main tools                     |
| ----------------------------- | ---------------------- | ----------------------------------------------- | ------------------------------ |
| **0 Setup**                   | `00_env/conda‑env.yml` | Re‑creates software stack                       | *mamba*, *conda‑lock*          |
| **1 Data harvest**            | `01_genomes/`          | Download & standardise 194 genomes              | *wget*, *Biopython*            |
| **2 Orthology**               | `02_orthology/`        | FastOMA + BUSCO single‑copy extraction          | *FastOMA v2.2*, *BUSCO v5.7*   |
| **3 Phylogenomics**           | `03_species_tree/`     | MAFFT + trimAl + IQ‑TREE supermatrix            | *MAFFT*, *trimAl*, *IQ‑TREE 2* |
| **4 Synteny**                 | `04_synteny/`          | MCScanX on 10 k genome pairs, colinearity stats | *MCScanX*                      |
| **5 Phylogenetic regression** | `05_stats/`            | PGLS / D‑PGLS models, 10 000 perms              | *R nlme*, *caper*              |
| **6 Duplication analysis**    | `06_duplication/`      | Dup/Ret ratios from FastOMA output              | *pandas*                       |
| **7 Enrichment**              | `07_enrichment/`       | Batch ShinyGO queries                           | *R httr*                       |
| **8 Selection tests**         | `08_paml/`             | PAML site‑models (M1a vs M2a)                   | *PAML 4*                       |
| **9 Figures**                 | `09_figures/`          | Jupyter / Rmd notebooks for all plots           | *matplotlib*, *ggplot2*        |

> **Tip :** Every script is numbered in execution order and wrapped in a [Snakemake](https://snakemake.readthedocs.io) workflow in `workflow/`.

---

## 🚀 Quick start

```bash
# clone & enter
$ git clone https://github.com/your‑org/hymenoptera‑eusociality‑scripts.git
$ cd hymenoptera‑eusociality‑scripts/scripts

# create the locked environment (≈ 5 min)
$ mamba env create -f 00_env/conda-env.yml
$ conda activate eusociality-mm

# fetch genomes (≈ 120 GB)
$ bash 01_genomes/download_genomes.sh

# launch the full pipeline (edit -j to CPU count)
$ snakemake -s workflow/Snakefile -j 32 --rerun-incomplete
```

All heavy output is written to `results/` (ignored by Git).

---

## 🗄️ Data & storage

* Raw genome FASTA/GFF files (> 100 GB) are **not** tracked; they are pulled on‑demand.
* Random seed fixed at `SEED=20250708` for deterministic alignments & trees.
* A **toy dataset** with 5 species lives in `toy_data/` for CI and tutorials.

---

## 🖥️ Running on HPC (SLURM)

See `cluster/slurm_jobscript.sh` for a template. Heavy rules iteratively parallelise over species/loci via the Snakemake `--cluster` profile.

---

## 📜 Citing this work

Please cite the companion paper **“Eusociality and the Tempo of Genome Rearrangement in Hymenoptera”** and the individual software listed in `CITATION.cff`.

---

## 🪪 License

Code is released under the **MIT License**. Third‑party software keeps its original license.

---

## 🙋 Questions & contact

Open an issue or email *[ankush.sharma@uga.edu](mailto:ankush.sharma@uga.edu)*.
