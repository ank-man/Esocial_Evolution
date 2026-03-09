# 🐝 Scripts for “Eusociality and Genome Architecture in Hymenoptera”

This repository contains the **exact standalone scripts** referenced in the Materials & Methods section of our study on how eusocial behaviour shapes genome organisation across Hymenoptera.  Each script can be executed directly from the command line—no workflow manager or Conda environment required.

---

## 📁 Directory layout

| Step                          | Folder / key script | Purpose                                         | Main tools                     |
| ----------------------------- | ------------------- | ----------------------------------------------- | ------------------------------ |
| **1 Data harvest**            | `01_genomes/`       | Download & standardise 126 genomes              | *wget*, *Biopython*            |
| **2 Orthology**               | `02_orthology/`     | FastOMA + BUSCO single‑copy extraction          | *FastOMA v2.2*, *BUSCO v5.7*   |
| **3 Phylogenomics**           | `03_species_tree/`  | MAFFT + trimAl + IQ‑TREE supermatrix            | *MAFFT*, *trimAl*, *IQ‑TREE 2* |
| **4 Synteny**                 | `04_synteny/`       | MCScanX on 10 k genome pairs, colinearity stats | *MCScanX*                      |
| **5 Phylogenetic regression** | `05_stats/`         | PGLS / D‑PGLS models, 10 000 perms              | *R nlme*, *caper*              |
| **6 Duplication analysis**    | `06_duplication/`   | Dup/Ret ratios from FastOMA output              | *pandas*                       |
| **7 Enrichment**              | `07_enrichment/`    | Batch ShinyGO queries                           | *R httr*                       |
| **8 Selection tests**         | `08_paml/`          | PAML site‑models (M1a vs M2a)                   | *PAML 4*                       |
| **9 Figures**                 | `09_figures/`       | Jupyter / Rmd notebooks for all plots           | *matplotlib*, *ggplot2*        |

> **Tip :** Scripts are **numbered** in execution order. Run them sequentially or cherry‑pick individual steps.

---

## 🚀 Quick start

1. **Clone the repository**

   ```bash
   git clone https://github.com/ank-man/Esocial_Evolution.git
   cd hymenoptera‑eusociality‑scripts/scripts
   ```
2. **Install prerequisites**   
   Make sure the following tools are in your `$PATH` (versions in parenthesis are those we used):

   * wget, curl
   * Python ≥ 3.10 with Biopython, pandas, numpy
   * R ≥ 4.2 with packages **nlme**, **caper**, **httr**, **ggplot2**
   * FastOMA 2.2, BUSCO 5.7, MAFFT 7.475, trimAl 1.4, IQ‑TREE 2.2, MCScanX, PAML 4.10
3. **Fetch the genomes** (≈ 120 GB)

   ```bash
   bash 01_genomes/download_genomes.sh
   ```
4. **Run the remaining steps** in order, e.g.

   ```bash
   bash 02_orthology/run_fastoma.sh
   python 02_orthology/extract_busco_singletons.py
   bash 03_species_tree/align_mafft.py
   # ...continue following folder numbering
   ```

---


## 🧪 Simulation and validation experiments

To run the full set of simulation/validation experiments (synteny nulls, phylogenetic permutations, gene-family evolution simulation, and Hamilton-style linkage simulation):

```bash
python run_validation_experiments.py
```

Outputs are written to:

- `results/validation_experiments/01_synteny_null.json`
- `results/validation_experiments/02_phylo_permutation.json`
- `results/validation_experiments/03_gene_family_simulations.json`
- `results/validation_experiments/04_hamilton_synteny_simulation.json`

Use `--seed` and `--outdir` to customize reproducibility and output location.


For publication-quality figures from these outputs:

```bash
Rscript 09_figures/plot_validation_experiments.R results/validation_experiments results/validation_experiments/figures
```


---

## 🗄️ Data & storage

* Raw genome FASTA/GFF files (> 100 GB) are **not** tracked; the download script pulls them on‑demand.
* Random seed fixed at `SEED=20250708` for deterministic alignments & trees.
* A **toy dataset** with 5 species lives in `toy_data/` for CI and tutorials.

---

## 🖥️ Running on HPC (SLURM)

Each heavy step (FastOMA, MCScanX, IQ‑TREE, PAML) can be submitted to a cluster. Example wrappers live in `cluster/`; edit CPU, memory, and partition flags, then `sbatch`.

---

## 📜 Citing this work

Please cite the companion paper **“Eusociality and the Tempo of Genome Rearrangement in Hymenoptera”** and the individual software listed in `CITATION.cff`.

---

## 🪪 License

All custom code is released under the **MIT License**. Third‑party software retains its original license.

---

## 🙋 Questions & contact

Open an issue or email *[ankush.sharma@uga.edu](mailto:ankush.sharma@uga.edu)*.
