#!/usr/bin/env python3
"""Run simulation and validation experiments described for Hymenoptera analyses.

This script implements lightweight, reproducible approximations of:
1) Null simulations for synteny statistics (inversions + chromosome fusions).
2) Phylogenetic-style permutation tests for sociality labels.
3) Neutral/hypothesis-driven gene family evolution simulations.
4) A synteny-linked Hamilton's-rule forward simulation.

Outputs are written under ./results/validation_experiments.
"""
from __future__ import annotations

import argparse
import json
import math
import random
from dataclasses import dataclass
from pathlib import Path
from statistics import mean
from typing import Dict, List, Sequence, Tuple


SEED = 20250708


def set_seed(seed: int) -> None:
    random.seed(seed)


# -----------------------
# 1) Synteny null simulation
# -----------------------

def build_genome(n_chr: int, genes_per_chr: int) -> List[List[int]]:
    gid = 1
    genome: List[List[int]] = []
    for _ in range(n_chr):
        chrom = list(range(gid, gid + genes_per_chr))
        gid += genes_per_chr
        genome.append(chrom)
    return genome


def apply_inversion(genome: List[List[int]], min_len: int, max_len: int) -> None:
    non_empty = [i for i, c in enumerate(genome) if len(c) >= min_len]
    if not non_empty:
        return
    ci = random.choice(non_empty)
    chrom = genome[ci]
    inv_len = random.randint(min_len, min(max_len, len(chrom)))
    start = random.randint(0, len(chrom) - inv_len)
    chrom[start : start + inv_len] = reversed(chrom[start : start + inv_len])


def apply_fusion(genome: List[List[int]]) -> None:
    non_empty = [i for i, c in enumerate(genome) if c]
    if len(non_empty) < 2:
        return
    a, b = sorted(random.sample(non_empty, 2))
    genome[a] = genome[a] + genome[b]
    genome[b] = []


def adjacency_set(genome: Sequence[Sequence[int]]) -> set[Tuple[int, int]]:
    pairs = set()
    for chrom in genome:
        for i in range(len(chrom) - 1):
            x, y = chrom[i], chrom[i + 1]
            pairs.add((min(x, y), max(x, y)))
    return pairs


def synteny_recovery_score(ref: Sequence[Sequence[int]], test: Sequence[Sequence[int]]) -> float:
    a = adjacency_set(ref)
    b = adjacency_set(test)
    if not a:
        return 0.0
    return len(a & b) / len(a)


def simulate_synteny_null(
    n_reps: int = 200,
    n_chr_values: Sequence[int] = (8, 16),
    genes_per_chr: int = 120,
    inversion_counts: Sequence[int] = (10, 30, 60),
    inversion_sizes: Sequence[Tuple[int, int]] = ((5, 15), (20, 40), (50, 80)),
    fusion_counts: Sequence[int] = (0, 2, 4),
) -> List[Dict[str, float]]:
    out: List[Dict[str, float]] = []
    for n_chr in n_chr_values:
        for inv_n in inversion_counts:
            for inv_min, inv_max in inversion_sizes:
                for fus_n in fusion_counts:
                    scores = []
                    for _ in range(n_reps):
                        ref = build_genome(n_chr=n_chr, genes_per_chr=genes_per_chr)
                        test = [c[:] for c in ref]
                        for _ in range(inv_n):
                            apply_inversion(test, inv_min, inv_max)
                        for _ in range(fus_n):
                            apply_fusion(test)
                        scores.append(synteny_recovery_score(ref, test))
                    out.append(
                        {
                            "n_chr": n_chr,
                            "inversions": inv_n,
                            "inv_min": inv_min,
                            "inv_max": inv_max,
                            "fusions": fus_n,
                            "mean_recovery": mean(scores),
                            "min_recovery": min(scores),
                            "max_recovery": max(scores),
                        }
                    )
    return out


# -----------------------
# 2) Permutation test (with tree-like constrained shuffling)
# -----------------------

@dataclass
class SpeciesDatum:
    species: str
    clade: str
    social: int  # 1=social, 0=solitary
    block_size: float


def synthesize_species_data() -> List[SpeciesDatum]:
    """Toy dataset with clade structure and moderate sociality effect."""
    clades = {
        "Ants": ["A1", "A2", "A3", "A4", "A5", "A6"],
        "Bees": ["B1", "B2", "B3", "B4", "B5", "B6"],
        "Wasps": ["W1", "W2", "W3", "W4", "W5", "W6"],
    }
    social_species = {"A1", "A2", "A3", "A4", "B1", "B2", "B3", "W1", "W2"}
    data: List[SpeciesDatum] = []
    for clade, sps in clades.items():
        clade_shift = {"Ants": 0.08, "Bees": 0.03, "Wasps": -0.04}[clade]
        for sp in sps:
            social = 1 if sp in social_species else 0
            noise = random.uniform(-0.05, 0.05)
            block_size = 0.55 + clade_shift + social * 0.12 + noise
            data.append(SpeciesDatum(species=sp, clade=clade, social=social, block_size=block_size))
    return data


def social_effect(data: Sequence[SpeciesDatum]) -> float:
    social_vals = [d.block_size for d in data if d.social == 1]
    solitary_vals = [d.block_size for d in data if d.social == 0]
    return mean(social_vals) - mean(solitary_vals)


def phylo_permutation_pvalue(data: Sequence[SpeciesDatum], n_perm: int = 5000) -> Dict[str, float]:
    observed = social_effect(data)
    clade_to_labels: Dict[str, List[int]] = {}
    by_clade: Dict[str, List[SpeciesDatum]] = {}
    for d in data:
        clade_to_labels.setdefault(d.clade, []).append(d.social)
        by_clade.setdefault(d.clade, []).append(d)

    extreme = 0
    permuted_effects = []
    for _ in range(n_perm):
        perm_data: List[SpeciesDatum] = []
        for clade, members in by_clade.items():
            labels = clade_to_labels[clade][:]
            random.shuffle(labels)
            for m, lbl in zip(members, labels):
                perm_data.append(SpeciesDatum(m.species, m.clade, lbl, m.block_size))
        eff = social_effect(perm_data)
        permuted_effects.append(eff)
        if abs(eff) >= abs(observed):
            extreme += 1

    pval = (extreme + 1) / (n_perm + 1)
    return {
        "observed_effect": observed,
        "perm_pvalue_two_sided": pval,
        "perm_mean_effect": mean(permuted_effects),
        "n_permutations": n_perm,
    }


# -----------------------
# 3) Gene family evolution simulation
# -----------------------

def simulate_gene_families(
    n_families: int = 800,
    n_steps: int = 30,
    p_dup: float = 0.09,
    p_loss: float = 0.07,
    eusocial_loss_modifier: float = 1.0,
) -> Dict[str, float]:
    """Simple birth/death walk across social vs solitary branches."""
    # 10 eusocial and 10 solitary species sampled at tips.
    eusocial_sizes = []
    solitary_sizes = []

    for _ in range(n_families):
        start = random.randint(5, 20)
        size_social = start
        size_solitary = start
        for _ in range(n_steps):
            # social trajectory
            if random.random() < p_dup:
                size_social += 1
            if size_social > 1 and random.random() < p_loss * eusocial_loss_modifier:
                size_social -= 1
            # solitary trajectory
            if random.random() < p_dup:
                size_solitary += 1
            if size_solitary > 1 and random.random() < p_loss:
                size_solitary -= 1
        eusocial_sizes.append(size_social)
        solitary_sizes.append(size_solitary)

    return {
        "n_families": n_families,
        "mean_family_size_eusocial": mean(eusocial_sizes),
        "mean_family_size_solitary": mean(solitary_sizes),
        "delta_social_minus_solitary": mean(eusocial_sizes) - mean(solitary_sizes),
        "p_dup": p_dup,
        "p_loss_base": p_loss,
        "eusocial_loss_modifier": eusocial_loss_modifier,
    }


# -----------------------
# 4) Synteny-linked Hamilton simulation
# -----------------------

def logistic(x: float) -> float:
    return 1.0 / (1.0 + math.exp(-x))


def simulate_social_allele(
    generations: int = 120,
    pop_size: int = 1500,
    init_freq: float = 0.05,
    b: float = 0.10,
    c: float = 0.03,
    relatedness: float = 0.4,
    recombination_rate: float = 0.02,
) -> Dict[str, float]:
    """One-locus approximation where effective relatedness declines with recombination."""
    p = init_freq
    history = [p]
    for _ in range(generations):
        r_syn = relatedness * (1.0 - recombination_rate)
        selection_score = b * r_syn - c
        w_a = 1.0 + 0.6 * selection_score
        w_A = 1.0
        # viability selection
        p_sel = (p * w_a) / (p * w_a + (1 - p) * w_A)
        # drift (binomial via Bernoulli sum to avoid numpy dependency)
        count = 0
        for _ in range(pop_size):
            if random.random() < p_sel:
                count += 1
        p = count / pop_size
        history.append(p)

    return {
        "final_freq": history[-1],
        "max_freq": max(history),
        "min_freq": min(history),
        "mean_freq": mean(history),
        "generations": generations,
        "recombination_rate": recombination_rate,
        "hamilton_margin": b * relatedness * (1.0 - recombination_rate) - c,
    }


def run_hamilton_grid() -> List[Dict[str, float]]:
    out = []
    for rr in (0.01, 0.05, 0.15, 0.30):
        summary = simulate_social_allele(recombination_rate=rr)
        out.append(summary)
    return out


def write_json(path: Path, payload) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2))


def main() -> None:
    parser = argparse.ArgumentParser(description="Run all simulation/validation experiments")
    parser.add_argument("--seed", type=int, default=SEED)
    parser.add_argument("--outdir", type=Path, default=Path("results/validation_experiments"))
    args = parser.parse_args()

    set_seed(args.seed)
    outdir: Path = args.outdir
    outdir.mkdir(parents=True, exist_ok=True)

    synteny = simulate_synteny_null()
    write_json(outdir / "01_synteny_null.json", synteny)

    species_data = synthesize_species_data()
    perm = phylo_permutation_pvalue(species_data)
    write_json(outdir / "02_phylo_permutation.json", perm)

    neutral = simulate_gene_families(eusocial_loss_modifier=1.0)
    hypothesis = simulate_gene_families(eusocial_loss_modifier=0.7)
    write_json(
        outdir / "03_gene_family_simulations.json",
        {"neutral_model": neutral, "hypothesis_model": hypothesis},
    )

    hamilton = run_hamilton_grid()
    write_json(outdir / "04_hamilton_synteny_simulation.json", hamilton)

    summary = {
        "seed": args.seed,
        "outputs": [
            "01_synteny_null.json",
            "02_phylo_permutation.json",
            "03_gene_family_simulations.json",
            "04_hamilton_synteny_simulation.json",
        ],
    }
    write_json(outdir / "summary.json", summary)

    print(f"Wrote validation outputs to: {outdir}")


if __name__ == "__main__":
    main()
