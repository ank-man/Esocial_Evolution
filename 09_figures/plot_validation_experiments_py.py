#!/usr/bin/env python3
"""Create publication-ready SVG figures from validation experiment JSON outputs.

Stdlib-only implementation for environments without R/Matplotlib.
"""
from __future__ import annotations
import json
from pathlib import Path
from statistics import mean

W, H = 1200, 800
M = dict(l=90, r=40, t=70, b=90)
PLOT_W = W - M['l'] - M['r']
PLOT_H = H - M['t'] - M['b']

COLORS = {
    'dark': '#222222',
    'grid': '#DDDDDD',
    'blue': '#1f77b4',
    'orange': '#ff7f0e',
    'green': '#2ca02c',
    'red': '#d62728',
    'purple': '#9467bd',
}


def sx(x, xmin, xmax):
    return M['l'] + (x - xmin) / (xmax - xmin) * PLOT_W if xmax > xmin else M['l']


def sy(y, ymin, ymax):
    return M['t'] + (1 - (y - ymin) / (ymax - ymin)) * PLOT_H if ymax > ymin else M['t'] + PLOT_H / 2


def svg_header(title: str):
    return [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}">',
        '<rect width="100%" height="100%" fill="white"/>',
        f'<text x="{W/2}" y="35" text-anchor="middle" font-size="28" font-family="Arial" fill="{COLORS["dark"]}">{title}</text>',
    ]


def draw_axes(parts, xlabel, ylabel, xmin, xmax, ymin, ymax, xticks, yticks):
    parts.append(f'<line x1="{M["l"]}" y1="{M["t"]+PLOT_H}" x2="{M["l"]+PLOT_W}" y2="{M["t"]+PLOT_H}" stroke="{COLORS["dark"]}"/>')
    parts.append(f'<line x1="{M["l"]}" y1="{M["t"]}" x2="{M["l"]}" y2="{M["t"]+PLOT_H}" stroke="{COLORS["dark"]}"/>')
    for t in xticks:
        x = sx(t, xmin, xmax)
        parts.append(f'<line x1="{x}" y1="{M["t"]+PLOT_H}" x2="{x}" y2="{M["t"]+PLOT_H+6}" stroke="{COLORS["dark"]}"/>')
        parts.append(f'<text x="{x}" y="{M["t"]+PLOT_H+24}" text-anchor="middle" font-size="16" font-family="Arial">{t}</text>')
    for t in yticks:
        y = sy(t, ymin, ymax)
        parts.append(f'<line x1="{M["l"]-6}" y1="{y}" x2="{M["l"]}" y2="{y}" stroke="{COLORS["dark"]}"/>')
        parts.append(f'<line x1="{M["l"]}" y1="{y}" x2="{M["l"]+PLOT_W}" y2="{y}" stroke="{COLORS["grid"]}"/>')
        parts.append(f'<text x="{M["l"]-10}" y="{y+5}" text-anchor="end" font-size="16" font-family="Arial">{t:.2f}</text>')
    parts.append(f'<text x="{M["l"]+PLOT_W/2}" y="{H-25}" text-anchor="middle" font-size="18" font-family="Arial">{xlabel}</text>')
    parts.append(f'<text x="25" y="{M["t"]+PLOT_H/2}" transform="rotate(-90 25 {M["t"]+PLOT_H/2})" text-anchor="middle" font-size="18" font-family="Arial">{ylabel}</text>')


def fig_synteny(data, out):
    filt = [d for d in data if d['n_chr'] == 8 and d['fusions'] == 0]
    groups = {}
    for d in filt:
        k = f"{d['inv_min']}-{d['inv_max']}"
        groups.setdefault(k, []).append((d['inversions'], d['mean_recovery']))
    for k in groups:
        groups[k].sort()

    xmin, xmax = 10, 60
    ymin, ymax = min(d['mean_recovery'] for d in filt)-0.01, max(d['mean_recovery'] for d in filt)+0.01
    parts = svg_header('Figure 1. Synteny recovery under simulated rearrangements')
    draw_axes(parts, 'Number of inversions', 'Mean recovery', xmin, xmax, ymin, ymax, [10,30,60], [round(ymin,2), round((ymin+ymax)/2,2), round(ymax,2)])
    palette=[COLORS['blue'],COLORS['orange'],COLORS['green']]
    for i,(k,pts) in enumerate(sorted(groups.items())):
        color=palette[i%len(palette)]
        path=' '.join(f"{'M' if j==0 else 'L'} {sx(x,xmin,xmax):.1f} {sy(y,ymin,ymax):.1f}" for j,(x,y) in enumerate(pts))
        parts.append(f'<path d="{path}" fill="none" stroke="{color}" stroke-width="3"/>')
        for x,y in pts:
            parts.append(f'<circle cx="{sx(x,xmin,xmax):.1f}" cy="{sy(y,ymin,ymax):.1f}" r="5" fill="{color}"/>')
        parts.append(f'<text x="{M["l"]+PLOT_W-180}" y="{M["t"]+30+i*24}" font-size="16" font-family="Arial" fill="{color}">{k}</text>')
    parts.append('</svg>')
    out.write_text('\n'.join(parts))


def fig_permutation(perm, out):
    vals=[perm['observed_effect'], perm['perm_mean_effect']]
    labels=['Observed effect','Permuted mean']
    colors=[COLORS['purple'], COLORS['orange']]
    ymin, ymax = 0, max(vals)*1.25
    parts=svg_header('Figure 2. Phylogenetic permutation summary')
    draw_axes(parts,'Metric','Effect size',0,3,ymin,ymax,[1,2],[0, ymax/2, ymax])
    barw=180
    for i,v in enumerate(vals, start=1):
        x=sx(i,0,3)
        y=sy(v,ymin,ymax)
        y0=sy(0,ymin,ymax)
        parts.append(f'<rect x="{x-barw/2:.1f}" y="{y:.1f}" width="{barw}" height="{y0-y:.1f}" fill="{colors[i-1]}"/>')
        parts.append(f'<text x="{x}" y="{y-10:.1f}" text-anchor="middle" font-size="16" font-family="Arial">{v:.3f}</text>')
        parts.append(f'<text x="{x}" y="{y0+40:.1f}" text-anchor="middle" font-size="16" font-family="Arial">{labels[i-1]}</text>')
    parts.append(f'<text x="{W/2}" y="{H-55}" text-anchor="middle" font-size="16" font-family="Arial">Empirical two-sided p = {perm["perm_pvalue_two_sided"]:.4g} ({perm["n_permutations"]} permutations)</text>')
    parts.append('</svg>')
    out.write_text('\n'.join(parts))


def fig_gene_family(gf, out):
    neutral=gf['neutral_model']['delta_social_minus_solitary']
    hypo=gf['hypothesis_model']['delta_social_minus_solitary']
    vals=[neutral,hypo]
    labels=['Neutral','Hypothesis']
    ymin,ymax=0,max(vals)*1.25
    parts=svg_header('Figure 3. Gene-family simulation effect sizes')
    draw_axes(parts,'Model','Δ mean family size (social - solitary)',0,3,ymin,ymax,[1,2],[0,ymax/2,ymax])
    for i,v in enumerate(vals, start=1):
        x=sx(i,0,3); y=sy(v,ymin,ymax); y0=sy(0,ymin,ymax)
        c=[COLORS['blue'],COLORS['green']][i-1]
        parts.append(f'<rect x="{x-90:.1f}" y="{y:.1f}" width="180" height="{y0-y:.1f}" fill="{c}"/>')
        parts.append(f'<text x="{x}" y="{y-10:.1f}" text-anchor="middle" font-size="16" font-family="Arial">{v:.3f}</text>')
        parts.append(f'<text x="{x}" y="{y0+40:.1f}" text-anchor="middle" font-size="16" font-family="Arial">{labels[i-1]}</text>')
    parts.append('</svg>')
    out.write_text('\n'.join(parts))


def fig_hamilton(ham, out):
    ham=sorted(ham, key=lambda d:d['recombination_rate'])
    xs=[d['recombination_rate'] for d in ham]
    ys=[d['final_freq'] for d in ham]
    xmin,xmax=min(xs),max(xs)
    ymin,ymax=0,max(ys)*1.15 if max(ys)>0 else 1
    parts=svg_header('Figure 4. Hamilton-style linkage simulation')
    draw_axes(parts,'Recombination rate','Final social-allele frequency',xmin,xmax,ymin,ymax,xs,[0,ymax/2,ymax])
    path=' '.join(f"{'M' if i==0 else 'L'} {sx(x,xmin,xmax):.1f} {sy(y,ymin,ymax):.1f}" for i,(x,y) in enumerate(zip(xs,ys)))
    parts.append(f'<path d="{path}" fill="none" stroke="{COLORS["red"]}" stroke-width="3"/>')
    for d in ham:
        x,y=d['recombination_rate'],d['final_freq']
        parts.append(f'<circle cx="{sx(x,xmin,xmax):.1f}" cy="{sy(y,ymin,ymax):.1f}" r="5" fill="{COLORS["red"]}"/>')
        parts.append(f'<text x="{sx(x,xmin,xmax):.1f}" y="{sy(y,ymin,ymax)-10:.1f}" text-anchor="middle" font-size="13" font-family="Arial">m={d["hamilton_margin"]:.3f}</text>')
    parts.append('</svg>')
    out.write_text('\n'.join(parts))


def main():
    import argparse
    ap=argparse.ArgumentParser()
    ap.add_argument('in_dir', nargs='?', default='results/validation_experiments')
    ap.add_argument('out_dir', nargs='?', default='results/validation_experiments/figures')
    args=ap.parse_args()
    in_dir=Path(args.in_dir)
    out_dir=Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    sy=json.loads((in_dir/'01_synteny_null.json').read_text())
    perm=json.loads((in_dir/'02_phylo_permutation.json').read_text())
    gf=json.loads((in_dir/'03_gene_family_simulations.json').read_text())
    ham=json.loads((in_dir/'04_hamilton_synteny_simulation.json').read_text())

    fig_synteny(sy, out_dir/'figure_synteny_null_recovery.svg')
    fig_permutation(perm, out_dir/'figure_phylo_permutation_summary.svg')
    fig_gene_family(gf, out_dir/'figure_gene_family_delta.svg')
    fig_hamilton(ham, out_dir/'figure_hamilton_linkage.svg')

    conclusion = f"""# Simulation conclusions

1. **Permutation test supports a sociality-associated effect** in this simulated setup.
   - Observed social-solitary effect: {perm['observed_effect']:.3f}
   - Two-sided empirical p-value: {perm['perm_pvalue_two_sided']:.4g}

2. **Gene-family simulations are sensitive to the modeled hypothesis.**
   - Neutral model Δ(social-solitary): {gf['neutral_model']['delta_social_minus_solitary']:.3f}
   - Hypothesis model Δ(social-solitary): {gf['hypothesis_model']['delta_social_minus_solitary']:.3f}
   - Interpretation: reduced loss in social lineages produces a clearly larger divergence.

3. **Synteny null recovery declines with increasing inversion load** (in this toy setup), while inversion-size/fusion effects are comparatively minor across the tested grid.

4. **Hamilton-style simulation is stochastic at single-replicate scale.**
   Final frequencies do not monotonically track recombination in one run; robust inference needs replicate runs per parameter and confidence intervals.
"""
    (out_dir/'simulation_conclusion.md').write_text(conclusion)

    # compact csv summary
    csv_lines=[
        'metric,value',
        f"observed_effect,{perm['observed_effect']}",
        f"perm_pvalue_two_sided,{perm['perm_pvalue_two_sided']}",
        f"delta_neutral,{gf['neutral_model']['delta_social_minus_solitary']}",
        f"delta_hypothesis,{gf['hypothesis_model']['delta_social_minus_solitary']}"
    ]
    (out_dir/'validation_figure_summary_table.csv').write_text('\n'.join(csv_lines)+'\n')
    print(f'Wrote figures and conclusion to: {out_dir}')


if __name__ == '__main__':
    main()
