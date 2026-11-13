"""
LIONESS Project - Step 3: Analyse the sample-specific networks
==============================================================

Research question
    Do pancreatic tumours show different gene-gene co-expression (network
    wiring) from matched normal pancreas, and can LIONESS sample-specific
    networks detect this?

Analyses (all with the 1000 most-variable genes; script 04 varies this)
    A. Differential edges: paired test of tumour vs matched-normal edge weights
       (Benjamini-Hochberg FDR 5%) + label-swap permutation null.
    B. Is the "rewiring" just differential expression?  We compare edge
       changes with what you would expect if two genes were only shifted up or
       down together.
    C. Hub genes of the differential network.
    D. Do the sample networks separate tumour from normal?  (PCA + patient-
       grouped cross-validation, compared with using expression directly.)

Outputs: results/*.png and results/*.csv (small, committed to the repo).
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import GroupKFold, cross_val_score
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

from lioness_utils import (RESULTS_DIR, benjamini_hochberg, differential_edges,
                           load_data, paired_differences, paired_t_stats,
                           upper_indices)

N_GENES = 1000
TAG = f"{N_GENES}genes"


def main():
    expr, meta = load_data()
    genes = (RESULTS_DIR / f"genes_{TAG}.txt").read_text().split("\n")
    nets = np.load(RESULTS_DIR / f"lioness_networks_{TAG}.npy")
    sub = expr.loc[genes]
    iu = upper_indices(len(genes))
    y = (meta["condition"] == "tumor").astype(int).values
    groups = meta["patient"].values
    lines = []  # text collected for results/summary_main_analysis.txt

    def log(msg=""):
        print(msg)
        lines.append(msg)

    # ------------------------------------------------------------------ A
    log("=== A. Differential edges (paired, tumour vs matched normal) ===")
    res = differential_edges(nets, meta, n_perm=200)
    n_edges = nets.shape[1]
    up = int((res["sig"] & (res["mean_diff"] > 0)).sum())
    down = int((res["sig"] & (res["mean_diff"] < 0)).sum())
    log(f"Edges tested: {n_edges}")
    log(f"Significant at FDR 5%: {res['n_sig']} ({100 * res['n_sig'] / n_edges:.1f}%)"
        f"  [stronger in tumour: {up}, weaker in tumour: {down}]")
    log(f"Label-swap null: {res['null_counts'].mean():.1f} edges pass the same "
        f"|t|>={res['crit']:.2f} cut-off on average "
        f"(range {res['null_counts'].min()}-{res['null_counts'].max()})")
    nc = res["null_counts"]
    log(f"  median {np.median(nc):.0f}, 95th percentile {np.percentile(nc, 95):.0f} "
        "(very spread out: neighbouring edges share genes, so they are not independent)")
    log(f"Empirical FDR (mean null count / observed) = {res['emp_fdr']:.4f}; "
        f"using the 95th percentile instead: {np.percentile(nc, 95) / res['n_sig']:.4f}")

    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.hist(res["t"], bins=100, color="steelblue", alpha=0.8, log=True)
    ax.axvline(res["crit"], color="red", ls="--", label="FDR 5% cut-off")
    ax.axvline(-res["crit"], color="red", ls="--")
    ax.set_xlabel("paired t-statistic (tumour - normal edge weight)")
    ax.set_ylabel("number of edges (log scale)")
    ax.set_title("Edge-level differences between tumour and matched normal")
    ax.legend()
    fig.tight_layout()
    fig.savefig(RESULTS_DIR / "fig1_edge_tstat_distribution.png", dpi=150)
    plt.close(fig)

    # ------------------------------------------------------------------ B
    log("\n=== B. Is the rewiring just differential expression? ===")
    Dg, _, _ = paired_differences(sub.values.T, meta)  # patients x genes
    gene_t = paired_t_stats(Dg)
    gene_q = benjamini_hochberg(2 * stats.t.sf(np.abs(gene_t), df=Dg.shape[0] - 1))
    gene_de = gene_q <= 0.05
    log(f"Genes differentially expressed (paired t, FDR 5%): "
        f"{int(gene_de.sum())} of {len(genes)}")
    both_de = gene_de[iu[0]] & gene_de[iu[1]]
    log(f"Significant edges where BOTH genes are DE: "
        f"{100 * both_de[res['sig']].mean():.1f}%  "
        f"(all edges: {100 * both_de.mean():.1f}%)")
    # If two genes simply move up/down together between conditions, their joint
    # correlation changes in proportion to the product of their shifts.
    expected = gene_t[iu[0]] * gene_t[iu[1]]
    rho, _ = stats.spearmanr(res["mean_diff"], expected)
    log(f"Spearman correlation between mean edge change and the product of the "
        f"two genes' DE t-statistics: rho = {rho:.3f}")

    rng = np.random.default_rng(0)
    pick = rng.choice(n_edges, 20000, replace=False)
    fig, ax = plt.subplots(figsize=(6, 5))
    ax.scatter(expected[pick], res["mean_diff"][pick], s=3, alpha=0.3)
    ax.set_xlabel("product of the two genes' DE t-statistics")
    ax.set_ylabel("mean LIONESS edge change (tumour - normal)")
    ax.set_title(f"Edge change vs differential expression (rho = {rho:.2f})")
    fig.tight_layout()
    fig.savefig(RESULTS_DIR / "fig2_edge_change_vs_expression.png", dpi=150)
    plt.close(fig)

    # ------------------------------------------------------------------ C
    log("\n=== C. Hub genes of the differential network ===")
    degree = np.bincount(np.r_[iu[0][res["sig"]], iu[1][res["sig"]]],
                         minlength=len(genes))
    hubs = pd.DataFrame({"gene": genes, "n_significant_edges": degree,
                         "expression_t": gene_t, "DE_fdr5": gene_de})
    hubs = hubs.sort_values("n_significant_edges", ascending=False)
    hubs.head(25).to_csv(RESULTS_DIR / "top_hub_genes.csv", index=False)
    log(hubs.head(10).to_string(index=False))
    log(f"Top-100 hubs that are also DE: {hubs.head(100)['DE_fdr5'].mean() * 100:.0f}%")

    sig_idx = np.where(res["sig"])[0]
    order = sig_idx[np.argsort(-np.abs(res["t"][sig_idx]))][:5000]
    pd.DataFrame({"gene1": np.array(genes)[iu[0][order]],
                  "gene2": np.array(genes)[iu[1][order]],
                  "mean_edge_change": res["mean_diff"][order],
                  "t": res["t"][order], "q": res["q"][order]}
                 ).to_csv(RESULTS_DIR / "top_differential_edges.csv", index=False)

    # ------------------------------------------------------------------ D
    log("\n=== D. Do sample networks separate tumour from normal? ===")
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.5))
    for ax, data, name in [(axes[0], sub.values.T, "gene expression"),
                           (axes[1], nets, "LIONESS edge weights")]:
        pc = PCA(n_components=2).fit(data)
        Zp = pc.transform(data)
        for lab, col in [(0, "tab:blue"), (1, "tab:red")]:
            ax.scatter(Zp[y == lab, 0], Zp[y == lab, 1], c=col,
                       label=["normal", "tumour"][lab], s=25)
        ax.set_xlabel(f"PC1 ({100 * pc.explained_variance_ratio_[0]:.0f}%)")
        ax.set_ylabel(f"PC2 ({100 * pc.explained_variance_ratio_[1]:.0f}%)")
        ax.set_title(f"PCA of {name}")
        ax.legend()
    fig.tight_layout()
    fig.savefig(RESULTS_DIR / "fig3_pca_expression_vs_networks.png", dpi=150)
    plt.close(fig)

    # How spread out are the networks inside each group? (mean absolute edge
    # weight per sample; large values = sample strongly different from average)
    strength = np.abs(nets).mean(axis=1)
    log(f"Mean |edge weight| per sample: normal median {np.median(strength[y == 0]):.2f} "
        f"(max {strength[y == 0].max():.2f}), tumour median {np.median(strength[y == 1]):.2f} "
        f"(max {strength[y == 1].max():.2f})")

    # Patient-grouped cross-validation: the tumour and normal sample of the
    # same patient are never split between training and test data.
    cv = GroupKFold(n_splits=6)
    clf = make_pipeline(StandardScaler(), LogisticRegression(C=0.01, max_iter=2000))
    acc_expr = cross_val_score(clf, sub.values.T, y, groups=groups, cv=cv).mean()
    acc_net = cross_val_score(clf, nets, y, groups=groups, cv=cv).mean()
    y_perm = np.random.default_rng(1).permutation(y)
    acc_ctrl = cross_val_score(clf, nets, y_perm, groups=groups, cv=cv).mean()
    log("Cross-validated accuracy (6-fold, grouped by patient):")
    log(f"  gene expression (1000 genes):     {acc_expr:.3f}")
    log(f"  LIONESS edge weights (499500):    {acc_net:.3f}")
    log(f"  control, shuffled labels (edges): {acc_ctrl:.3f}")

    (RESULTS_DIR / "summary_main_analysis.txt").write_text("\n".join(lines))
    print("\nSaved figures and tables to results/")


if __name__ == "__main__":
    main()
