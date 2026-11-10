"""
LIONESS Project - Step 4: Modifications beyond the tutorial
===========================================================

The tutorial only builds the networks. Here we add three things and evaluate
each one:

  1. PARAMETER SWEEP (number of genes). Pearson LIONESS weights of a gene pair
     depend only on those two genes, so changing the number of genes does not
     change any edge weight - it only changes how many edges are tested. We
     check what that does to the number of significant edges.

  2. PARAMETER SWEEP (cohort size and edge-selection rule). LIONESS needs
     enough samples. We repeatedly draw random subsets of 9, 18 or 27 patients
     (from the 36) and ask how well the top edges agree with the top edges
     found with all patients (Jaccard index). We also vary the sparsity
     cut-off (top 0.1%, 1%, 5% of edges ranked by |t|).

  3. METHOD COMPARISON with the classical alternative taught in network
     analysis: group-wise differential co-expression
         (Pearson r in tumour samples) - (Pearson r in normal samples),
     tested with Fisher's z and a label-swap permutation null.
     Same data, same genes, same 'top-k edges' rule, same stability test.

Outputs: results/fig4_*.png, results/fig5_*.png, results/*.csv, results/summary_sweep_comparison.txt
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

from lioness_utils import (RESULTS_DIR, benjamini_hochberg, differential_edges,
                           lioness_pearson, load_data, select_variable_genes,
                           upper_indices)

rng = np.random.default_rng(2024)
lines = []


def log(msg=""):
    print(msg)
    lines.append(msg)


# ---------------------------------------------------------------------------
# Classical differential co-expression (comparison method)
# ---------------------------------------------------------------------------
def classic_diff_coexpression(X, tumor_idx, normal_idx, iu):
    """Fisher-z difference of Pearson r between tumour and normal samples."""
    r_t = np.corrcoef(X[:, tumor_idx])[iu]
    r_n = np.corrcoef(X[:, normal_idx])[iu]
    z_t = np.arctanh(np.clip(r_t, -0.999999, 0.999999))
    z_n = np.arctanh(np.clip(r_n, -0.999999, 0.999999))
    se = np.sqrt(1.0 / (len(tumor_idx) - 3) + 1.0 / (len(normal_idx) - 3))
    return (z_t - z_n) / se, r_t - r_n


def classic_with_null(X, meta, iu, n_perm=100, fdr=0.05):
    """Classic differential co-expression + within-pair label-swap null."""
    tumor = np.where(meta["condition"] == "tumor")[0]
    normal = np.array([np.where((meta["patient"] == meta["patient"].iloc[t])
                                & (meta["condition"] == "normal"))[0][0] for t in tumor])
    zstat, dr = classic_diff_coexpression(X, tumor, normal, iu)
    # the two samples of one patient are NOT independent, so we do not trust
    # the normal approximation; BH on permutation-calibrated cut-off instead.
    p = 2 * stats.norm.sf(np.abs(zstat))
    q = benjamini_hochberg(p)
    sig = q <= fdr
    crit = np.abs(zstat[sig]).min() if sig.any() else np.inf
    null = []
    for _ in range(n_perm):
        swap = rng.random(len(tumor)) < 0.5
        t_idx = np.where(swap, normal, tumor)
        n_idx = np.where(swap, tumor, normal)
        zs, _ = classic_diff_coexpression(X, t_idx, n_idx, iu)
        null.append(int((np.abs(zs) >= crit).sum()))
    return dict(stat=zstat, dr=dr, sig=sig, q=q, crit=crit,
                null_counts=np.array(null), n_sig=int(sig.sum()))


def top_k_set(stat, k_frac):
    """Indices of the top k_frac of edges ranked by |statistic|."""
    k = max(1, int(round(k_frac * stat.size)))
    return set(np.argpartition(-np.abs(stat), k - 1)[:k].tolist())


def jaccard(a, b):
    return len(a & b) / len(a | b) if (a | b) else np.nan


def main():
    expr, meta = load_data()
    meta = meta.reset_index(drop=True)

    # ======================================================================
    # 1. Parameter sweep: number of genes
    # ======================================================================
    log("=== 1. Sweep: number of genes (full cohort, 36 patients) ===")
    rows = []
    for n_genes in [250, 500, 1000, 1500]:
        sub = select_variable_genes(expr, n_genes)
        _, nets = lioness_pearson(sub.values, verbose=False)
        res = differential_edges(nets, meta, n_perm=100)
        rows.append(dict(n_genes=n_genes, n_edges=nets.shape[1],
                         n_significant=res["n_sig"],
                         pct_significant=100 * res["n_sig"] / nets.shape[1],
                         null_mean=res["null_counts"].mean(),
                         null_p95=np.percentile(res["null_counts"], 95),
                         emp_fdr_mean=res["emp_fdr"]))
        del nets
    sweep_genes = pd.DataFrame(rows)
    sweep_genes.to_csv(RESULTS_DIR / "sweep_n_genes.csv", index=False)
    log(sweep_genes.round(3).to_string(index=False))

    # ======================================================================
    # 2. Parameter sweep: cohort size and sparsity cut-off + 3. comparison
    # ======================================================================
    N_GENES = 500
    sub = select_variable_genes(expr, N_GENES)
    X = sub.values
    iu = upper_indices(N_GENES)
    patients = np.array(sorted(meta["patient"].unique()))
    tumor_all = np.array([np.where((meta.patient == p) & (meta.condition == "tumor"))[0][0]
                          for p in patients])
    normal_all = np.array([np.where((meta.patient == p) & (meta.condition == "normal"))[0][0]
                           for p in patients])

    def run_methods(patient_subset):
        """Return |statistic| vectors for LIONESS-paired-t and classic method
        using only the given patients."""
        cols = np.r_[tumor_all[patient_subset], normal_all[patient_subset]]
        Xs = X[:, cols]
        m = len(patient_subset)
        meta_s = pd.DataFrame({"patient": np.r_[patients[patient_subset], patients[patient_subset]],
                               "condition": ["tumor"] * m + ["normal"] * m})
        _, nets = lioness_pearson(Xs, verbose=False)
        res = differential_edges(nets, meta_s, n_perm=0)
        z, dr = classic_diff_coexpression(Xs, np.arange(m), np.arange(m, 2 * m), iu)
        return res["t"], res["n_sig"], z, nets, meta_s

    log("\n=== 2/3. Cohort-size sweep and method comparison "
        f"({N_GENES} genes; reference = all 36 patients) ===")
    all_idx = np.arange(len(patients))
    t_ref, nsig_ref, z_ref, nets_ref, meta_ref = run_methods(all_idx)
    ks = [0.001, 0.01, 0.05]
    ref_sets = {("LIONESS", k): top_k_set(t_ref, k) for k in ks}
    ref_sets.update({("Classic", k): top_k_set(z_ref, k) for k in ks})

    # agreement of the two methods on the full cohort
    rho_methods, _ = stats.spearmanr(np.abs(t_ref), np.abs(z_ref))
    log(f"Full cohort: LIONESS-paired-t significant edges = {nsig_ref}")
    for k in ks:
        log(f"  top {100 * k:g}% edges: Jaccard(LIONESS, classic) = "
            f"{jaccard(ref_sets[('LIONESS', k)], ref_sets[('Classic', k)]):.3f}")
    log(f"  Spearman correlation of |statistics| between methods: {rho_methods:.3f}")

    stab = []
    n_rep = 10
    for m in [9, 18, 27]:
        for rep in range(n_rep):
            subset = np.sort(rng.choice(len(patients), size=m, replace=False))
            t_s, nsig, z_s, _, _ = run_methods(subset)
            for k in ks:
                stab.append(dict(m_patients=m, rep=rep, k_pct=100 * k,
                                 method="LIONESS (paired t)",
                                 jaccard=jaccard(top_k_set(t_s, k), ref_sets[("LIONESS", k)]),
                                 n_sig_fdr5=nsig))
                stab.append(dict(m_patients=m, rep=rep, k_pct=100 * k,
                                 method="Classic diff. co-expression",
                                 jaccard=jaccard(top_k_set(z_s, k), ref_sets[("Classic", k)]),
                                 n_sig_fdr5=np.nan))
    stab = pd.DataFrame(stab)
    stab.to_csv(RESULTS_DIR / "stability_subsampling.csv", index=False)
    table = (stab.groupby(["method", "k_pct", "m_patients"])["jaccard"]
             .agg(["mean", "std"]).round(3))
    log("\nStability: Jaccard overlap with full-cohort top edges "
        f"(mean/std over {n_rep} random patient subsets)")
    log(table.to_string())
    nsig_table = (stab[stab.method.str.startswith("LIONESS")]
                  .groupby("m_patients")["n_sig_fdr5"].agg(["mean", "std"]).round(0))
    log("\nLIONESS: edges significant at FDR 5% (of "
        f"{len(iu[0])}) when using only m patients:")
    log(nsig_table.to_string())

    # --- classic method: significance and null on full data (1000 genes) ---
    sub1000 = select_variable_genes(expr, 1000)
    iu1000 = upper_indices(1000)
    cres = classic_with_null(sub1000.values, meta, iu1000, n_perm=100)
    log("\nClassic differential co-expression, 1000 genes, all 36 patients:")
    log(f"  significant edges (BH FDR 5%, normal approximation): {cres['n_sig']}")
    log(f"  label-swap null: mean {cres['null_counts'].mean():.0f} edges pass the same cut-off "
        f"(empirical FDR {cres['null_counts'].mean() / max(cres['n_sig'], 1):.3f})")

    # ======================================================================
    # Figures
    # ======================================================================
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    axes[0].bar(sweep_genes["n_genes"].astype(str), sweep_genes["pct_significant"],
                color="steelblue")
    axes[0].set_xlabel("number of genes")
    axes[0].set_ylabel("% of edges significant (FDR 5%)")
    axes[0].set_title("Sweep 1: number of genes")
    sub_k = stab[stab.k_pct == 1.0]
    for method, col in [("LIONESS (paired t)", "tab:red"), ("Classic diff. co-expression", "tab:green")]:
        g = sub_k[sub_k.method == method].groupby("m_patients")["jaccard"]
        axes[1].errorbar(g.mean().index, g.mean().values, yerr=g.std().values,
                         marker="o", capsize=3, label=method, color=col)
    axes[1].set_xticks([9, 18, 27])
    axes[1].set_xlabel("number of patients used")
    axes[1].set_ylabel("Jaccard overlap with full-cohort top 1% edges")
    axes[1].set_title("Sweep 2: stability vs cohort size")
    axes[1].legend()
    fig.tight_layout()
    fig.savefig(RESULTS_DIR / "fig4_sweeps.png", dpi=150)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(5.5, 5))
    pick = rng.choice(len(t_ref), 20000, replace=False)
    ax.scatter(np.abs(z_ref[pick]), np.abs(t_ref[pick]), s=3, alpha=0.3)
    ax.set_xlabel("|z| classic differential co-expression")
    ax.set_ylabel("|t| LIONESS paired test")
    ax.set_title(f"Method comparison ({N_GENES} genes), rho = {rho_methods:.2f}")
    fig.tight_layout()
    fig.savefig(RESULTS_DIR / "fig5_method_comparison.png", dpi=150)
    plt.close(fig)

    (RESULTS_DIR / "summary_sweep_comparison.txt").write_text("\n".join(lines))
    print("\nSaved figures and tables to results/")


if __name__ == "__main__":
    main()
