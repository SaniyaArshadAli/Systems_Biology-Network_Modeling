"""
LIONESS Project - Step 2: Build sample-specific networks with LIONESS
=====================================================================

Base network : Pearson correlation between genes (computed over samples).
LIONESS      : w_q = N * w_all - (N - 1) * w_(all without q)

Usage
    python 02_run_lioness.py                 # default: 1000 most variable genes
    python 02_run_lioness.py --n-genes 500

Output (in results/, not committed because of its size)
    lioness_networks_<n>genes.npy   (samples x edges, float32)
    aggregate_network_<n>genes.npy  (edges,)
    genes_<n>genes.txt              (gene order used for the edges)

Edges are stored as the upper triangle of the gene x gene matrix (pairs i < j).
The original version of this project stored the full symmetric matrix including
the diagonal, which double-counted every edge; storing unique pairs fixes that.
"""

import argparse

import numpy as np

from lioness_utils import RESULTS_DIR, lioness_pearson, load_data, select_variable_genes


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--n-genes", type=int, default=1000,
                        help="number of most-variable genes to use")
    args = parser.parse_args()

    expr, meta = load_data()
    print(f"Loaded {expr.shape[0]} genes x {expr.shape[1]} samples "
          f"({(meta.condition == 'tumor').sum()} tumour, "
          f"{(meta.condition == 'normal').sum()} normal)")

    sub = select_variable_genes(expr, args.n_genes)
    print(f"Using the {sub.shape[0]} most variable genes")

    # LIONESS uses ALL samples jointly; the group labels are never used here.
    agg, networks = lioness_pearson(sub.values)
    print(f"Edges per network: {agg.size}  (unique gene pairs)")

    # quick sanity checks
    assert np.isfinite(networks).all(), "non-finite values in LIONESS networks"
    print(f"Weight range of sample networks: {networks.min():.2f} to {networks.max():.2f}")
    print(f"Aggregate Pearson range:         {agg.min():.2f} to {agg.max():.2f}")

    RESULTS_DIR.mkdir(exist_ok=True)
    tag = f"{args.n_genes}genes"
    np.save(RESULTS_DIR / f"lioness_networks_{tag}.npy", networks)
    np.save(RESULTS_DIR / f"aggregate_network_{tag}.npy", agg.astype(np.float32))
    (RESULTS_DIR / f"genes_{tag}.txt").write_text("\n".join(sub.index))
    print(f"Saved networks for {tag} to results/")


if __name__ == "__main__":
    main()
