"""
Small correctness tests for the LIONESS implementation.
Run with:  python -m pytest tests/    (or simply: python tests/test_lioness.py)
"""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import pandas as pd  # noqa: E402

from lioness_utils import (benjamini_hochberg, lioness_pearson,  # noqa: E402
                           lioness_pearson_numpy, upper_indices)


def test_lioness_matches_definition():
    """Each sample network must equal N*w_all - (N-1)*w_without_q, computed
    here with a slow, obviously-correct loop over gene pairs."""
    rng = np.random.default_rng(0)
    X = rng.normal(size=(6, 12))  # 6 genes, 12 samples
    agg, nets = lioness_pearson(X, verbose=False)
    iu = upper_indices(6)
    N = X.shape[1]
    for q in range(N):
        rest = np.delete(X, q, axis=1)
        for e, (i, j) in enumerate(zip(*iu)):
            w_all = np.corrcoef(X[i], X[j])[0, 1]
            w_rest = np.corrcoef(rest[i], rest[j])[0, 1]
            assert np.isclose(nets[q, e], N * w_all - (N - 1) * w_rest, atol=1e-5)


def test_edges_are_unique_pairs():
    X = np.random.default_rng(1).normal(size=(10, 8))
    agg, nets = lioness_pearson(X, verbose=False)
    assert agg.size == 10 * 9 // 2
    assert nets.shape == (8, 45)


def test_netzoopy_matches_numpy_reference():
    """netZooPy's LIONESS must agree with an independent numpy implementation,
    and must keep the gene order of the input (no alphabetical re-sorting)."""
    rng = np.random.default_rng(2)
    X = rng.normal(size=(15, 10))
    genes = ["ZZZ", "AAA", "MMM"] + [f"G{i:02d}" for i in range(12)]
    df = pd.DataFrame(X, index=genes, columns=[f"S{j}" for j in range(10)])
    agg_nz, nets_nz = lioness_pearson(df, verbose=False)
    agg_np, nets_np = lioness_pearson_numpy(X, verbose=False)
    assert np.allclose(agg_nz, agg_np, atol=1e-10)
    assert np.allclose(nets_nz, nets_np, atol=1e-4)


def test_benjamini_hochberg():
    p = np.array([0.01, 0.04, 0.03, 0.005])
    q = benjamini_hochberg(p)
    assert np.allclose(q, [0.02, 0.04, 0.04, 0.02])


if __name__ == "__main__":
    test_lioness_matches_definition()
    test_edges_are_unique_pairs()
    test_netzoopy_matches_numpy_reference()
    test_benjamini_hochberg()
    print("all tests passed")
