"""
Shared helper functions for the LIONESS project.

Everything that is used by more than one script lives here so that the
numbered scripts (01 ... 04) stay short and readable.

LIONESS in one line (Kuijjer et al. 2019, iScience):
    w_q = N * w_all - (N - 1) * w_(all without q)
where w is ANY network-inference function (here: Pearson correlation),
N is the number of samples and q is the sample we want a network for.

The sample-specific networks are computed with netZooPy (the Python
implementation of the Network Zoo, https://netzoo.github.io/), using
PANDA in "correlation only" mode (no motif / PPI prior) as the base network
and netZooPy's Lioness class on top of it. A plain-numpy version of the same
equation (`lioness_pearson_numpy`) is kept only as an independent reference
for the unit tests.
"""

import contextlib
import io
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Paths (relative to this file, so the project works on any machine)
# ---------------------------------------------------------------------------
ROOT = Path(__file__).resolve().parent
DATA_DIR = ROOT / "data"
RESULTS_DIR = ROOT / "results"
EXPR_FILE = DATA_DIR / "expression_gene_level.csv.gz"
META_FILE = DATA_DIR / "sample_metadata.csv"


def load_data():
    """Load the gene-level expression matrix (genes x samples) and metadata."""
    if not EXPR_FILE.exists():
        raise FileNotFoundError(
            f"{EXPR_FILE} not found. Run `python 01_download_data.py` first."
        )
    expr = pd.read_csv(EXPR_FILE, index_col=0)
    meta = pd.read_csv(META_FILE, index_col=0)
    # make sure the sample order is identical in both tables
    expr = expr[meta.index]
    return expr, meta


def select_variable_genes(expr, n_genes):
    """Keep the n_genes genes with the highest variance across all samples.

    Pearson LIONESS needs an N x N gene-gene matrix per sample, so using all
    ~20,000 genes is not practical. Note: variance is computed over tumour AND
    normal samples, so this filter favours genes that differ between the two
    groups. This is a (documented) bias of the project.
    """
    variances = expr.var(axis=1)
    top = variances.sort_values(ascending=False).index[:n_genes]
    return expr.loc[top]


def upper_indices(n_genes):
    """Row/column indices of the upper triangle (gene pairs i < j)."""
    return np.triu_indices(n_genes, k=1)


def pearson_edges(X, iu):
    """Pearson correlation between all gene pairs.

    X  : genes x samples array
    iu : output of upper_indices -> only unique pairs (no diagonal, no
         duplicates) are returned as a 1-D vector.
    """
    return np.corrcoef(X)[iu]


def lioness_pearson(X, verbose=True, ncores=1):
    """Run LIONESS with Pearson correlation as the base network (netZooPy).

    netZooPy's PANDA class is run without a motif or PPI prior, in which case
    its network is the Pearson co-expression matrix of the expression data.
    netZooPy's Lioness class then applies
        w_q = N * w_all - (N - 1) * w_(all without q)
    for every sample q. Each sample network is written by netZooPy to a
    temporary folder (save_single=True) and read back one at a time, keeping
    only the unique gene pairs (upper triangle) as float32. This avoids
    netZooPy's default edge x sample table with gene-name labels, which needs
    several GB of RAM for 1000+ genes.

    Parameters
    ----------
    X : array or DataFrame, genes x samples
    verbose : print progress
    ncores : number of CPU cores passed to netZooPy's Lioness

    Returns
    -------
    agg      : (E,) aggregate network (all samples), E = genes*(genes-1)/2
    networks : (samples, E) float32 array, one row = one sample-specific network
    """
    from netZooPy.lioness import Lioness
    from netZooPy.panda import Panda

    if isinstance(X, pd.DataFrame):
        df = X
    else:
        X = np.asarray(X, dtype=np.float64)
        df = pd.DataFrame(X, index=[f"gene{i}" for i in range(X.shape[0])],
                          columns=[f"sample{j}" for j in range(X.shape[1])])
    df = df.astype(np.float64)
    n_genes, n_samples = df.shape
    iu = upper_indices(n_genes)

    # netZooPy prints a lot of timing messages; keep the console readable.
    quiet = contextlib.redirect_stdout(io.StringIO())
    with tempfile.TemporaryDirectory(prefix="lioness_netzoopy_") as tmp, quiet:
        # modeProcess="legacy" keeps the genes in the given order (the
        # "union"/"intersection" modes sort them alphabetically).
        panda = Panda(df, None, None, modeProcess="legacy", remove_missing=False,
                      keep_expression_matrix=True, save_memory=False)
        agg = np.asarray(panda.panda_network.values)[iu]
        Lioness(panda, computing="cpu", precision="double", ncores=ncores,
                save_dir=tmp, save_fmt="npy", save_single=True,
                ignore_final=True)
        # netZooPy names the files lioness.<sample>.<index>.npy
        networks = np.empty((n_samples, agg.size), dtype=np.float32)
        for q in range(n_samples):
            path = Path(tmp) / f"lioness.{panda.expression_samples[q]}.{q}.npy"
            networks[q] = np.load(path)[iu]
    if verbose:
        print(f"  netZooPy LIONESS: {n_samples}/{n_samples} sample networks done")
    return agg, networks


def lioness_pearson_numpy(X, verbose=True):
    """Reference implementation of Pearson LIONESS in plain numpy.

    Not used by the analysis scripts (they use netZooPy via
    `lioness_pearson`); kept so the tests can check netZooPy against an
    independent implementation of the LIONESS equation.

    Parameters
    ----------
    X : array, genes x samples

    Returns
    -------
    agg      : (E,) aggregate network (all samples), E = genes*(genes-1)/2
    networks : (samples, E) float32 array, one row = one sample-specific network
    """
    n_genes, n_samples = X.shape
    iu = upper_indices(n_genes)
    agg = pearson_edges(X, iu)

    networks = np.empty((n_samples, agg.size), dtype=np.float32)
    for q in range(n_samples):
        # leave sample q out and recompute the network
        loo = pearson_edges(np.delete(X, q, axis=1), iu)
        # LIONESS equation
        networks[q] = n_samples * agg - (n_samples - 1) * loo
        if verbose and (q + 1) % 10 == 0:
            print(f"  LIONESS: {q + 1}/{n_samples} samples done")
    return agg, networks


def benjamini_hochberg(p):
    """Benjamini-Hochberg FDR-adjusted p-values (q-values)."""
    p = np.asarray(p, dtype=float)
    n = p.size
    order = np.argsort(p)
    ranked = p[order] * n / (np.arange(n) + 1)
    # enforce monotonicity from the largest p-value downwards
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    q = np.empty(n)
    q[order] = np.minimum(ranked, 1.0)
    return q


def paired_differences(W, meta):
    """Return tumour-minus-normal differences for each patient.

    W    : (samples, features) array (edge weights or expression values)
    meta : metadata with columns 'patient' and 'condition'
    Output shape: (patients, features)
    """
    tumor, normal = [], []
    for patient, grp in meta.groupby("patient"):
        t = np.where((meta["patient"] == patient) & (meta["condition"] == "tumor"))[0]
        n = np.where((meta["patient"] == patient) & (meta["condition"] == "normal"))[0]
        if len(t) == 1 and len(n) == 1:
            tumor.append(t[0])
            normal.append(n[0])
    return W[tumor] - W[normal], np.array(tumor), np.array(normal)


def differential_edges(W, meta, n_perm=200, fdr=0.05, seed=0):
    """Find edges whose weight differs between tumour and normal (paired).

    Because every patient contributes one tumour and one matched normal sample,
    we test the per-patient difference (tumour - normal) of each edge with a
    paired t-test and control the false discovery rate with Benjamini-Hochberg.

    As an extra check that does not rely on the t-distribution (LIONESS
    networks are computed from overlapping data, so they are not perfectly
    independent), we build a null distribution by randomly swapping the
    tumour/normal labels inside each patient pair and counting how many edges
    would still pass the same |t| cut-off.

    Returns a dict with t-statistics, q-values, mean differences, null counts.
    """
    from scipy import stats

    D, _, _ = paired_differences(W.astype(np.float64), meta)
    n = D.shape[0]
    mean_diff = D.mean(axis=0)
    t = np.nan_to_num(paired_t_stats(D))
    p = 2 * stats.t.sf(np.abs(t), df=n - 1)
    q = benjamini_hochberg(p)
    sig = q <= fdr
    crit = np.abs(t[sig]).min() if sig.any() else np.inf

    # label-swap null: flipping the sign of a patient's difference is the same
    # as swapping that patient's tumour/normal labels. The sum of squares of D
    # does not change under sign flips, so the new t is cheap to compute.
    rng = np.random.default_rng(seed)
    sumsq = (D ** 2).sum(axis=0)
    null_counts = []
    for _ in range(n_perm):
        s = rng.choice([-1.0, 1.0], size=n)
        m = (s @ D) / n
        var = (sumsq - n * m ** 2) / (n - 1)
        tp = m / np.sqrt(np.maximum(var, 1e-12) / n)
        null_counts.append(int((np.abs(tp) >= crit).sum()))
    null_counts = np.array(null_counts)
    n_sig = int(sig.sum())
    emp_fdr = null_counts.mean() / n_sig if (n_sig and n_perm) else np.nan
    return dict(t=t, p=p, q=q, sig=sig, mean_diff=mean_diff, crit=crit,
                n_sig=n_sig, null_counts=null_counts, emp_fdr=emp_fdr)


def paired_t_stats(D):
    """Paired t-statistic for every column of the difference matrix D."""
    n = D.shape[0]
    mean = D.mean(axis=0)
    sd = D.std(axis=0, ddof=1)
    sd = np.where(sd == 0, np.nan, sd)
    return mean / (sd / np.sqrt(n))
