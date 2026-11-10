"""
LIONESS Project - Step 1: Download and prepare the dataset
===========================================================

Dataset : GEO series GSE15471 (Badea et al., 2008)
          "Whole-tissue gene expression study of pancreatic ductal
          adenocarcinoma" - Affymetrix Human Genome U133 Plus 2.0 arrays (GPL570)
          36 patients, each with one tumour sample and one matched normal
          pancreas sample (plus a few technical replicate arrays).
Source  : NCBI GEO (public, no login needed)

What this script does
  1. Downloads the series matrix (RMA-normalised, log2 expression per probe)
     and the GPL570 probe annotation from the NCBI GEO FTP server.
  2. Reads sample information (patient ID, tumour / normal).
  3. Removes technical replicate arrays ("_rep") so every patient has exactly
     one tumour and one normal sample.
  4. Converts probes -> gene symbols (see comments below).
  5. Saves a gene x sample table and a metadata table in data/.

The raw downloads are stored in data/raw/ and are NOT committed to GitHub
(see .gitignore). Re-running this script recreates everything.
"""

import gzip
import io
import urllib.request
from pathlib import Path

import pandas as pd

from lioness_utils import DATA_DIR, EXPR_FILE, META_FILE

RAW_DIR = DATA_DIR / "raw"

SERIES_URL = ("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE15nnn/GSE15471/"
              "matrix/GSE15471_series_matrix.txt.gz")
ANNOT_URL = ("https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPLnnn/GPL570/"
             "annot/GPL570.annot.gz")


def download(url, dest):
    """Download url -> dest (skipped if the file already exists)."""
    if dest.exists() and dest.stat().st_size > 0:
        print(f"  already downloaded: {dest.name}")
        return
    print(f"  downloading {url}")
    # NCBI is happy with a plain python user agent
    with urllib.request.urlopen(url, timeout=120) as resp, open(dest, "wb") as out:
        out.write(resp.read())
    print(f"  saved {dest.name} ({dest.stat().st_size / 1e6:.1f} MB)")


def read_series_matrix(path):
    """Parse a GEO series matrix file into (expression table, sample table)."""
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as fh:
        lines = fh.read().splitlines()

    # --- sample information lines start with "!Sample_" -------------------
    def values(line):
        return [v.strip('"') for v in line.split("\t")[1:]]

    info = {"geo_accession": None, "title": None, "patient": None, "condition": None}
    for line in lines:
        if line.startswith("!Sample_geo_accession"):
            info["geo_accession"] = values(line)
        elif line.startswith("!Sample_title"):
            info["title"] = values(line)
        elif line.startswith("!Sample_characteristics_ch1"):
            vals = values(line)
            if vals[0].startswith("patient:"):
                info["patient"] = [v.split(":")[1].strip() for v in vals]
            elif vals[0].startswith("sample:"):
                info["condition"] = [v.split(":")[1].strip() for v in vals]
    meta = pd.DataFrame(info).set_index("geo_accession")

    # --- expression table sits between two marker lines --------------------
    start = lines.index("!series_matrix_table_begin") + 1
    end = lines.index("!series_matrix_table_end")
    expr = pd.read_csv(io.StringIO("\n".join(lines[start:end])), sep="\t", index_col=0)
    expr.columns = [c.strip('"') for c in expr.columns]
    expr.index = [str(i).strip('"') for i in expr.index]
    return expr, meta


def read_platform_annotation(path):
    """Return a Series: probe ID -> gene symbol from the GPL570 annotation."""
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as fh:
        lines = fh.read().splitlines()
    start = lines.index("!platform_table_begin") + 1
    end = lines.index("!platform_table_end")
    table = pd.read_csv(io.StringIO("\n".join(lines[start:end])), sep="\t",
                        usecols=["ID", "Gene symbol"], index_col="ID",
                        dtype=str)
    return table["Gene symbol"]


def probes_to_genes(expr, probe2gene):
    """Collapse probe-level data to gene level.

    Choices (kept simple on purpose):
      * probes without a symbol, or mapping to several genes ("///"), are
        dropped because their signal cannot be assigned to one gene;
      * Affymetrix control probes (AFFX-) are dropped;
      * when several probes measure the same gene we keep the probe with the
        highest mean expression across samples (a common, simple rule).
    """
    symbols = probe2gene.reindex(expr.index)
    keep = symbols.notna() & ~symbols.fillna("").str.contains("///") \
        & ~expr.index.to_series().str.startswith("AFFX")
    expr = expr[keep.values]
    symbols = symbols[keep.values]

    # best probe per gene = probe with the highest average expression
    order = expr.mean(axis=1).sort_values(ascending=False).index
    expr = expr.loc[order]
    symbols = symbols.loc[order]
    expr = expr[~symbols.duplicated().values]
    expr.index = symbols[~symbols.duplicated()].values
    return expr


def main():
    RAW_DIR.mkdir(parents=True, exist_ok=True)
    print("Step 1/3: downloading from NCBI GEO ...")
    series_path = RAW_DIR / "GSE15471_series_matrix.txt.gz"
    annot_path = RAW_DIR / "GPL570.annot.gz"
    download(SERIES_URL, series_path)
    download(ANNOT_URL, annot_path)

    print("Step 2/3: reading files ...")
    expr_probe, meta = read_series_matrix(series_path)
    probe2gene = read_platform_annotation(annot_path)
    print(f"  probe-level matrix: {expr_probe.shape[0]} probes x "
          f"{expr_probe.shape[1]} arrays")

    # Drop technical replicate arrays (titles end with "_rep"): they are the
    # same biological sample measured twice and would otherwise be counted as
    # independent samples.
    is_rep = meta["title"].str.endswith("_rep")
    print(f"  removing {is_rep.sum()} technical replicate arrays")
    meta = meta[~is_rep]
    expr_probe = expr_probe[meta.index]

    print("Step 3/3: probes -> genes ...")
    expr = probes_to_genes(expr_probe, probe2gene)
    print(f"  gene-level matrix: {expr.shape[0]} genes x {expr.shape[1]} samples")

    # sort samples: all normal first, then tumour, ordered by patient
    meta = meta.sort_values(["condition", "patient"], ascending=[False, True])
    expr = expr[meta.index]

    # sanity checks - every patient must have exactly one normal and one tumour
    counts = meta.groupby(["patient", "condition"]).size().unstack(fill_value=0)
    assert (counts["normal"] == 1).all() and (counts["tumor"] == 1).all(), \
        "unexpected sample layout"
    assert not expr.isna().any().any(), "missing values in expression matrix"

    expr.to_csv(EXPR_FILE)
    meta[["title", "patient", "condition"]].to_csv(META_FILE)
    print(f"\nSaved {EXPR_FILE.relative_to(DATA_DIR.parent)}")
    print(f"Saved {META_FILE.relative_to(DATA_DIR.parent)}")
    print(meta["condition"].value_counts().to_string())
    print("\nReady for 02_run_lioness.py")


if __name__ == "__main__":
    main()
