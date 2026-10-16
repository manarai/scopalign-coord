#!/usr/bin/env python3
"""Content fingerprint of staged data files (platform-independent).

h5ad: sha256 over shape, CSR X (indptr int64, indices int32, data float64,
sorted indices), obs_names, var_names, obs['qc_pass'] (uint8, if present),
obs['cell_type'] (if present). Other files: sha256 of the bytes.
Usage: python checker/fingerprint.py <file> [<file> ...]
"""
import hashlib, sys
from pathlib import Path

def fp_h5ad(path):
    import anndata as ad, numpy as np
    from scipy import sparse
    a = ad.read_h5ad(path)
    h = hashlib.sha256()
    h.update(np.asarray(a.shape, dtype=np.int64).tobytes())
    X = sparse.csr_matrix(a.X); X.sort_indices()
    h.update(X.indptr.astype(np.int64).tobytes()); h.update(X.indices.astype(np.int32).tobytes())
    h.update(X.data.astype(np.float64).tobytes())
    h.update("\n".join(map(str, a.obs_names)).encode()); h.update("\n".join(map(str, a.var_names)).encode())
    if "qc_pass" in a.obs:
        h.update(a.obs["qc_pass"].astype(bool).to_numpy().astype(np.uint8).tobytes())
    if "cell_type" in a.obs:
        h.update("\n".join(a.obs["cell_type"].astype(str)).encode())
    return h.hexdigest()

def fp(path):
    p = Path(path)
    return fp_h5ad(p) if p.suffix == ".h5ad" else hashlib.sha256(p.read_bytes()).hexdigest()

if __name__ == "__main__":
    for f in sys.argv[1:]:
        print(f"{fp(f)}  {f}")
