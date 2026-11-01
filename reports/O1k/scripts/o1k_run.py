#!/usr/bin/env python
"""Reproducible runner for O1k: synthetic positive control for O1g E5b.

Copy this file to reports/O1k/scripts/o1k_run.py in the O1k worktree and run
it with that worktree's virtual environment.  It deliberately treats O1g as
read-only: the only O1g access is the check-stage import and loading of its
cached Test R fits.
"""

import argparse
import hashlib
import importlib.util
import json
import os
import subprocess
import sys
import tempfile
from concurrent.futures import ProcessPoolExecutor
from dataclasses import asdict
from datetime import datetime, timezone
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import torch
from scipy.stats import ortho_group
from sklearn.metrics import silhouette_score
from sklearn.neighbors import NearestNeighbors

from scjdo2.energy_drift import EnergyDriftAdapter, EnergyDriftConfig, _EnergyMLP


# This layout is valid after this script is copied to reports/O1k/scripts.
ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "reports" / "O1k"
CACHE = OUT / "cache"
ARTIFACTS = OUT / "artifacts"

# O1g is intentionally not inferred from the current worktree.  It is an
# issued, read-only input to this diagnostic, as specified by O1k.
O1G = Path.home() / "Documents" / "scopalign-manus-O1g" / "reports" / "O1g"
O1G_SCRIPT = O1G / "scripts" / "o1g_run.py"

D = 10
K = 4
GAMMAS = (1.0, 1.15, 1.33, 1.6)
NS = (1000, 3000, 10000)
REPS = (0, 1, 2)
FIRST5 = (20261014, 20261016, 20261017, 20261018, 20261019)
KNOWN_FAILURES = [
    "tests/test_phase_b_pitchfork.py::test_pitchfork_saddle_localisation[0.1]",
    "tests/test_phase_b_pitchfork.py::test_pitchfork_saddle_localisation[0.3]",
    "tests/test_phase_b_pitchfork.py::test_pitchfork_saddle_localisation[0.5]",
    "tests/test_phase_b_pitchfork.py::test_pitchfork_half_split_reproducibility",
]
GATED_KEYS = (
    "n_cases",
    "start_overlap",
    "loss_copy_check",
    "symmetry_max",
    "negctrl_err_shape_A",
    "errA_max",
    "errB_g115_n3000",
    "errB_g133_n3000",
    "errB_g160_n3000",
    "errB_g115_n10000",
    "errB_g133_n10000",
    "errB_g160_n10000",
    "frac_hat_below_star_B",
    "fe_u_3000",
    "fe_J_3000",
    "o1k_outcome",
)
CASE_COLUMNS = (
    "rung",
    "gamma",
    "n",
    "r",
    "err",
    "err_shape",
    "err_pos",
    "L_hat",
    "L_star",
    "best_start",
    "n_near_best",
    "fe_u_X",
    "fe_u_Y",
    "fe_J_X",
    "fe_J_Y",
)
WORKER_PROBLEM = None


def chunks(handle):
    """Yield a file in fixed-size blocks for SHA-256."""
    while True:
        block = handle.read(1 << 20)
        if not block:
            return
        yield block


def sha(path):
    """Return the SHA-256 of a file, refusing silently missing inputs."""
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(path)
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in chunks(handle):
            digest.update(block)
    return digest.hexdigest()


def atomic_path(path, suffix=".tmp"):
    """Create a same-directory temporary path suitable for os.replace."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    handle = tempfile.NamedTemporaryFile(
        dir=path.parent, prefix=f".{path.name}.", suffix=suffix, delete=False
    )
    handle.close()
    return Path(handle.name)


def atomic_json(path, value):
    tmp = atomic_path(path)
    try:
        tmp.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
        os.replace(tmp, path)
    finally:
        tmp.unlink(missing_ok=True)


def atomic_dataframe(path, frame):
    tmp = atomic_path(path)
    try:
        frame.to_csv(tmp, index=False)
        os.replace(tmp, path)
    finally:
        tmp.unlink(missing_ok=True)


def atomic_npz(path, **arrays):
    tmp = atomic_path(path, suffix=".npz")
    try:
        np.savez(tmp, **arrays)
        os.replace(tmp, path)
    finally:
        tmp.unlink(missing_ok=True)


def atomic_torch_save(path, value):
    tmp = atomic_path(path, suffix=".pt")
    try:
        torch.save(value, tmp)
        os.replace(tmp, path)
    finally:
        tmp.unlink(missing_ok=True)


def make_meta():
    """Create output directories and preserve the first stage start time."""
    OUT.mkdir(parents=True, exist_ok=True)
    CACHE.mkdir(parents=True, exist_ok=True)
    ARTIFACTS.mkdir(parents=True, exist_ok=True)
    gitignore = OUT / ".gitignore"
    if not gitignore.exists():
        gitignore.write_text("cache/\n__pycache__/\nrun.log\nrun.pid\n")
    meta_path = OUT / "run_meta.json"
    if not meta_path.exists():
        atomic_json(meta_path, {"started": datetime.now(timezone.utc).isoformat()})
    return json.loads(meta_path.read_text())


def gamma_tag(gamma):
    return f"g{gamma:.2f}".replace(".", "p")


def case_id(gamma, n, r):
    return f"{gamma_tag(gamma)}_n{int(n)}_r{int(r)}"


def case_specs():
    """Return the 36 issued cases in a deterministic, independent order."""
    for gamma in GAMMAS:
        for n in NS:
            for r in REPS:
                yield float(gamma), int(n), int(r), case_id(gamma, n, r)


def landscape_arrays():
    """Construct all fixed landscape quantities using the exact draw order."""
    rng = np.random.default_rng(20264200)
    means = np.zeros((K, D), dtype=np.float64)
    # O1k explicitly specifies one standard_normal(3) call for each well.
    for k in range(K):
        means[k, :3] = 3.0 * rng.standard_normal(3)

    rotations = np.stack(
        [np.asarray(ortho_group.rvs(D, random_state=20264210 + k), dtype=np.float64) for k in range(K)],
        axis=0,
    )
    sigma2 = np.stack(
        [np.asarray(gamma ** (5.5 - np.arange(1, D + 1)), dtype=np.float64) for gamma in GAMMAS],
        axis=0,
    )
    q_stars = np.stack(
        [np.asarray(ortho_group.rvs(D, random_state=20264100 + r), dtype=np.float64) for r in REPS],
        axis=0,
    )
    return {
        "means": means,
        "rotations": rotations,
        "gammas": np.asarray(GAMMAS, dtype=np.float64),
        "sigma2": sigma2,
        "q_stars": q_stars,
    }


def ensure_landscape():
    """Write the required landscape artifact once; values remain deterministic."""
    arrays = landscape_arrays()
    path = ARTIFACTS / "landscape.npz"
    if not path.exists():
        atomic_npz(path, **arrays)
    return arrays


def sigma2_for_gamma(gamma):
    return np.asarray(gamma ** (5.5 - np.arange(1, D + 1)), dtype=np.float64)


def sample_case(gamma, n, r, landscape=None):
    """Generate a case exactly in O1k's stated label/normal draw order."""
    if landscape is None:
        landscape = landscape_arrays()
    means = landscape["means"]
    rotations = landscape["rotations"]
    sigma = np.sqrt(sigma2_for_gamma(gamma))
    q_star = np.asarray(landscape["q_stars"][int(r)], dtype=np.float64)

    rng = np.random.default_rng(20264300 + int(r))
    x_labels = rng.integers(0, K, int(n))
    x_eps = rng.standard_normal((int(n), D))
    y_labels = rng.integers(0, K, int(n))
    y_eps = rng.standard_normal((int(n), D))

    def draw(labels, eps):
        scaled = eps * sigma[None, :]
        return means[labels] + np.einsum("nij,nj->ni", rotations[labels], scaled)

    x = draw(x_labels, x_eps)
    y_prime = draw(y_labels, y_eps)
    y = y_prime @ q_star.T
    return x, y, x_labels, y_labels, q_star


class ExactEnergy(torch.nn.Module):
    """The exact float64 Gaussian-mixture energy in a specified row frame."""

    def __init__(self, means, rotations, sigma2, coordinate_to_base=None):
        super().__init__()
        means = np.asarray(means, dtype=np.float64)
        rotations = np.asarray(rotations, dtype=np.float64)
        sigma2 = np.asarray(sigma2, dtype=np.float64)
        if coordinate_to_base is None:
            coordinate_to_base = np.eye(means.shape[1], dtype=np.float64)
        coordinate_to_base = np.asarray(coordinate_to_base, dtype=np.float64)
        precision = np.einsum(
            "kij,j,klj->kil", rotations, 1.0 / sigma2, rotations
        )
        # The einsum above is R diag(1/sigma2) R^T.  Symmetrising only removes
        # round-off asymmetry and does not alter the defined covariance.
        precision = 0.5 * (precision + np.swapaxes(precision, 1, 2))
        log_norm = -0.5 * (
            means.shape[1] * np.log(2.0 * np.pi) + np.sum(np.log(sigma2))
        ) - np.log(means.shape[0])
        self.register_buffer("means", torch.as_tensor(means, dtype=torch.float64))
        self.register_buffer("precision", torch.as_tensor(precision, dtype=torch.float64))
        self.register_buffer("coordinate_to_base", torch.as_tensor(coordinate_to_base, dtype=torch.float64))
        self.register_buffer("log_norm", torch.as_tensor(log_norm, dtype=torch.float64))

    def forward(self, x):
        z = x @ self.coordinate_to_base
        delta = z[:, None, :] - self.means[None, :, :]
        quadratic = torch.einsum("nki,kij,nkj->nk", delta, self.precision, delta)
        log_density = self.log_norm - 0.5 * quadratic
        return -torch.logsumexp(log_density, dim=1)


def exact_models(gamma, q_star, landscape=None):
    """Return the one-model X/Y exact-energy ensembles specified for rung A."""
    if landscape is None:
        landscape = landscape_arrays()
    base = ExactEnergy(
        landscape["means"], landscape["rotations"], sigma2_for_gamma(gamma)
    ).eval()
    # For a row y, the column convention E(Q*^T y) is y @ Q*.
    rotated = ExactEnergy(
        landscape["means"],
        landscape["rotations"],
        sigma2_for_gamma(gamma),
        coordinate_to_base=q_star,
    ).eval()
    return [base], [rotated]


# These field/ensemble functions are a faithful dynamic-dimension copy of O1g.
def fields(model, x, need_j, track, d=None):
    if d is None:
        d = int(x.shape[1])
    z = x.clone().requires_grad_(True)
    e = model(z)
    g = torch.autograd.grad(
        e.sum(), z, create_graph=(need_j or track), retain_graph=(need_j or track)
    )[0]
    u = -g
    if not need_j:
        return (u, None) if track else (u.detach(), None)
    rows = []
    for j in range(int(d)):
        rows.append(
            -torch.autograd.grad(
                g[:, j].sum(),
                z,
                create_graph=track,
                retain_graph=(track or j < int(d) - 1),
            )[0]
        )
    jacobian = torch.stack(rows, dim=1)
    return (u, jacobian) if track else (u.detach(), jacobian.detach())


def ensemble(models, x, need_j, track, d=None):
    us = []
    js = []
    for model in models:
        u, jacobian = fields(model, x, need_j, track, d=d)
        us.append(u)
        if need_j:
            js.append(jacobian)
    return torch.stack(us).mean(0), torch.stack(js).mean(0) if need_j else None


def fit_workers():
    """The issued runner uses no more than three module-level worker processes."""
    return min(3, max(1, (os.cpu_count() or 1) - 1))


def start_workers():
    return min(3, fit_workers())


def cfg(seed):
    """O1g's issued ESM/Hutchinson scJDO2 configuration, unchanged."""
    return EnergyDriftConfig(
        hidden=96,
        depth=3,
        n_epochs=1500,
        seed=int(seed),
        batch_size=256,
        loss="esm",
        lr=5e-4,
        beta=1.0,
        esm_laplacian="hutchinson",
        device="cpu",
    )


def local_fit(case, side, seed):
    return CACHE / "drifts" / f"{case}_{side}_{int(seed)}.pt"


def load_saved(path):
    return torch.load(path, map_location="cpu", weights_only=False)


def fit_row(case, side, seed, saved):
    return {
        "fit_id": f"{case}_{side}_{int(seed)}",
        "case": case,
        "side": side,
        "seed": int(seed),
        "median_asym": float(saved["symmetry"]["median_asym"]),
        "max_asym": float(saved["symmetry"]["max_asym"]),
    }


def fit_worker(item):
    """Module-level, spawn-safe worker for one cached scJDO2 fit."""
    torch.set_num_threads(1)
    case, side, seed, z, output_path = item
    output_path = Path(output_path)
    if output_path.exists():
        return fit_row(case, side, seed, load_saved(output_path))
    output_path.parent.mkdir(parents=True, exist_ok=True)
    drift = EnergyDriftAdapter(D, config=cfg(seed)).fit(
        np.asarray(z, dtype=np.float32), verbose=True
    )
    symmetry = drift.symmetry_check(np.asarray(z, dtype=np.float32), n_cells=100)
    saved = {
        "state_dict": drift._model.state_dict(),
        "config": asdict(drift.config),
        "symmetry": {key: float(value) for key, value in symmetry.items()},
    }
    atomic_torch_save(output_path, saved)
    return fit_row(case, side, seed, saved)


def save_symmetry(rows):
    """Atomically merge resumable fit rows, keyed by their deterministic id."""
    path = OUT / "fit_symmetry.csv"
    existing = pd.read_csv(path).to_dict("records") if path.exists() else []
    by_id = {row["fit_id"]: row for row in existing}
    by_id.update({row["fit_id"]: row for row in rows})
    ordered = pd.DataFrame(
        [by_id[key] for key in sorted(by_id)],
        columns=("fit_id", "case", "side", "seed", "median_asym", "max_asym"),
    )
    atomic_dataframe(path, ordered)


def fit_stage():
    """Fit (or load) the specified 36 x 2 x 5 rung-B drift models."""
    make_meta()
    landscape = ensure_landscape()
    rows = []
    pending = []
    for gamma, n, r, case in case_specs():
        x, y, _, _, _ = sample_case(gamma, n, r, landscape)
        for side, z in (("X", x), ("Y", y)):
            for seed in FIRST5:
                path = local_fit(case, side, seed)
                if path.exists():
                    rows.append(fit_row(case, side, seed, load_saved(path)))
                else:
                    pending.append((case, side, seed, np.asarray(z, dtype=np.float32), str(path)))
    if pending:
        with ProcessPoolExecutor(max_workers=fit_workers()) as executor:
            rows.extend(list(executor.map(fit_worker, pending)))
    save_symmetry(rows)
    print(
        json.dumps(
            {"fit_rows": len(rows), "new_fits": len(pending), "fit_workers": fit_workers()},
            sort_keys=True,
        ),
        flush=True,
    )


def build_model(saved, d=D):
    model = _EnergyMLP(
        int(d), int(saved["config"]["hidden"]), int(saved["config"]["depth"])
    )
    model.load_state_dict(saved["state_dict"])
    model = model.to(dtype=torch.float64, device="cpu").eval()
    for parameter in model.parameters():
        parameter.requires_grad_(False)
    return model


def get_models(case):
    """Read the five cached rung-B models per copy; never refit in estimate."""
    mx = []
    my = []
    rows = []
    for seed in FIRST5:
        for side, bucket in (("X", mx), ("Y", my)):
            path = local_fit(case, side, seed)
            if not path.exists():
                raise RuntimeError(
                    f"Missing required rung-B fit after --stage fits: {path}"
                )
            saved = load_saved(path)
            bucket.append(build_model(saved, D))
            rows.append(fit_row(case, side, seed, saved))
    save_symmetry(rows)
    return mx, my


def starts(d=D):
    """O1k's twenty fixed starts: identity then random states 1..19."""
    yield 0, np.eye(int(d), dtype=np.float64)
    for i in range(1, 20):
        yield i, np.asarray(ortho_group.rvs(int(d), random_state=20264000 + i), dtype=np.float64)


def start_overlap():
    """Check that none of the twenty issued starts is any issued Q*."""
    best = 0
    trial_starts = [q for _, q in starts(D)]
    for r in REPS:
        q_star = np.asarray(ortho_group.rvs(D, random_state=20264100 + r), dtype=np.float64)
        for q0 in trial_starts:
            best = max(best, int(np.linalg.norm(q0 - q_star, ord="fro") < 1e-8))
    return int(best)


class Problem:
    """O1g E5b loss with the dimension explicitly carried by the instance."""

    def __init__(self, mx, my, x, y, r, d=D):
        self.mx = mx
        self.my = my
        self.d = int(d)
        self.x = torch.as_tensor(np.asarray(x, dtype=np.float64), dtype=torch.float64)
        self.y = torch.as_tensor(np.asarray(y, dtype=np.float64), dtype=torch.float64)
        rng = np.random.default_rng(20261040 + int(r))
        self.xa = self.x[
            rng.choice(len(x), size=min(256, len(x)), replace=False)
        ]
        self.yb = self.y[
            rng.choice(len(y), size=min(256, len(y)), replace=False)
        ]
        self.ux, _ = ensemble(mx, self.x, False, False, d=self.d)
        self.uy, _ = ensemble(my, self.y, False, False, d=self.d)
        _, self.jx = ensemble(mx, self.xa, True, False, d=self.d)
        _, self.jy = ensemble(my, self.yb, True, False, d=self.d)
        self.du = torch.sum(self.ux * self.ux)
        self.dv = torch.sum(self.uy * self.uy)
        self.dx = torch.sum(self.jx * self.jx)
        self.dy = torch.sum(self.jy * self.jy)

    def loss(self, q, track):
        uy, _ = ensemble(self.my, self.x @ q.T, False, track, d=self.d)
        ux, _ = ensemble(self.mx, self.y @ q, False, track, d=self.d)
        dl = (
            torch.sum((uy - self.ux @ q.T) ** 2) / self.du
            + torch.sum((ux - self.uy @ q) ** 2) / self.dv
        ) / 2
        _, jy = ensemble(self.my, self.xa @ q.T, True, track, d=self.d)
        _, jx = ensemble(self.mx, self.yb @ q, True, track, d=self.d)
        qj = torch.einsum("ij,njk,lk->nil", q, self.jx, q)
        qtj = torch.einsum("ji,njk,kl->nil", q, self.jy, q)
        hl = (
            torch.sum((jy - qj) ** 2) / self.dx
            + torch.sum((jx - qtj) ** 2) / self.dy
        ) / 2
        return dl + hl, dl, hl


class Closure:
    """Unchanged O1g LBFGS manifold closure."""

    def __init__(self, problem, q0, skew):
        self.problem = problem
        self.q0 = q0
        self.skew = skew

    def __call__(self):
        self.skew.grad = None
        q = self.q0 @ torch.matrix_exp(self.skew - self.skew.T)
        value, _, _ = self.problem.loss(q, True)
        value.backward()
        return value


def init_start_worker(problem):
    """Module-level initializer required for ProcessPoolExecutor start workers."""
    global WORKER_PROBLEM
    torch.set_num_threads(1)
    WORKER_PROBLEM = problem


def optimise_start(item):
    """Module-level optimiser worker for one one of the twenty fixed starts."""
    torch.set_num_threads(1)
    start, q0_numpy = item
    q0 = torch.as_tensor(q0_numpy, dtype=torch.float64)
    skew = torch.zeros((WORKER_PROBLEM.d, WORKER_PROBLEM.d), dtype=torch.float64, requires_grad=True)
    optimiser = torch.optim.LBFGS(
        [skew],
        lr=1.0,
        max_iter=500,
        tolerance_grad=1e-10,
        tolerance_change=1e-12,
        line_search_fn="strong_wolfe",
    )
    optimiser.step(Closure(WORKER_PROBLEM, q0, skew))
    q = (q0 @ torch.matrix_exp(skew - skew.T)).detach()
    value, dl, hl = WORKER_PROBLEM.loss(q, False)
    return {
        "start": int(start),
        "L_hat": float(value.detach()),
        "D": float(dl.detach()),
        "H": float(hl.detach()),
        "Q": q.cpu().numpy(),
    }


def sortkey(answer):
    return answer["L_hat"], answer["start"]


def estimate(problem):
    """Run every issued start using at most three module-level workers."""
    with ProcessPoolExecutor(
        max_workers=start_workers(),
        initializer=init_start_worker,
        initargs=(problem,),
    ) as executor:
        answers = list(executor.map(optimise_start, starts(problem.d)))
    return min(answers, key=sortkey), answers


def relative_field_error(fitted, exact, z, need_j, d=D, batch_size=128):
    """Relative squared Frobenius error, evaluated in bounded deterministic batches."""
    numerator = 0.0
    denominator = 0.0
    z = np.asarray(z, dtype=np.float64)
    for begin in range(0, len(z), int(batch_size)):
        batch = torch.as_tensor(z[begin : begin + int(batch_size)], dtype=torch.float64)
        estimate_field = ensemble(fitted, batch, need_j, False, d=d)[1 if need_j else 0]
        exact_field = ensemble(exact, batch, need_j, False, d=d)[1 if need_j else 0]
        numerator += float(torch.sum((estimate_field - exact_field) ** 2))
        denominator += float(torch.sum(exact_field**2))
    return float(np.sqrt(numerator / denominator)) if denominator else float("inf")


def case_cache_path(rung, case):
    return CACHE / "cases" / f"{rung}_{case}.npz"


def cache_complete(path):
    """A cache is complete only when every emitted row field and Q are present."""
    path = Path(path)
    if not path.exists():
        return False
    try:
        with np.load(path, allow_pickle=False) as data:
            return "Q" in data.files and all(column in data.files for column in CASE_COLUMNS)
    except (OSError, ValueError, KeyError):
        return False


def value_from_npz(value):
    value = np.asarray(value)
    return value.item() if value.shape == () else value


def load_case_cache(path):
    with np.load(path, allow_pickle=False) as data:
        row = {column: value_from_npz(data[column]) for column in CASE_COLUMNS}
        q_hat = np.asarray(data["Q"], dtype=np.float64)
    row["rung"] = str(row["rung"])
    row["gamma"] = float(row["gamma"])
    row["n"] = int(row["n"])
    row["r"] = int(row["r"])
    for column in ("best_start", "n_near_best"):
        row[column] = int(row[column])
    for column in (
        "err",
        "err_shape",
        "err_pos",
        "L_hat",
        "L_star",
        "fe_u_X",
        "fe_u_Y",
        "fe_J_X",
        "fe_J_Y",
    ):
        row[column] = float(row[column])
    return row, q_hat


def save_case_cache(path, row, q_hat):
    arrays = {column: np.asarray(row[column]) for column in CASE_COLUMNS}
    arrays["Q"] = np.asarray(q_hat, dtype=np.float64)
    atomic_npz(path, **arrays)


def write_cases_csv():
    rows = []
    for rung in ("A", "B"):
        for gamma, n, r, case in case_specs():
            path = case_cache_path(rung, case)
            if cache_complete(path):
                row, _ = load_case_cache(path)
                rows.append(row)
    frame = pd.DataFrame(rows, columns=CASE_COLUMNS)
    if not frame.empty:
        frame = frame.sort_values(["rung", "gamma", "n", "r"], kind="stable")
    atomic_dataframe(OUT / "cases.csv", frame)
    return frame


def write_q_hats():
    q_hats = np.full((2, len(GAMMAS), len(NS), len(REPS), D, D), np.nan, dtype=np.float64)
    present = np.zeros((2, len(GAMMAS), len(NS), len(REPS)), dtype=bool)
    for rung_index, rung in enumerate(("A", "B")):
        for gamma_index, gamma in enumerate(GAMMAS):
            for n_index, n in enumerate(NS):
                for r_index, r in enumerate(REPS):
                    path = case_cache_path(rung, case_id(gamma, n, r))
                    if cache_complete(path):
                        _, q_hat = load_case_cache(path)
                        q_hats[rung_index, gamma_index, n_index, r_index] = q_hat
                        present[rung_index, gamma_index, n_index, r_index] = True
    atomic_npz(
        ARTIFACTS / "Q_hats.npz",
        Q_hats=q_hats,
        completed=present,
        rungs=np.asarray(("A", "B")),
        gammas=np.asarray(GAMMAS, dtype=np.float64),
        ns=np.asarray(NS, dtype=np.int64),
        reps=np.asarray(REPS, dtype=np.int64),
    )


def emit_finished_case(row):
    """Append exactly one JSON progress record after a case's atomic checkpoint."""
    progress = {
        "event": "case_finished",
        "rung": row["rung"],
        "gamma": float(row["gamma"]),
        "n": int(row["n"]),
        "r": int(row["r"]),
        "err": float(row["err"]),
        "err_shape": float(row["err_shape"]),
        "L_hat": float(row["L_hat"]),
        "L_star": float(row["L_star"]),
    }
    line = json.dumps(progress, sort_keys=True)
    with (OUT / "run.log").open("a") as handle:
        handle.write(line + "\n")
    print(line, flush=True)


def estimate_one_case(rung, gamma, n, r, case, landscape):
    path = case_cache_path(rung, case)
    if cache_complete(path):
        return False

    x, y, _, _, q_star = sample_case(gamma, n, r, landscape)
    true_x, true_y = exact_models(gamma, q_star, landscape)
    if rung == "A":
        mx, my = true_x, true_y
    elif rung == "B":
        mx, my = get_models(case)
    else:
        raise ValueError(rung)

    problem = Problem(mx, my, x, y, r, d=D)
    best, answers = estimate(problem)
    q_hat = np.asarray(best["Q"], dtype=np.float64)
    q_star_tensor = torch.as_tensor(q_star, dtype=torch.float64)
    l_star, _, _ = problem.loss(q_star_tensor, False)

    p_shape = np.diag(np.asarray([0.0] * 3 + [1.0] * 7, dtype=np.float64))
    p_pos = np.diag(np.asarray([1.0] * 3 + [0.0] * 7, dtype=np.float64))
    delta = q_hat - q_star
    if rung == "B":
        fe_u_x = relative_field_error(mx, true_x, x, False, d=D)
        fe_u_y = relative_field_error(my, true_y, y, False, d=D)
        # Problem's xa/yb are exactly the pre-registered 256-cell H subsets.
        fe_j_x = relative_field_error(mx, true_x, problem.xa.cpu().numpy(), True, d=D)
        fe_j_y = relative_field_error(my, true_y, problem.yb.cpu().numpy(), True, d=D)
    else:
        fe_u_x = float("nan")
        fe_u_y = float("nan")
        fe_j_x = float("nan")
        fe_j_y = float("nan")

    row = {
        "rung": rung,
        "gamma": float(gamma),
        "n": int(n),
        "r": int(r),
        "err": float(np.linalg.norm(delta, ord="fro") / np.sqrt(20.0)),
        "err_shape": float(np.linalg.norm(delta @ p_shape, ord="fro") / np.sqrt(14.0)),
        "err_pos": float(np.linalg.norm(delta @ p_pos, ord="fro") / np.sqrt(6.0)),
        "L_hat": float(best["L_hat"]),
        "L_star": float(l_star.detach()),
        "best_start": int(best["start"]),
        "n_near_best": int(sum(answer["L_hat"] <= 1.01 * best["L_hat"] for answer in answers)),
        "fe_u_X": float(fe_u_x),
        "fe_u_Y": float(fe_u_y),
        "fe_J_X": float(fe_j_x),
        "fe_J_Y": float(fe_j_y),
    }
    save_case_cache(path, row, q_hat)
    emit_finished_case(row)
    return True


def estimate_stage():
    """Estimate all rung-A cases followed by all rung-B cases, resumably."""
    make_meta()
    landscape = ensure_landscape()
    new_cases = 0
    for rung in ("A", "B"):
        for gamma, n, r, case in case_specs():
            new_cases += int(estimate_one_case(rung, gamma, n, r, case, landscape))
            # Both public summaries are reconstructed from only atomically complete cases.
            write_cases_csv()
            write_q_hats()
    completed = write_cases_csv()
    write_q_hats()
    print(
        json.dumps(
            {"completed_cases": int(len(completed)), "new_cases": int(new_cases)},
            sort_keys=True,
        ),
        flush=True,
    )


def import_o1g_by_absolute_path():
    """Import O1g's own runner solely for its reference Problem/loss implementation."""
    source = O1G_SCRIPT.resolve()
    if not source.is_file():
        raise FileNotFoundError(
            f"O1g runner required for the copy check is absent: {source}"
        )
    name = "o1g_run_o1k_copy_check"
    spec = importlib.util.spec_from_file_location(name, source)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not import O1g runner at {source}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def o1g_test_r_inputs(o1g):
    """Reproduce only O1g's cached R_n2759_r0 selection without calling its writers."""
    z = np.load(o1g.LAT1).astype(np.float32)
    if len(z) != 5518:
        raise RuntimeError(f"O1g LAT1 must contain 5518 cells, found {len(z)}")
    n = 2759
    rng = np.random.default_rng(20261030)
    permutation = rng.permutation(5518)
    q_star = np.asarray(ortho_group.rvs(o1g.D, random_state=20261031), dtype=np.float64)
    x = z[np.sort(permutation[:n])]
    y = (z[np.sort(permutation[n : 2 * n])].astype(np.float64) @ q_star.T).astype(np.float32)
    return x, y, q_star


def o1g_models_readonly(o1g, case):
    """Load the O1g cached ensemble without calling get_models/save_sym."""
    mx = []
    my = []
    read_paths = [Path(o1g.__file__).resolve(), Path(o1g.LAT1).resolve()]
    for seed in o1g.SEEDS20:
        for side, bucket in (("X", mx), ("Y", my)):
            path = Path(o1g.local_fit(case, side, seed)).resolve()
            if not path.is_file():
                raise RuntimeError(
                    "O1k copy check requires O1g's cached Test R fit, but it is missing: "
                    f"{path}"
                )
            saved = o1g.load_saved(path)
            bucket.append(o1g.build_model(saved))
            read_paths.append(path)
    unique = []
    seen = set()
    for path in read_paths:
        if path not in seen:
            seen.add(path)
            unique.append(path)
    return mx, my, unique


def copy_check():
    """Compare copied E5b loss with O1g itself at I and its true Test-R Q*."""
    o1g = import_o1g_by_absolute_path()
    if int(o1g.D) != 30:
        raise RuntimeError(f"Expected O1g dimension 30 for copy check, found {o1g.D}")
    case = "R_n2759_r0"
    x, y, q_star = o1g_test_r_inputs(o1g)
    mx, my, paths = o1g_models_readonly(o1g, case)
    copied = Problem(mx, my, x, y, 0, d=int(o1g.D))
    reference = o1g.Problem(mx, my, x, y, 0)

    comparisons = {}
    max_difference = 0.0
    for label, q in (("I", np.eye(o1g.D, dtype=np.float64)), ("Q_star", q_star)):
        q_tensor = torch.as_tensor(q, dtype=torch.float64)
        copied_value = float(copied.loss(q_tensor, False)[0].detach())
        o1g_value = float(reference.loss(q_tensor, False)[0].detach())
        difference = abs(copied_value - o1g_value)
        comparisons[label] = {
            "copied_loss": copied_value,
            "o1g_loss": o1g_value,
            "absolute_difference": difference,
        }
        max_difference = max(max_difference, difference)
    return {
        "case": case,
        "D": int(o1g.D),
        "loss_copy_check": float(max_difference),
        "losses": comparisons,
        "readonly_files": [
            {"path": str(path), "sha256": sha(path)} for path in paths
        ],
    }


def check_stage():
    """Run all issued non-fit checks and create the landscape artifact."""
    make_meta()
    ensure_landscape()
    record = copy_check()
    record["start_overlap"] = int(start_overlap())
    atomic_json(OUT / "check.json", record)
    print(
        json.dumps(
            {
                "loss_copy_check": record["loss_copy_check"],
                "start_overlap": record["start_overlap"],
            },
            sort_keys=True,
        ),
        flush=True,
    )


def require_completed_inputs():
    """Read report inputs only after the full issued estimate has completed."""
    check_path = OUT / "check.json"
    cases_path = OUT / "cases.csv"
    symmetry_path = OUT / "fit_symmetry.csv"
    missing = [str(path) for path in (check_path, cases_path, symmetry_path) if not path.is_file()]
    if missing:
        raise RuntimeError("--stage report requires completed earlier stages; missing: " + ", ".join(missing))
    check = json.loads(check_path.read_text())
    cases = pd.read_csv(cases_path)
    symmetry = pd.read_csv(symmetry_path)
    if len(cases) != 72:
        raise RuntimeError(f"--stage report requires 72 completed cases, found {len(cases)}")
    if len(symmetry) != 360:
        raise RuntimeError(f"--stage report requires 360 fit rows, found {len(symmetry)}")
    if not (ARTIFACTS / "Q_hats.npz").is_file():
        raise RuntimeError("--stage report requires artifacts/Q_hats.npz")
    return check, cases, symmetry


def gamma_mask(frame, gamma):
    return np.isclose(frame["gamma"].to_numpy(dtype=float), float(gamma), atol=1e-12, rtol=0.0)


def scalar_median(frame, column, *, rung=None, gamma=None, n=None):
    subset = frame
    if rung is not None:
        subset = subset.loc[subset["rung"] == rung]
    if gamma is not None:
        subset = subset.loc[gamma_mask(subset, gamma)]
    if n is not None:
        subset = subset.loc[subset["n"] == int(n)]
    if len(subset) != 3:
        raise RuntimeError(
            f"Expected three replicate rows for median {column}, got {len(subset)}"
        )
    return float(subset[column].median())


def gated_values(check, cases, symmetry):
    """Compute exactly the sixteen pre-registered gates and no other gate."""
    a = cases.loc[cases["rung"] == "A"]
    b = cases.loc[cases["rung"] == "B"]
    err_a_medians = [
        scalar_median(cases, "err", rung="A", gamma=gamma, n=n)
        for gamma in (1.15, 1.33, 1.6)
        for n in NS
    ]
    b_positive = b.loc[~gamma_mask(b, 1.0)]
    if len(b_positive) != 27:
        raise RuntimeError(f"Expected 27 positive-anisotropy rung-B cases, found {len(b_positive)}")
    b_n3000 = b.loc[b["n"] == 3000]
    if len(b_n3000) != 12:
        raise RuntimeError(f"Expected 12 rung-B n=3000 cases, found {len(b_n3000)}")
    fe_u = np.concatenate(
        [b_n3000["fe_u_X"].to_numpy(dtype=float), b_n3000["fe_u_Y"].to_numpy(dtype=float)]
    )
    fe_j = np.concatenate(
        [b_n3000["fe_J_X"].to_numpy(dtype=float), b_n3000["fe_J_Y"].to_numpy(dtype=float)]
    )
    if not np.isfinite(fe_u).all() or not np.isfinite(fe_j).all():
        raise RuntimeError("Rung-B field-fit errors must be finite before reporting")

    values = {
        "n_cases": int(len(cases)),
        "start_overlap": int(check["start_overlap"]),
        "loss_copy_check": float(check["loss_copy_check"]),
        "symmetry_max": float(symmetry["median_asym"].max()),
        "negctrl_err_shape_A": scalar_median(cases, "err_shape", rung="A", gamma=1.0, n=3000),
        "errA_max": float(max(err_a_medians)),
        "errB_g115_n3000": scalar_median(cases, "err", rung="B", gamma=1.15, n=3000),
        "errB_g133_n3000": scalar_median(cases, "err", rung="B", gamma=1.33, n=3000),
        "errB_g160_n3000": scalar_median(cases, "err", rung="B", gamma=1.6, n=3000),
        "errB_g115_n10000": scalar_median(cases, "err", rung="B", gamma=1.15, n=10000),
        "errB_g133_n10000": scalar_median(cases, "err", rung="B", gamma=1.33, n=10000),
        "errB_g160_n10000": scalar_median(cases, "err", rung="B", gamma=1.6, n=10000),
        "frac_hat_below_star_B": float(
            np.mean(
                b_positive["L_hat"].to_numpy(dtype=float)
                <= b_positive["L_star"].to_numpy(dtype=float)
            )
        ),
        "fe_u_3000": float(np.median(fe_u)),
        "fe_J_3000": float(np.median(fe_j)),
    }
    validity_failures = []
    if values["n_cases"] != 72:
        validity_failures.append("n_cases != 72")
    if values["start_overlap"] != 0:
        validity_failures.append("start_overlap != 0")
    if values["loss_copy_check"] > 1e-10:
        validity_failures.append("loss_copy_check > 1e-10")
    if values["symmetry_max"] > 1e-4:
        validity_failures.append("symmetry_max > 1e-4")
    if values["negctrl_err_shape_A"] < 0.30:
        validity_failures.append("negctrl_err_shape_A < 0.30")

    if validity_failures:
        # "REVISE" is the explicit first branch of the issued decision rule.
        outcome = "REVISE"
    elif values["errA_max"] > 0.10:
        outcome = "test_fault"
    else:
        gamma_min = next(
            (
                gamma
                for gamma, key in (
                    (1.15, "errB_g115_n3000"),
                    (1.33, "errB_g133_n3000"),
                    (1.6, "errB_g160_n3000"),
                )
                if values[key] <= 0.10
            ),
            None,
        )
        if gamma_min in (1.15, 1.33):
            outcome = "pbmc_range_ok"
        elif gamma_min == 1.6:
            outcome = "anisotropy_limited"
        else:
            outcome = "fit_limited"
    values["o1k_outcome"] = outcome
    if tuple(values) != GATED_KEYS:
        raise AssertionError("gated key order changed")
    return values, validity_failures


def ladder_figure(cases):
    """Plot issued median/replicate ladders for error and shape-subspace error."""
    fig, axes = plt.subplots(2, len(NS), figsize=(4.7 * len(NS), 7.0), sharex=True)
    colors = {"A": "C0", "B": "C1"}
    labels = {"A": "rung A (exact)", "B": "rung B (fitted)"}
    x = np.arange(len(GAMMAS), dtype=float)
    for col, n in enumerate(NS):
        for row_index, metric in enumerate(("err", "err_shape")):
            axis = axes[row_index, col]
            for rung in ("A", "B"):
                medians = []
                for index, gamma in enumerate(GAMMAS):
                    subset = cases.loc[
                        (cases["rung"] == rung)
                        & (cases["n"] == int(n))
                        & gamma_mask(cases, gamma)
                    ]
                    values = subset[metric].to_numpy(dtype=float)
                    axis.scatter(
                        np.full(len(values), x[index]),
                        values,
                        color=colors[rung],
                        alpha=0.55,
                        s=24,
                        zorder=2,
                    )
                    medians.append(float(np.median(values)))
                axis.plot(x, medians, marker="o", color=colors[rung], label=labels[rung], zorder=3)
            axis.axhline(0.10, color="black", linestyle="--", linewidth=1.0, label="0.10" if col == 0 and row_index == 0 else None)
            axis.set_xticks(x)
            axis.set_xticklabels([str(gamma) for gamma in GAMMAS])
            axis.set_xlabel("anisotropy $\\gamma$")
            axis.set_title(f"n = {n}")
            axis.grid(axis="y", alpha=0.25)
            if col == 0:
                axis.set_ylabel("rotation error" if metric == "err" else "shape-subspace error")
            if row_index == 0 and col == 0:
                axis.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT / "fig_ladder.png", dpi=180)
    plt.close(fig)


def umap_coordinates(rep):
    """O1g's seeded scanpy UMAP settings, applied to one synthetic view."""
    import anndata as ad
    import scanpy as sc

    adata = ad.AnnData(X=np.zeros((len(rep), 1), dtype=np.float32))
    adata.obsm["view"] = np.asarray(rep, dtype=np.float64)
    sc.pp.neighbors(adata, n_neighbors=15, use_rep="view", random_state=20261014)
    sc.tl.umap(adata, min_dist=0.5, random_state=20261014)
    return np.asarray(adata.obsm["X_umap"], dtype=np.float64)


def mixing(rep, dataset_labels):
    neighbors = NearestNeighbors(n_neighbors=31).fit(rep).kneighbors(rep, return_distance=False)[:, 1:]
    return float(np.mean(np.mean(dataset_labels[neighbors] != dataset_labels[:, None], axis=1)) / 0.5)


def umap_figure():
    """Create the prescribed gamma=1.33, n=3000, r=0 rung-B illustrative views."""
    gamma, n, r = 1.33, 3000, 0
    cache_path = case_cache_path("B", case_id(gamma, n, r))
    if not cache_complete(cache_path):
        raise RuntimeError(f"Missing rung-B cache required for UMAP view: {cache_path}")
    _, q_hat = load_case_cache(cache_path)
    landscape = landscape_arrays()
    x, y, x_labels, y_labels, q_star = sample_case(gamma, n, r, landscape)
    datasets = np.asarray(["X"] * len(x) + ["Y"] * len(y), dtype=str)
    wells = np.concatenate([x_labels, y_labels]).astype(int)
    views = (
        ("no rotation", np.eye(D, dtype=np.float64)),
        ("rung B estimate", q_hat),
        ("known answer", q_star),
    )
    fig, axes = plt.subplots(3, 2, figsize=(10.5, 14.0))
    palette = plt.get_cmap("tab10")
    for index, (name, q) in enumerate(views):
        rep = np.vstack([x @ q.T, y])
        coords = umap_coordinates(rep)
        mix = mixing(rep, datasets)
        sil = float(silhouette_score(rep, wells, metric="euclidean"))

        by_copy = axes[index, 0]
        for label, color in (("X", "C0"), ("Y", "C1")):
            mask = datasets == label
            by_copy.scatter(coords[mask, 0], coords[mask, 1], s=4, alpha=0.7, color=color, label=label, linewidths=0)
        by_copy.set_title(f"{name}: copy; mix={mix:.3f}, sil={sil:.3f}")
        by_copy.legend(markerscale=2, fontsize=8)
        by_copy.set_xticks([])
        by_copy.set_yticks([])

        by_well = axes[index, 1]
        for well in range(K):
            mask = wells == well
            by_well.scatter(coords[mask, 0], coords[mask, 1], s=4, alpha=0.7, color=palette(well), label=f"well {well}", linewidths=0)
        by_well.set_title(f"{name}: wells; mix={mix:.3f}, sil={sil:.3f}")
        by_well.legend(markerscale=2, fontsize=8, ncol=2)
        by_well.set_xticks([])
        by_well.set_yticks([])
    fig.tight_layout()
    fig.savefig(OUT / "fig_umap.png", dpi=180)
    plt.close(fig)


def write_lockfile():
    """Capture the active launch environment atomically for the report."""
    path = OUT / "requirements.lock"
    tmp = atomic_path(path)
    try:
        with tmp.open("w") as handle:
            subprocess.run([sys.executable, "-m", "pip", "freeze"], check=True, stdout=handle)
        os.replace(tmp, path)
    finally:
        tmp.unlink(missing_ok=True)


def git_commit():
    return subprocess.check_output(
        ["git", "-C", str(ROOT), "rev-parse", "HEAD"], text=True
    ).strip()


def checkpoint_stage(stage):
    """Commit and push tracked O1k stage products when invoked through --stage all.

    Caches and detached-run files remain excluded by .gitignore.  This is an
    execution checkpoint only: it does not alter any generated quantity.
    """
    if os.environ.get("O1K_AUTO_CHECKPOINT", "1") != "1":
        return
    subprocess.run(["git", "-C", str(ROOT), "add", "reports/O1k"], check=True)
    changed = subprocess.run(
        ["git", "-C", str(ROOT), "diff", "--cached", "--quiet"]
    ).returncode
    if changed == 0:
        return
    if changed != 1:
        raise RuntimeError(f"git diff --cached failed after O1k {stage} stage: {changed}")
    subprocess.run(
        ["git", "-C", str(ROOT), "commit", "-m", f"task: O1k {stage} stage"],
        check=True,
    )
    subprocess.run(["git", "-C", str(ROOT), "push", "origin", "task/O1k"], check=True)


def regeneration_commands():
    command = 'cd "$HOME/Documents/scopalign-manus-O1k" && .venv/bin/python reports/O1k/scripts/o1k_run.py'
    lock_command = 'cd "$HOME/Documents/scopalign-manus-O1k" && .venv/bin/python -m pip freeze > reports/O1k/requirements.lock'
    result = {
        "artifacts/landscape.npz": command + " --stage check",
        "check.json": command + " --stage check",
        "fit_symmetry.csv": command + " --stage fits",
        "cases.csv": command + " --stage estimate",
        "artifacts/Q_hats.npz": command + " --stage estimate",
        "fig_ladder.png": command + " --stage report",
        "fig_umap.png": command + " --stage report",
        "requirements.lock": lock_command,
        "report.json": command + " --stage report",
        "scripts/o1k_run.py": command + " --stage all",
    }
    for key in GATED_KEYS:
        if key in ("start_overlap", "loss_copy_check"):
            result[key] = command + " --stage check"
        elif key == "symmetry_max":
            result[key] = command + " --stage fits"
        else:
            result[key] = command + " --stage report"
    return result


def checksum_entries(paths):
    return [{"path": str(path), "sha256": sha(OUT / path)} for path in paths]


def report_stage():
    """Build figures, lockfile, and the final checker-facing O1k report."""
    meta = make_meta()
    check, cases, symmetry = require_completed_inputs()
    gated, validity_failures = gated_values(check, cases, symmetry)
    ladder_figure(cases)
    umap_figure()
    write_lockfile()

    deliverable_paths = (
        "requirements.lock",
        "cases.csv",
        "fit_symmetry.csv",
        "fig_ladder.png",
        "fig_umap.png",
        "scripts/o1k_run.py",
    )
    artifact_paths = (
        "artifacts/landscape.npz",
        "artifacts/Q_hats.npz",
        "check.json",
    )
    report = {
        "task_id": "O1k",
        "executor": "manus",
        "branch": "task/O1k",
        "commit": git_commit(),
        "seed": 20261014,
        "started": meta["started"],
        "finished": datetime.now(timezone.utc).isoformat(),
        # Status describes execution completeness.  A failed validity decision
        # remains recorded solely in the pre-registered o1k_outcome gate.
        "status": "complete",
        "gated": gated,
        "deliverables": checksum_entries(deliverable_paths),
        "artifacts": checksum_entries(artifact_paths),
        "regeneration": regeneration_commands(),
        "deviations": [],
        "deviations_requested": [],
        "observations": [],
        "known_failures": KNOWN_FAILURES,
        "lockfile_sha256": sha(OUT / "requirements.lock"),
    }
    if set(report["gated"]) != set(GATED_KEYS) or len(report["gated"]) != 16:
        raise AssertionError("O1k report must gate exactly the sixteen issued keys")
    atomic_json(OUT / "report.json", report)
    print(json.dumps(gated, sort_keys=True), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--stage", required=True, choices=("check", "fits", "estimate", "report", "all")
    )
    stage = parser.parse_args().stage
    if stage in ("check", "all"):
        check_stage()
        if stage == "all":
            checkpoint_stage("check")
    if stage in ("fits", "all"):
        fit_stage()
        if stage == "all":
            checkpoint_stage("fits")
    if stage in ("estimate", "all"):
        estimate_stage()
        if stage == "all":
            checkpoint_stage("estimate")
    if stage in ("report", "all"):
        report_stage()
        if stage == "all":
            checkpoint_stage("report")


if __name__ == "__main__":
    main()
