#!/usr/bin/env python
"""Reproducible runner for O1l: O1k E5b with axis-sign starts and a basin diagnostic.

Run from ~/Documents/scopalign-manus-O1l with --stage check, estimate, basin,
report, or all.  O1k's landscape, case cache, fitted models, results and runner
are strictly read-only inputs at the absolute O1k path declared below.
"""

import argparse
import hashlib
import json
import multiprocessing
import os
import subprocess
import sys
import tempfile
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime, timezone
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import torch
from scipy.linalg import expm
from scipy.stats import ortho_group
from sklearn.metrics import silhouette_score
from sklearn.neighbors import NearestNeighbors

from scjdo2.energy_drift import _EnergyMLP

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "reports" / "O1l"
CACHE = OUT / "cache"
ARTIFACTS = OUT / "artifacts"
O1K = Path.home() / "Documents" / "scopalign-manus-O1k" / "reports" / "O1k"
O1K_CACHE = O1K / "cache"
O1K_ARTIFACTS = O1K / "artifacts"

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
    "inputs_match",
    "l_star_repro_max",
    "negctrl_err_shape_A_s",
    "errA_s_max",
    "errB_s_g115_n3000",
    "errB_s_g133_n3000",
    "errB_s_g160_n3000",
    "errB_s_g115_n10000",
    "errB_s_g133_n10000",
    "errB_s_g160_n10000",
    "basin_radius_A",
    "basin_radius_B",
    "o1l_outcome",
)
CASE_COLUMNS = (
    "rung", "gamma", "n", "r", "err", "err_shape", "err_pos", "L_hat",
    "L_star", "screen_rank_best",
) + tuple(f"cov_eig_X_{j}" for j in range(D)) + tuple(f"cov_eig_Y_{j}" for j in range(D))
WORKER_PROBLEM = None


def chunks(handle):
    while True:
        block = handle.read(1 << 20)
        if not block:
            return
        yield block


def sha(path):
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(path)
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in chunks(handle):
            digest.update(block)
    return digest.hexdigest()


def atomic_path(path, suffix=".tmp"):
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


def make_meta():
    OUT.mkdir(parents=True, exist_ok=True)
    CACHE.mkdir(parents=True, exist_ok=True)
    ARTIFACTS.mkdir(parents=True, exist_ok=True)
    gitignore = OUT / ".gitignore"
    ignored = "cache/\n__pycache__/\nrun.log\nrun.pid\n.cache/\n.mplconfig/\n"
    if not gitignore.exists() or gitignore.read_text() != ignored:
        gitignore.write_text(ignored)
    meta_path = OUT / "run_meta.json"
    if not meta_path.exists():
        atomic_json(meta_path, {"started": datetime.now(timezone.utc).isoformat()})
    return json.loads(meta_path.read_text())


def gamma_tag(gamma):
    return f"g{gamma:.2f}".replace(".", "p")


def case_id(gamma, n, r):
    return f"{gamma_tag(gamma)}_n{int(n)}_r{int(r)}"


def case_specs():
    for gamma in GAMMAS:
        for n in NS:
            for r in REPS:
                yield float(gamma), int(n), int(r), case_id(gamma, n, r)


def landscape_arrays():
    """O1k's fixed landscape, with exactly its random draw order."""
    rng = np.random.default_rng(20264200)
    means = np.zeros((K, D), dtype=np.float64)
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
    return {"means": means, "rotations": rotations, "gammas": np.asarray(GAMMAS), "sigma2": sigma2, "q_stars": q_stars}


def load_o1k_landscape():
    path = O1K_ARTIFACTS / "landscape.npz"
    if not path.is_file():
        raise FileNotFoundError(f"Required read-only O1k landscape is absent: {path}")
    with np.load(path, allow_pickle=False) as data:
        result = {key: np.asarray(data[key]) for key in data.files}
    expected = landscape_arrays()
    if set(result) != set(expected) or any(not np.array_equal(result[key], expected[key]) for key in expected):
        raise RuntimeError("O1k landscape does not equal its prescribed deterministic landscape")
    return result


def sigma2_for_gamma(gamma):
    return np.asarray(gamma ** (5.5 - np.arange(1, D + 1)), dtype=np.float64)


def sample_case(gamma, n, r, landscape):
    """O1k's exact cells: fixed independent label then normal draws per copy."""
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
        return means[labels] + np.einsum("nij,nj->ni", rotations[labels], eps * sigma[None, :])

    x = draw(x_labels, x_eps)
    y_prime = draw(y_labels, y_eps)
    y = y_prime @ q_star.T
    return x, y, x_labels, y_labels, q_star


class ExactEnergy(torch.nn.Module):
    """O1k's exact float64 Gaussian-mixture energy in a row coordinate frame."""

    def __init__(self, means, rotations, sigma2, coordinate_to_base=None):
        super().__init__()
        means = np.asarray(means, dtype=np.float64)
        rotations = np.asarray(rotations, dtype=np.float64)
        sigma2 = np.asarray(sigma2, dtype=np.float64)
        if coordinate_to_base is None:
            coordinate_to_base = np.eye(means.shape[1], dtype=np.float64)
        precision = np.einsum("kij,j,klj->kil", rotations, 1.0 / sigma2, rotations)
        precision = 0.5 * (precision + np.swapaxes(precision, 1, 2))
        log_norm = -0.5 * (means.shape[1] * np.log(2.0 * np.pi) + np.sum(np.log(sigma2))) - np.log(means.shape[0])
        self.register_buffer("means", torch.as_tensor(means, dtype=torch.float64))
        self.register_buffer("precision", torch.as_tensor(precision, dtype=torch.float64))
        self.register_buffer("coordinate_to_base", torch.as_tensor(np.asarray(coordinate_to_base, dtype=np.float64), dtype=torch.float64))
        self.register_buffer("log_norm", torch.as_tensor(log_norm, dtype=torch.float64))

    def forward(self, x):
        z = x @ self.coordinate_to_base
        delta = z[:, None, :] - self.means[None, :, :]
        quadratic = torch.einsum("nki,kij,nkj->nk", delta, self.precision, delta)
        return -torch.logsumexp(self.log_norm - 0.5 * quadratic, dim=1)


def exact_models(gamma, q_star, landscape):
    base = ExactEnergy(landscape["means"], landscape["rotations"], sigma2_for_gamma(gamma)).eval()
    rotated = ExactEnergy(landscape["means"], landscape["rotations"], sigma2_for_gamma(gamma), coordinate_to_base=q_star).eval()
    return [base], [rotated]


def fields(model, x, need_j, track, d=None):
    if d is None:
        d = int(x.shape[1])
    z = x.clone().requires_grad_(True)
    e = model(z)
    g = torch.autograd.grad(e.sum(), z, create_graph=(need_j or track), retain_graph=(need_j or track))[0]
    u = -g
    if not need_j:
        return (u, None) if track else (u.detach(), None)
    rows = []
    for j in range(int(d)):
        rows.append(-torch.autograd.grad(g[:, j].sum(), z, create_graph=track, retain_graph=(track or j < int(d) - 1))[0])
    jacobian = torch.stack(rows, dim=1)
    return (u, jacobian) if track else (u.detach(), jacobian.detach())


def ensemble(models, x, need_j, track, d=None):
    us, js = [], []
    for model in models:
        u, jacobian = fields(model, x, need_j, track, d=d)
        us.append(u)
        if need_j:
            js.append(jacobian)
    return torch.stack(us).mean(0), (torch.stack(js).mean(0) if need_j else None)


def build_model(saved, d=D):
    model = _EnergyMLP(int(d), int(saved["config"]["hidden"]), int(saved["config"]["depth"]))
    model.load_state_dict(saved["state_dict"])
    model = model.to(dtype=torch.float64, device="cpu").eval()
    for parameter in model.parameters():
        parameter.requires_grad_(False)
    return model


def get_models(case):
    """Load O1k's five cached models per side, with no writes or fitting."""
    mx, my = [], []
    for seed in FIRST5:
        for side, bucket in (("X", mx), ("Y", my)):
            path = O1K_CACHE / "drifts" / f"{case}_{side}_{int(seed)}.pt"
            if not path.is_file():
                raise FileNotFoundError(f"Required read-only O1k fit is absent: {path}")
            bucket.append(build_model(torch.load(path, map_location="cpu", weights_only=False), D))
    return mx, my


class Problem:
    """O1k's unchanged E5b D + H loss and O1k's H subset selection."""

    def __init__(self, mx, my, x, y, r, d=D):
        self.mx, self.my, self.d = mx, my, int(d)
        self.x = torch.as_tensor(np.asarray(x, dtype=np.float64), dtype=torch.float64)
        self.y = torch.as_tensor(np.asarray(y, dtype=np.float64), dtype=torch.float64)
        rng = np.random.default_rng(20261040 + int(r))
        self.xa = self.x[rng.choice(len(x), size=min(256, len(x)), replace=False)]
        self.yb = self.y[rng.choice(len(y), size=min(256, len(y)), replace=False)]
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
        dl = (torch.sum((uy - self.ux @ q.T) ** 2) / self.du + torch.sum((ux - self.uy @ q) ** 2) / self.dv) / 2
        _, jy = ensemble(self.my, self.xa @ q.T, True, track, d=self.d)
        _, jx = ensemble(self.mx, self.yb @ q, True, track, d=self.d)
        qj = torch.einsum("ij,njk,lk->nil", q, self.jx, q)
        qtj = torch.einsum("ji,njk,kl->nil", q, self.jy, q)
        hl = (torch.sum((jy - qj) ** 2) / self.dx + torch.sum((jx - qtj) ** 2) / self.dy) / 2
        return dl + hl, dl, hl


class Closure:
    """O1k's unchanged LBFGS manifold closure."""

    def __init__(self, problem, q0, skew):
        self.problem, self.q0, self.skew = problem, q0, skew

    def __call__(self):
        self.skew.grad = None
        q = self.q0 @ torch.matrix_exp(self.skew - self.skew.T)
        value, _, _ = self.problem.loss(q, True)
        value.backward()
        return value


def start_workers():
    return min(3, max(1, (os.cpu_count() or 1) - 1))


def init_start_worker(problem):
    global WORKER_PROBLEM
    torch.set_num_threads(1)
    WORKER_PROBLEM = problem


def optimise_start(item):
    """O1k's exact LBFGS configuration, applied to one issued start."""
    torch.set_num_threads(1)
    start, q0_numpy = item
    q0 = torch.as_tensor(q0_numpy, dtype=torch.float64)
    skew = torch.zeros((WORKER_PROBLEM.d, WORKER_PROBLEM.d), dtype=torch.float64, requires_grad=True)
    optimiser = torch.optim.LBFGS([skew], lr=1.0, max_iter=500, tolerance_grad=1e-10, tolerance_change=1e-12, line_search_fn="strong_wolfe")
    optimiser.step(Closure(WORKER_PROBLEM, q0, skew))
    q = (q0 @ torch.matrix_exp(skew - skew.T)).detach()
    value, dl, hl = WORKER_PROBLEM.loss(q, False)
    return {"start": int(start), "L_hat": float(value.detach()), "D": float(dl.detach()), "H": float(hl.detach()), "Q": q.cpu().numpy()}


def refine(problem, items):
    """Refine eight selected starts with three forked workers (at most three)."""
    # The default macOS spawn context serializes this torch-backed Problem through
    # PyTorch's prohibited shared-memory helper.  Fork preserves the same loss,
    # selected starts, and LBFGS settings without that helper.
    with ProcessPoolExecutor(
        max_workers=start_workers(),
        mp_context=multiprocessing.get_context("fork"),
        initializer=init_start_worker,
        initargs=(problem,),
    ) as executor:
        answers = list(executor.map(optimise_start, items))
    return min(answers, key=lambda answer: (answer["L_hat"], answer["start"])), answers


def optimise_single(problem, q0):
    """The same LBFGS routine for a single basin-diagnostic start."""
    global WORKER_PROBLEM
    WORKER_PROBLEM = problem
    return optimise_start((0, np.asarray(q0, dtype=np.float64)))


def case_cache_path(rung, case):
    return CACHE / "cases" / f"{rung}_{case}.npz"


def basin_cache_path(rung, case):
    return CACHE / "basin" / f"{rung}_{case}.npz"


def value_from_npz(value):
    value = np.asarray(value)
    return value.item() if value.shape == () else value


def case_cache_complete(path):
    path = Path(path)
    if not path.is_file():
        return False
    try:
        with np.load(path, allow_pickle=False) as data:
            required = set(CASE_COLUMNS) | {"Q", "screen_index", "screen_L"}
            return required.issubset(data.files) and data["Q"].shape == (D, D) and data["screen_L"].shape == (1 << D,)
    except (OSError, ValueError, KeyError):
        return False


def load_case_cache(path):
    with np.load(path, allow_pickle=False) as data:
        row = {column: value_from_npz(data[column]) for column in CASE_COLUMNS}
        q_hat = np.asarray(data["Q"], dtype=np.float64)
        screen_index = np.asarray(data["screen_index"], dtype=np.int64)
        screen_l = np.asarray(data["screen_L"], dtype=np.float64)
    row["rung"] = str(row["rung"])
    row["gamma"] = float(row["gamma"])
    row["n"] = int(row["n"])
    row["r"] = int(row["r"])
    row["screen_rank_best"] = int(row["screen_rank_best"])
    for column in CASE_COLUMNS:
        if column.startswith("cov_eig") or column in ("err", "err_shape", "err_pos", "L_hat", "L_star"):
            row[column] = float(row[column])
    return row, q_hat, screen_index, screen_l


def save_case_cache(path, row, q_hat, screen_index, screen_l):
    arrays = {column: np.asarray(row[column]) for column in CASE_COLUMNS}
    arrays["Q"] = np.asarray(q_hat, dtype=np.float64)
    arrays["screen_index"] = np.asarray(screen_index, dtype=np.int64)
    arrays["screen_L"] = np.asarray(screen_l, dtype=np.float64)
    atomic_npz(path, **arrays)


def write_cases_csv():
    rows = []
    for rung in ("A", "B"):
        for gamma, n, r, case in case_specs():
            path = case_cache_path(rung, case)
            if case_cache_complete(path):
                row, _, _, _ = load_case_cache(path)
                rows.append(row)
    frame = pd.DataFrame(rows, columns=CASE_COLUMNS)
    if not frame.empty:
        frame = frame.sort_values(["rung", "gamma", "n", "r"], kind="stable")
    atomic_dataframe(OUT / "cases.csv", frame)
    return frame


def write_screen_csv():
    rows = []
    for rung in ("A", "B"):
        for gamma, n, r, case in case_specs():
            path = case_cache_path(rung, case)
            if not case_cache_complete(path):
                continue
            _, _, indices, losses = load_case_cache(path)
            rows.extend({"rung": rung, "gamma": gamma, "n": n, "r": r, "case": case, "s_index": int(i), "L": float(losses[int(i)])} for i in indices)
    frame = pd.DataFrame(rows, columns=("rung", "gamma", "n", "r", "case", "s_index", "L"))
    if not frame.empty:
        frame = frame.sort_values(["rung", "gamma", "n", "r", "s_index"], kind="stable")
    atomic_dataframe(ARTIFACTS / "screen.csv", frame)
    return frame


def write_q_hats():
    q_hats = np.full((2, len(GAMMAS), len(NS), len(REPS), D, D), np.nan, dtype=np.float64)
    present = np.zeros((2, len(GAMMAS), len(NS), len(REPS)), dtype=bool)
    for ri, rung in enumerate(("A", "B")):
        for gi, gamma in enumerate(GAMMAS):
            for ni, n in enumerate(NS):
                for pi, r in enumerate(REPS):
                    path = case_cache_path(rung, case_id(gamma, n, r))
                    if case_cache_complete(path):
                        _, q_hat, _, _ = load_case_cache(path)
                        q_hats[ri, gi, ni, pi] = q_hat
                        present[ri, gi, ni, pi] = True
    atomic_npz(ARTIFACTS / "Q_hats.npz", Q_hats=q_hats, completed=present, rungs=np.asarray(("A", "B")), gammas=np.asarray(GAMMAS), ns=np.asarray(NS), reps=np.asarray(REPS))


def sign_start(index, ux, uy):
    signs = np.asarray([1.0 if ((int(index) >> j) & 1) == 0 else -1.0 for j in range(D)], dtype=np.float64)
    return uy @ np.diag(signs) @ ux.T


def emit_finished_case(row):
    print(json.dumps({"event": "case_finished", "rung": row["rung"], "gamma": float(row["gamma"]), "n": int(row["n"]), "r": int(row["r"]), "err": float(row["err"]), "err_shape": float(row["err_shape"]), "L_hat": float(row["L_hat"]), "L_star": float(row["L_star"])}, sort_keys=True), flush=True)


def estimate_one_case(rung, gamma, n, r, case, landscape):
    path = case_cache_path(rung, case)
    if case_cache_complete(path):
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
    eig_x, ux = np.linalg.eigh(np.cov(x, rowvar=False))
    eig_y, uy = np.linalg.eigh(np.cov(y, rowvar=False))
    eig_x, ux = eig_x[::-1], ux[:, ::-1]
    eig_y, uy = eig_y[::-1], uy[:, ::-1]
    losses = np.empty(1 << D, dtype=np.float64)
    for s_index in range(1 << D):
        q0 = torch.as_tensor(sign_start(s_index, ux, uy), dtype=torch.float64)
        losses[s_index] = float(problem.loss(q0, False)[0].detach())
    screen_order = np.asarray(sorted(range(1 << D), key=lambda index: (float(losses[index]), int(index))), dtype=np.int64)
    selected = screen_order[:8]
    best, _ = refine(problem, [(int(index), sign_start(int(index), ux, uy)) for index in selected])
    q_hat = np.asarray(best["Q"], dtype=np.float64)
    q_star_tensor = torch.as_tensor(q_star, dtype=torch.float64)
    l_star = float(problem.loss(q_star_tensor, False)[0].detach())
    delta = q_hat - q_star
    p_shape = np.diag(np.asarray([0.0] * 3 + [1.0] * 7, dtype=np.float64))
    p_pos = np.diag(np.asarray([1.0] * 3 + [0.0] * 7, dtype=np.float64))
    rank = int(np.flatnonzero(screen_order == int(best["start"]))[0]) + 1
    row = {"rung": rung, "gamma": float(gamma), "n": int(n), "r": int(r), "err": float(np.linalg.norm(delta, ord="fro") / np.sqrt(20.0)), "err_shape": float(np.linalg.norm(delta @ p_shape, ord="fro") / np.sqrt(14.0)), "err_pos": float(np.linalg.norm(delta @ p_pos, ord="fro") / np.sqrt(6.0)), "L_hat": float(best["L_hat"]), "L_star": l_star, "screen_rank_best": rank}
    row.update({f"cov_eig_X_{j}": float(eig_x[j]) for j in range(D)})
    row.update({f"cov_eig_Y_{j}": float(eig_y[j]) for j in range(D)})
    save_case_cache(path, row, q_hat, np.arange(1 << D, dtype=np.int64), losses)
    emit_finished_case(row)
    return True


def estimate_stage():
    make_meta()
    landscape = load_o1k_landscape()
    new_cases = 0
    for rung in ("A", "B"):
        for gamma, n, r, case in case_specs():
            new_cases += int(estimate_one_case(rung, gamma, n, r, case, landscape))
            write_cases_csv()
            write_screen_csv()
            write_q_hats()
    completed = write_cases_csv()
    write_screen_csv()
    write_q_hats()
    print(json.dumps({"completed_cases": int(len(completed)), "new_cases": int(new_cases)}, sort_keys=True), flush=True)


def input_hashes(gamma, n, r, landscape):
    x, y, x_labels, y_labels, q_star = sample_case(gamma, n, r, landscape)
    hasher = hashlib.sha256()
    for value in (x, y, x_labels, y_labels, q_star):
        hasher.update(np.ascontiguousarray(value).tobytes())
    return hasher.hexdigest()


def check_stage():
    """Verify O1k's read-only case metadata, landscape, and all 360 cached fits."""
    make_meta()
    landscape = load_o1k_landscape()
    source_cases_path = O1K / "cases.csv"
    source_fits_path = O1K / "fit_symmetry.csv"
    if not source_cases_path.is_file() or not source_fits_path.is_file():
        raise FileNotFoundError("O1k cases.csv and fit_symmetry.csv are required read-only inputs")
    source_cases = pd.read_csv(source_cases_path)
    source_fits = pd.read_csv(source_fits_path)
    case_records = []
    case_metadata_match = True
    for gamma, n, r, case in case_specs():
        digest = input_hashes(gamma, n, r, landscape)
        for rung in ("A", "B"):
            path = O1K_CACHE / "cases" / f"{rung}_{case}.npz"
            ok = path.is_file()
            if ok:
                with np.load(path, allow_pickle=False) as data:
                    ok = str(data["rung"].item()) == rung and float(data["gamma"].item()) == float(gamma) and int(data["n"].item()) == int(n) and int(data["r"].item()) == int(r)
            row = source_cases.loc[(source_cases["rung"] == rung) & np.isclose(source_cases["gamma"], gamma) & (source_cases["n"] == n) & (source_cases["r"] == r)]
            ok = bool(ok and len(row) == 1)
            case_metadata_match = bool(case_metadata_match and ok)
            case_records.append({"rung": rung, "case": case, "case_archive": str(path), "input_sha256": digest, "metadata_matches": ok})
    fit_records = []
    expected_ids = set(source_fits["fit_id"].astype(str))
    fit_set_match = len(source_fits) == 360 and source_fits["fit_id"].is_unique
    for fit_id in sorted(expected_ids):
        path = O1K_CACHE / "drifts" / f"{fit_id}.pt"
        exists = path.is_file()
        fit_set_match = bool(fit_set_match and exists)
        fit_records.append({"fit_id": fit_id, "path": str(path), "sha256": sha(path) if exists else None})
    record = {"o1k_landscape": {"path": str(O1K_ARTIFACTS / "landscape.npz"), "sha256": sha(O1K_ARTIFACTS / "landscape.npz")}, "source_cases": {"path": str(source_cases_path), "sha256": sha(source_cases_path), "rows": int(len(source_cases))}, "source_fit_symmetry": {"path": str(source_fits_path), "sha256": sha(source_fits_path), "rows": int(len(source_fits))}, "case_inputs": case_records, "fit_sha256": fit_records, "inputs_match": int(case_metadata_match and fit_set_match)}
    atomic_json(ARTIFACTS / "inputs.json", record)
    print(json.dumps({"inputs_match": record["inputs_match"], "fit_sha256_rows": len(fit_records)}, sort_keys=True), flush=True)


def basin_cache_complete(path):
    path = Path(path)
    if not path.is_file():
        return False
    try:
        with np.load(path, allow_pickle=False) as data:
            return set(("i", "t", "theta", "final_err", "returned")).issubset(data.files) and len(data["i"]) == 20
    except (OSError, ValueError, KeyError):
        return False


def find_theta(q_star, a, target):
    """Bisection to the issued 1e-6 tolerance on the prescribed distance."""
    def distance(theta):
        return float(np.linalg.norm(q_star @ (expm(float(theta) * a) - np.eye(D)), ord="fro") / np.sqrt(20.0))
    low, high = 0.0, 1.0
    while distance(high) < target:
        high *= 2.0
        if high > 128.0:
            raise RuntimeError(f"Could not bracket basin target {target}")
    while high - low > 1e-6:
        mid = (low + high) / 2.0
        if distance(mid) < target:
            low = mid
        else:
            high = mid
    return (low + high) / 2.0


def write_basin_csv():
    rows = []
    for rung in ("A", "B"):
        for gamma in GAMMAS:
            for r in REPS:
                path = basin_cache_path(rung, case_id(gamma, 3000, r))
                if not basin_cache_complete(path):
                    continue
                with np.load(path, allow_pickle=False) as data:
                    for i, t, theta, final_err, returned in zip(data["i"], data["t"], data["theta"], data["final_err"], data["returned"]):
                        rows.append({"rung": rung, "gamma": float(gamma), "r": int(r), "i": int(i), "t": float(t), "theta": float(theta), "final_err": float(final_err), "returned": int(returned)})
    frame = pd.DataFrame(rows, columns=("rung", "gamma", "r", "i", "t", "theta", "final_err", "returned"))
    if not frame.empty:
        frame = frame.sort_values(["rung", "gamma", "r", "i", "t"], kind="stable")
    atomic_dataframe(OUT / "basin.csv", frame)
    return frame


def basin_one_case(rung, gamma, r, case, landscape):
    path = basin_cache_path(rung, case)
    if basin_cache_complete(path):
        return False
    x, y, _, _, q_star = sample_case(gamma, 3000, r, landscape)
    true_x, true_y = exact_models(gamma, q_star, landscape)
    mx, my = (true_x, true_y) if rung == "A" else get_models(case)
    problem = Problem(mx, my, x, y, r, d=D)
    rows = []
    for i in range(1, 6):
        generator = np.random.default_rng(20264400 + i)
        g = generator.standard_normal((D, D))
        skew = g - g.T
        a = skew / np.linalg.norm(skew, ord="fro")
        for target in (0.05, 0.10, 0.20, 0.40):
            theta = find_theta(q_star, a, target)
            q0 = q_star @ expm(theta * a)
            answer = optimise_single(problem, q0)
            final_err = float(np.linalg.norm(np.asarray(answer["Q"], dtype=np.float64) - q_star, ord="fro") / np.sqrt(20.0))
            rows.append((i, target, theta, final_err, int(final_err <= 0.10)))
    atomic_npz(path, i=np.asarray([row[0] for row in rows], dtype=np.int64), t=np.asarray([row[1] for row in rows], dtype=np.float64), theta=np.asarray([row[2] for row in rows], dtype=np.float64), final_err=np.asarray([row[3] for row in rows], dtype=np.float64), returned=np.asarray([row[4] for row in rows], dtype=np.int64))
    print(json.dumps({"event": "basin_case_finished", "rung": rung, "gamma": gamma, "r": r}, sort_keys=True), flush=True)
    return True


def basin_stage():
    make_meta()
    landscape = load_o1k_landscape()
    new_cases = 0
    for rung in ("A", "B"):
        for gamma in GAMMAS:
            for r in REPS:
                case = case_id(gamma, 3000, r)
                new_cases += int(basin_one_case(rung, gamma, r, case, landscape))
                write_basin_csv()
    completed = write_basin_csv()
    print(json.dumps({"completed_basin_rows": int(len(completed)), "new_basin_cases": int(new_cases)}, sort_keys=True), flush=True)


def scalar_median(frame, column, *, rung=None, gamma=None, n=None):
    subset = frame
    if rung is not None:
        subset = subset.loc[subset["rung"] == rung]
    if gamma is not None:
        subset = subset.loc[np.isclose(subset["gamma"].to_numpy(float), float(gamma), atol=1e-12, rtol=0.0)]
    if n is not None:
        subset = subset.loc[subset["n"] == int(n)]
    if len(subset) != 3:
        raise RuntimeError(f"Expected three replicate rows for {column}, found {len(subset)}")
    return float(subset[column].median())


def basin_radius(frame, rung, gamma, r):
    subset = frame.loc[(frame["rung"] == rung) & np.isclose(frame["gamma"].to_numpy(float), gamma, atol=1e-12, rtol=0.0) & (frame["r"] == int(r))]
    if len(subset) != 20:
        raise RuntimeError(f"Expected 20 basin rows for {rung}, {gamma}, {r}; found {len(subset)}")
    choices = []
    for target in (0.05, 0.10, 0.20, 0.40):
        at_target = subset.loc[np.isclose(subset["t"].to_numpy(float), target, atol=1e-12, rtol=0.0)]
        if len(at_target) != 5:
            raise RuntimeError("Each basin target requires five directions")
        if int(at_target["returned"].sum()) >= 4:
            choices.append(target)
    return float(max(choices) if choices else 0.0)


def require_report_inputs():
    paths = (ARTIFACTS / "inputs.json", OUT / "cases.csv", OUT / "basin.csv", ARTIFACTS / "Q_hats.npz")
    missing = [str(path) for path in paths if not path.is_file()]
    if missing:
        raise RuntimeError("--stage report requires completed earlier stages: " + ", ".join(missing))
    inputs = json.loads((ARTIFACTS / "inputs.json").read_text())
    cases = pd.read_csv(OUT / "cases.csv")
    basin = pd.read_csv(OUT / "basin.csv")
    if len(cases) != 72:
        raise RuntimeError(f"--stage report requires 72 estimated cases, found {len(cases)}")
    if len(basin) != 480:
        raise RuntimeError(f"--stage report requires 480 basin rows, found {len(basin)}")
    return inputs, cases, basin


def gated_values(inputs, cases, basin):
    source_cases = pd.read_csv(O1K / "cases.csv")
    merged = cases.merge(source_cases[["rung", "gamma", "n", "r", "L_star"]], on=["rung", "gamma", "n", "r"], how="left", suffixes=("", "_o1k"), validate="one_to_one")
    if len(merged) != 72 or merged["L_star_o1k"].isna().any():
        raise RuntimeError("Could not pair every O1l case with O1k L_star")
    err_a = [scalar_median(cases, "err", rung="A", gamma=gamma, n=n) for gamma in (1.15, 1.33, 1.6) for n in NS]
    radii_a = [basin_radius(basin, "A", gamma, r) for gamma in (1.15, 1.33, 1.6) for r in REPS]
    radii_b = [basin_radius(basin, "B", gamma, r) for gamma in (1.15, 1.33, 1.6) for r in REPS]
    values = {
        "n_cases": int(len(cases)),
        "inputs_match": int(inputs["inputs_match"]),
        "l_star_repro_max": float(np.max(np.abs(merged["L_star"].to_numpy(float) - merged["L_star_o1k"].to_numpy(float)))),
        "negctrl_err_shape_A_s": scalar_median(cases, "err_shape", rung="A", gamma=1.0, n=3000),
        "errA_s_max": float(max(err_a)),
        "errB_s_g115_n3000": scalar_median(cases, "err", rung="B", gamma=1.15, n=3000),
        "errB_s_g133_n3000": scalar_median(cases, "err", rung="B", gamma=1.33, n=3000),
        "errB_s_g160_n3000": scalar_median(cases, "err", rung="B", gamma=1.6, n=3000),
        "errB_s_g115_n10000": scalar_median(cases, "err", rung="B", gamma=1.15, n=10000),
        "errB_s_g133_n10000": scalar_median(cases, "err", rung="B", gamma=1.33, n=10000),
        "errB_s_g160_n10000": scalar_median(cases, "err", rung="B", gamma=1.6, n=10000),
        "basin_radius_A": float(np.median(radii_a)),
        "basin_radius_B": float(np.median(radii_b)),
    }
    validity_failures = []
    if values["n_cases"] != 72:
        validity_failures.append("n_cases != 72")
    if values["inputs_match"] != 1:
        validity_failures.append("inputs_match = 0")
    if values["l_star_repro_max"] > 1e-10:
        validity_failures.append("l_star_repro_max > 1e-10")
    if values["negctrl_err_shape_A_s"] < 0.30:
        validity_failures.append("negctrl_err_shape_A_s < 0.30")
    if validity_failures:
        outcome = "REVISE"
    elif values["errA_s_max"] > 0.10:
        outcome = "search_unfixed"
    else:
        gamma_min = next((gamma for gamma, key in ((1.15, "errB_s_g115_n3000"), (1.33, "errB_s_g133_n3000"), (1.6, "errB_s_g160_n3000")) if values[key] <= 0.10), None)
        outcome = "pbmc_range_ok" if gamma_min in (1.15, 1.33) else "anisotropy_limited" if gamma_min == 1.6 else "fit_limited"
    values["o1l_outcome"] = outcome
    if tuple(values) != GATED_KEYS:
        raise AssertionError("O1l gated key order changed")
    return values, validity_failures


def ladder_figure(cases):
    source = pd.read_csv(O1K / "cases.csv")
    fig, axes = plt.subplots(1, len(NS), figsize=(4.8 * len(NS), 4.3), sharey=True)
    for axis, n in zip(np.atleast_1d(axes), NS):
        for rung, color in (("A", "C0"), ("B", "C1")):
            current = [scalar_median(cases, "err", rung=rung, gamma=gamma, n=n) for gamma in GAMMAS]
            previous = [scalar_median(source, "err", rung=rung, gamma=gamma, n=n) for gamma in GAMMAS]
            axis.plot(GAMMAS, current, marker="o", color=color, label=f"O1l {rung}" if n == NS[0] else None)
            axis.plot(GAMMAS, previous, marker="x", linestyle="--", color=color, alpha=0.65, label=f"O1k random-start {rung}" if n == NS[0] else None)
        axis.axhline(0.10, color="black", linestyle=":", linewidth=1.0, label="0.10" if n == NS[0] else None)
        axis.set_title(f"n = {n}")
        axis.set_xlabel("anisotropy $\\gamma$")
        axis.set_xticks(GAMMAS)
        axis.grid(axis="y", alpha=0.25)
    axes[0].set_ylabel("median rotation error")
    axes[0].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT / "fig_ladder.png", dpi=180)
    plt.close(fig)


def basin_figure(basin):
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.2), sharey=True)
    for axis, rung in zip(axes, ("A", "B")):
        subset = basin.loc[basin["rung"] == rung]
        for gamma, color in zip(GAMMAS, ("C0", "C1", "C2", "C3")):
            g = subset.loc[np.isclose(subset["gamma"].to_numpy(float), gamma, atol=1e-12, rtol=0.0)]
            fractions = [float(g.loc[np.isclose(g["t"].to_numpy(float), target, atol=1e-12, rtol=0.0), "returned"].mean()) for target in (0.05, 0.10, 0.20, 0.40)]
            axis.plot((0.05, 0.10, 0.20, 0.40), fractions, marker="o", color=color, label=f"$\\gamma$={gamma}")
        axis.set_title(f"rung {rung}")
        axis.set_xlabel("start distance")
        axis.set_xticks((0.05, 0.10, 0.20, 0.40))
        axis.grid(alpha=0.25)
        axis.legend(fontsize=8)
    axes[0].set_ylabel("fraction returned")
    axes[0].set_ylim(-0.03, 1.03)
    fig.tight_layout()
    fig.savefig(OUT / "fig_basin.png", dpi=180)
    plt.close(fig)


def umap_coordinates(rep):
    import anndata as ad
    import scanpy as sc
    adata = ad.AnnData(X=np.zeros((len(rep), 1), dtype=np.float32))
    adata.obsm["view"] = np.asarray(rep, dtype=np.float64)
    sc.pp.neighbors(adata, n_neighbors=15, use_rep="view", random_state=20261014)
    sc.tl.umap(adata, min_dist=0.5, random_state=20261014)
    return np.asarray(adata.obsm["X_umap"], dtype=np.float64)


def mixing(rep, labels):
    neighbors = NearestNeighbors(n_neighbors=31).fit(rep).kneighbors(rep, return_distance=False)[:, 1:]
    return float(np.mean(np.mean(labels[neighbors] != labels[:, None], axis=1)) / 0.5)


def umap_figure(landscape):
    gamma, n, r = 1.33, 3000, 0
    path = case_cache_path("B", case_id(gamma, n, r))
    if not case_cache_complete(path):
        raise RuntimeError("Missing rung-B gamma=1.33 n=3000 r=0 estimate for UMAP")
    _, q_hat, _, _ = load_case_cache(path)
    x, y, x_labels, y_labels, q_star = sample_case(gamma, n, r, landscape)
    datasets = np.asarray(["X"] * len(x) + ["Y"] * len(y), dtype=str)
    wells = np.concatenate([x_labels, y_labels]).astype(int)
    views = (("no rotation", np.eye(D)), ("E5s estimate", q_hat), ("known answer", q_star))
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
        by_well = axes[index, 1]
        for well in range(K):
            mask = wells == well
            by_well.scatter(coords[mask, 0], coords[mask, 1], s=4, alpha=0.7, color=palette(well), label=f"well {well}", linewidths=0)
        by_well.set_title(f"{name}: wells; mix={mix:.3f}, sil={sil:.3f}")
        by_well.legend(markerscale=2, fontsize=8, ncol=2)
        for axis in (by_copy, by_well):
            axis.set_xticks([])
            axis.set_yticks([])
    fig.tight_layout()
    fig.savefig(OUT / "fig_umap.png", dpi=180)
    plt.close(fig)


def write_lockfile():
    tmp = atomic_path(OUT / "requirements.lock")
    try:
        with tmp.open("w") as handle:
            subprocess.run([sys.executable, "-m", "pip", "freeze"], check=True, stdout=handle)
        os.replace(tmp, OUT / "requirements.lock")
    finally:
        tmp.unlink(missing_ok=True)


def git_commit():
    return subprocess.check_output(["git", "-C", str(ROOT), "rev-parse", "HEAD"], text=True).strip()


def checkpoint_stage(stage):
    """Commit and push tracked O1l products after each stage of --stage all."""
    if os.environ.get("O1L_AUTO_CHECKPOINT", "1") != "1":
        return
    subprocess.run(["git", "-C", str(ROOT), "add", "reports/O1l"], check=True)
    changed = subprocess.run(["git", "-C", str(ROOT), "diff", "--cached", "--quiet"]).returncode
    if changed == 0:
        return
    if changed != 1:
        raise RuntimeError(f"git diff --cached failed after O1l {stage} stage: {changed}")
    subprocess.run(["git", "-C", str(ROOT), "commit", "-m", f"task: O1l {stage} stage"], check=True)
    subprocess.run(["git", "-C", str(ROOT), "push", "origin", "task/O1l"], check=True)


def regeneration_commands():
    command = 'cd "$HOME/Documents/scopalign-manus-O1l" && .venv/bin/python reports/O1l/scripts/o1l_run.py'
    lock_command = 'cd "$HOME/Documents/scopalign-manus-O1l" && .venv/bin/python -m pip freeze > reports/O1l/requirements.lock'
    result = {
        "artifacts/inputs.json": command + " --stage check",
        "cases.csv": command + " --stage estimate",
        "artifacts/screen.csv": command + " --stage estimate",
        "artifacts/Q_hats.npz": command + " --stage estimate",
        "basin.csv": command + " --stage basin",
        "fig_ladder.png": command + " --stage report",
        "fig_basin.png": command + " --stage report",
        "fig_umap.png": command + " --stage report",
        "requirements.lock": lock_command,
        "report.json": command + " --stage report",
        "scripts/o1l_run.py": command + " --stage all",
    }
    for key in GATED_KEYS:
        if key == "inputs_match":
            result[key] = command + " --stage check"
        elif key in ("n_cases", "l_star_repro_max", "negctrl_err_shape_A_s", "errA_s_max", "errB_s_g115_n3000", "errB_s_g133_n3000", "errB_s_g160_n3000", "errB_s_g115_n10000", "errB_s_g133_n10000", "errB_s_g160_n10000"):
            result[key] = command + " --stage estimate"
        elif key in ("basin_radius_A", "basin_radius_B"):
            result[key] = command + " --stage basin"
        else:
            result[key] = command + " --stage report"
    return result


def checksum_entries(paths):
    return [{"path": str(path), "sha256": sha(OUT / path)} for path in paths]


def report_stage():
    meta = make_meta()
    inputs, cases, basin = require_report_inputs()
    gated, validity_failures = gated_values(inputs, cases, basin)
    landscape = load_o1k_landscape()
    ladder_figure(cases)
    basin_figure(basin)
    umap_figure(landscape)
    write_lockfile()
    deliverable_paths = ("requirements.lock", "cases.csv", "basin.csv", "fig_ladder.png", "fig_basin.png", "fig_umap.png", "scripts/o1l_run.py")
    artifact_paths = ("artifacts/inputs.json", "artifacts/screen.csv", "artifacts/Q_hats.npz")
    report = {
        "task_id": "O1l", "executor": "manus", "branch": "task/O1l", "commit": git_commit(), "seed": 20261014,
        "started": meta["started"], "finished": datetime.now(timezone.utc).isoformat(), "status": "complete", "gated": gated,
        "deliverables": checksum_entries(deliverable_paths), "artifacts": checksum_entries(artifact_paths), "regeneration": regeneration_commands(),
        "deviations": [{"what": "Used a three-worker fork context for refinement and report-local matplotlib/fontconfig cache directories.", "why": "The local execution policy denied PyTorch's default spawned shared-memory helper and global font-cache writes; the fork context retains O1k's loss and LBFGS settings, and cache placement does not affect analysis values."}],
        "deviations_requested": [], "observations": [], "known_failures": KNOWN_FAILURES,
        "lockfile_sha256": sha(OUT / "requirements.lock"),
    }
    if set(report["gated"]) != set(GATED_KEYS) or len(report["gated"]) != len(GATED_KEYS):
        raise AssertionError("O1l report must contain exactly fourteen gated keys")
    atomic_json(OUT / "report.json", report)
    print(json.dumps({"gated": gated, "validity_failures": validity_failures}, sort_keys=True), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage", required=True, choices=("check", "estimate", "basin", "report", "all"))
    stage = parser.parse_args().stage
    if stage in ("check", "all"):
        check_stage()
        if stage == "all":
            checkpoint_stage("check")
    if stage in ("estimate", "all"):
        estimate_stage()
        if stage == "all":
            checkpoint_stage("estimate")
    if stage in ("basin", "all"):
        basin_stage()
        if stage == "all":
            checkpoint_stage("basin")
    if stage in ("report", "all"):
        report_stage()
        if stage == "all":
            checkpoint_stage("report")


if __name__ == "__main__":
    main()
