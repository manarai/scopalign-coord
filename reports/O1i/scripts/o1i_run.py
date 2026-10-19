#!/usr/bin/env python
"""O1i: pre-registered E5r operator-subspace rotation."""
import argparse
import hashlib
import importlib.util
import json
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime, timezone
from pathlib import Path

sys.dont_write_bytecode = True
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import torch
from scipy.stats import ortho_group

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "reports" / "O1i"
O1G_SCRIPT = Path("/Users/terooatt/Documents/scopalign-manus-O1g/reports/O1g/scripts/o1g_run.py")
O1H_TOPK = ROOT / "reports" / "O1h" / "topk.csv"
D = 30
N = 2759
REPS = (0, 1, 2)
KS = (3, 5, 8)
KNOWN = [
    "tests/test_phase_b_pitchfork.py::test_pitchfork_saddle_localisation[0.1]",
    "tests/test_phase_b_pitchfork.py::test_pitchfork_saddle_localisation[0.3]",
    "tests/test_phase_b_pitchfork.py::test_pitchfork_saddle_localisation[0.5]",
    "tests/test_phase_b_pitchfork.py::test_pitchfork_half_split_reproducibility",
]
WORKER_PROBLEM = None


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def import_o1g():
    if not O1G_SCRIPT.is_file():
        raise RuntimeError(f"Missing authorized O1g script: {O1G_SCRIPT}")
    spec = importlib.util.spec_from_file_location("o1g_o1i_reference", O1G_SCRIPT)
    if spec is None or spec.loader is None:
        raise RuntimeError("Could not import authorized O1g script")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


O1G = import_o1g()
NWORKERS = min(3, O1G.start_workers())


def ensure_out():
    OUT.mkdir(parents=True, exist_ok=True)


def fit_path(case, side, seed):
    local = O1G.local_fit(case, side, seed)
    if local.is_file():
        return local
    old = O1G.old_fit(case, side, seed)
    if seed in O1G.FIRST5 and old.is_file():
        return old
    raise RuntimeError(f"Missing cached fit for {case}_{side}_{seed}")


def record_fits(rows):
    path = OUT / "fits_used.csv"
    old = pd.read_csv(path).to_dict("records") if path.is_file() else []
    all_rows = {(str(x["case"]), str(x["side"]), int(x["seed"])): x for x in old}
    all_rows.update({(str(x["case"]), str(x["side"]), int(x["seed"])): x for x in rows})
    pd.DataFrame([all_rows[key] for key in sorted(all_rows)]).to_csv(path, index=False)


def load_models(case):
    mx, my, rows = [], [], []
    for seed in O1G.SEEDS20:
        for side, bucket in (("X", mx), ("Y", my)):
            path = fit_path(case, side, seed)
            saved = O1G.load_saved(path)
            bucket.append(O1G.build_model(saved))
            rows.append({
                "case": case,
                "side": side,
                "seed": int(seed),
                "path": str(path),
                "sha256": sha(path),
                "median_asym": float(saved["symmetry"]["median_asym"]),
                "max_asym": float(saved["symmetry"]["max_asym"]),
            })
    record_fits(rows)
    return mx, my


def case_r(rep):
    z = np.load(O1G.LAT1).astype(np.float32)
    if z.shape != (5518, D):
        raise RuntimeError(f"Expected O1g latent shape (5518, {D}), found {z.shape}")
    perm = np.random.default_rng(20261030 + rep).permutation(len(z))
    qstar = np.asarray(ortho_group.rvs(D, random_state=20261031 + rep), dtype=np.float64)
    x = z[np.sort(perm[:N])]
    y = (z[np.sort(perm[N:2 * N])].astype(np.float64) @ qstar.T).astype(np.float32)
    return x, y, qstar


def variance_order(x):
    return np.argsort(-np.var(np.asarray(x, dtype=np.float64), axis=0), kind="stable")


def polar(a):
    u, _, vt = np.linalg.svd(np.asarray(a, dtype=np.float64), full_matrices=False)
    return u @ vt


def principal_sine(vy, qstar, vx):
    singular = np.linalg.svd(vy.T @ qstar @ vx, compute_uv=False)
    return float(np.sqrt(max(0.0, 1.0 - float(np.min(singular)) ** 2)))


class CaseBase:
    """Frozen field values for one copy pair and O1g's fixed H subsets."""

    def __init__(self, mx, my, x, y, rep):
        self.mx, self.my = mx, my
        self.x = torch.as_tensor(np.asarray(x, dtype=np.float64), dtype=torch.float64)
        self.y = torch.as_tensor(np.asarray(y, dtype=np.float64), dtype=torch.float64)
        rng = np.random.default_rng(20261040 + rep)
        self.xa = self.x[rng.choice(len(x), size=min(256, len(x)), replace=False)]
        self.yb = self.y[rng.choice(len(y), size=min(256, len(y)), replace=False)]
        self.ux, _ = O1G.ensemble(self.mx, self.x, False, False)
        self.uy, _ = O1G.ensemble(self.my, self.y, False, False)
        _, self.jx = O1G.ensemble(self.mx, self.xa, True, False)
        _, self.jy = O1G.ensemble(self.my, self.yb, True, False)

    @staticmethod
    def basis(u, j, k):
        un = np.asarray(u.cpu(), dtype=np.float64)
        jn = np.asarray(j.cpu(), dtype=np.float64)
        flow = (un.T @ un) / len(un)
        curv = np.mean(jn @ np.swapaxes(jn, 1, 2), axis=0)
        s = flow / np.trace(flow) + curv / np.trace(curv)
        eig, vec = np.linalg.eigh(s)
        order = np.argsort(eig)[::-1]
        eig, vec = eig[order], vec[:, order]
        gap = float(eig[k - 1] / eig[k])
        return vec[:, :k], vec[:, k:], gap, eig


class RestrictedProblem:
    """E5r loss for a fixed k using only the operator-derived bases."""

    def __init__(self, base, k):
        self.base, self.k = base, int(k)
        self.mx, self.my = base.mx, base.my
        self.x, self.y, self.xa, self.yb = base.x, base.y, base.xa, base.yb
        self.ux, self.uy, self.jx, self.jy = base.ux, base.uy, base.jx, base.jy
        vx, ux_perp, self.gap_x, self.eig_x = base.basis(self.ux, self.jx, k)
        vy, uy_perp, self.gap_y, self.eig_y = base.basis(self.uy, self.jy, k)
        self.vx = torch.as_tensor(vx, dtype=torch.float64)
        self.vy = torch.as_tensor(vy, dtype=torch.float64)
        self.ux_perp = torch.as_tensor(ux_perp, dtype=torch.float64)
        self.uy_perp = torch.as_tensor(uy_perp, dtype=torch.float64)
        self.w = torch.as_tensor(polar(uy_perp.T @ ux_perp), dtype=torch.float64)
        self.px = self.vx @ self.vx.T
        self.py = self.vy @ self.vy.T
        self.du = torch.sum((self.ux @ self.px) ** 2)
        self.dv = torch.sum((self.uy @ self.py) ** 2)
        self.pjx = self.project_matrix(self.jx, self.px)
        self.pjy = self.project_matrix(self.jy, self.py)
        self.dx = torch.sum(self.pjx ** 2)
        self.dy = torch.sum(self.pjy ** 2)
        if min(float(self.du), float(self.dv), float(self.dx), float(self.dy)) <= 0:
            raise RuntimeError("Restricted-loss normalizer is non-positive")

    @staticmethod
    def project_matrix(j, p):
        return torch.einsum("ij,njk,kl->nil", p, j, p)

    def q_from_r(self, r):
        return self.vy @ r @ self.vx.T + self.uy_perp @ self.w @ self.ux_perp.T

    def loss_q(self, q, track):
        uyq, _ = O1G.ensemble(self.my, self.x @ q.T, False, track)
        uxq, _ = O1G.ensemble(self.mx, self.y @ q, False, track)
        d1 = (uyq - self.ux @ q.T) @ self.py
        d2 = (uxq - self.uy @ q) @ self.px
        dl = (torch.sum(d1 ** 2) / self.du + torch.sum(d2 ** 2) / self.dv) / 2.0
        _, jyq = O1G.ensemble(self.my, self.xa @ q.T, True, track)
        _, jxq = O1G.ensemble(self.mx, self.yb @ q, True, track)
        qjx = torch.einsum("ij,njk,lk->nil", q, self.jx, q)
        qtjy = torch.einsum("ji,njk,kl->nil", q, self.jy, q)
        h1 = self.project_matrix(jyq - qjx, self.py)
        h2 = self.project_matrix(jxq - qtjy, self.px)
        hl = (torch.sum(h1 ** 2) / self.dx + torch.sum(h2 ** 2) / self.dy) / 2.0
        return dl + hl, dl, hl


def candidates(k):
    ans = []
    for number in range(2 ** k):
        signs = np.asarray([1.0 if (number >> bit) & 1 else -1.0 for bit in range(k)], dtype=np.float64)
        ans.append((number, f"sign_{number:0{k}b}", np.diag(signs)))
    for index in range(20):
        ans.append((2 ** k + index, f"random_{index:02d}", np.asarray(ortho_group.rvs(k, random_state=20265000 + index), dtype=np.float64)))
    return ans


def init_worker(problem):
    global WORKER_PROBLEM
    torch.set_num_threads(1)
    WORKER_PROBLEM = problem


def screen_worker(item):
    order, name, r0 = item
    r = torch.as_tensor(r0, dtype=torch.float64)
    value, dpart, hpart = WORKER_PROBLEM.loss_q(WORKER_PROBLEM.q_from_r(r), False)
    return {"candidate_order": int(order), "candidate": name, "L0": float(value), "D0": float(dpart), "H0": float(hpart)}


class RestrictedClosure:
    def __init__(self, problem, r0, a):
        self.problem, self.r0, self.a = problem, r0, a

    def __call__(self):
        self.a.grad = None
        r = self.r0 @ torch.matrix_exp(self.a - self.a.T)
        value, _, _ = self.problem.loss_q(self.problem.q_from_r(r), True)
        value.backward()
        return value


def refine_worker(item):
    order, name, r0 = item
    torch.set_num_threads(1)
    r0t = torch.as_tensor(r0, dtype=torch.float64)
    a = torch.zeros((WORKER_PROBLEM.k, WORKER_PROBLEM.k), dtype=torch.float64, requires_grad=True)
    opt = torch.optim.LBFGS([a], lr=1.0, max_iter=500, tolerance_grad=1e-10, tolerance_change=1e-12, line_search_fn="strong_wolfe")
    opt.step(RestrictedClosure(WORKER_PROBLEM, r0t, a))
    r = (r0t @ torch.matrix_exp(a - a.T)).detach()
    q = WORKER_PROBLEM.q_from_r(r).detach()
    value, dpart, hpart = WORKER_PROBLEM.loss_q(q, False)
    return {"candidate_order": int(order), "candidate": name, "L_hat": float(value), "D": float(dpart), "H": float(hpart), "Q": q.cpu().numpy()}


def screen_rows(case, k, rep, problem):
    path = OUT / "screen.csv"
    old = pd.read_csv(path).to_dict("records") if path.is_file() else []
    present = {int(x["candidate_order"]): x for x in old if str(x["case"]) == case and int(x["k"]) == k and int(x["replicate"]) == rep}
    todo = [item for item in candidates(k) if item[0] not in present]
    if todo:
        with ProcessPoolExecutor(max_workers=NWORKERS, initializer=init_worker, initargs=(problem,)) as ex:
            computed = list(ex.map(screen_worker, todo))
        for row in computed:
            row.update({"case": case, "k": int(k), "replicate": int(rep)})
        old.extend(computed)
        pd.DataFrame(old).sort_values(["case", "k", "replicate", "candidate_order"]).to_csv(path, index=False)
        present.update({int(x["candidate_order"]): x for x in computed})
    rows = [present[item[0]] for item in candidates(k)]
    if len(rows) != 2 ** k + 20:
        raise RuntimeError("Incomplete sign/random screen")
    return sorted(rows, key=lambda x: (float(x["L0"]), int(x["candidate_order"])))[:8]


def solve(base, k, case, rep):
    problem = RestrictedProblem(base, k)
    selected = screen_rows(case, k, rep, problem)
    table = {order: (name, r0) for order, name, r0 in candidates(k)}
    starts = [(int(row["candidate_order"]), table[int(row["candidate_order"])][0], table[int(row["candidate_order"])][1]) for row in selected]
    with ProcessPoolExecutor(max_workers=NWORKERS, initializer=init_worker, initargs=(problem,)) as ex:
        refined = list(ex.map(refine_worker, starts))
    best = min(refined, key=lambda x: (x["L_hat"], x["candidate_order"]))
    return problem, best, refined


def store_q(name, q):
    path = OUT / "Q_hats.npz"
    arrays = {}
    if path.is_file():
        with np.load(path, allow_pickle=False) as old:
            arrays.update({key: np.asarray(old[key]) for key in old.files})
    arrays[name] = np.asarray(q, dtype=np.float64)
    np.savez(path, **arrays)


def write_records(path, rows, keys):
    pd.DataFrame(rows).sort_values(keys).to_csv(path, index=False)


def test_r_stage():
    ensure_out()
    path = OUT / "testR.csv"
    rows = pd.read_csv(path).to_dict("records") if path.is_file() else []
    done = {(int(x["k"]), int(x["replicate"])) for x in rows}
    for rep in REPS:
        missing = [k for k in KS if (k, rep) not in done]
        if not missing:
            continue
        x, y, qstar = case_r(rep)
        case = f"R_n{N}_r{rep}"
        mx, my = load_models(case)
        base = CaseBase(mx, my, x, y, rep)
        order = variance_order(x)
        for k in missing:
            problem, best, refined = solve(base, k, case, rep)
            qhat = np.asarray(best["Q"], dtype=np.float64)
            lstar, _, _ = problem.loss_q(torch.as_tensor(qstar, dtype=torch.float64), False)
            record = {
                "k": int(k), "n": int(N), "replicate": int(rep),
                "err_k": float(np.linalg.norm((qhat - qstar)[:, order[:k]], ord="fro") / np.sqrt(2.0 * k)),
                "err_30": float(np.linalg.norm(qhat - qstar, ord="fro") / np.sqrt(60.0)),
                "sub_k": principal_sine(np.asarray(problem.vy), qstar, np.asarray(problem.vx)),
                "gap_x": float(problem.gap_x), "gap_y": float(problem.gap_y),
                "L_hat": float(best["L_hat"]), "L_star_k": float(lstar),
                "n_good_starts": int(sum(row["L_hat"] <= 1.01 * best["L_hat"] for row in refined)),
                "best_start": int(best["candidate_order"]),
            }
            store_q(f"R_k{k}_r{rep}", qhat)
            rows.append(record)
            write_records(path, rows, ["k", "replicate"])
            done.add((k, rep))
            print(json.dumps({"case": case, **record}, sort_keys=True), flush=True)


def load_f_operators():
    ops = {}
    for side in ("H1", "H2"):
        for seed in O1G.FIRST5:
            path = O1G.CACHE / "operators" / f"F_{side}_{seed}.npy"
            if not path.is_file():
                raise RuntimeError(f"Missing O1g read-only Test F operator cache: {path}")
            ops[(side, seed)] = np.load(path)
    return ops


def views_local(z1, z2, q_e5b, q_e5r, qload, h1, h2, ops):
    l1 = h1.obs["cell_type"].astype(str).to_numpy()
    l2 = h2.obs["cell_type"].astype(str).to_numpy()
    ds = np.asarray(["H1"] * len(z1) + ["H2"] * len(z2), dtype=str)
    types = np.concatenate([l1, l2]).astype(str)
    tested = sorted(t for t in set(l1) & set(l2) if np.sum(l1 == t) >= 50 and np.sum(l2 == t) >= 50)
    mask = np.isin(types, tested)
    m1 = {seed: ops[("H1", seed)] for seed in O1G.FIRST5}
    m2 = {seed: ops[("H2", seed)] for seed in O1G.FIRST5}
    reps = {
        "no_rotation": (np.vstack([z1, z2]), np.eye(D)),
        "E5b_estimate": (np.vstack([z1 @ q_e5b.T, z2]), q_e5b),
        "E5r_k5": (np.vstack([z1 @ q_e5r.T, z2]), q_e5r),
        "known_answer": (np.vstack([z1 @ qload.T, z2]), qload),
    }
    rows, coords = [], {}
    for name, (rep, q) in reps.items():
        mix = O1G.mixing(rep, ds)
        sil = float(O1G.silhouette_score(rep[mask], types[mask], metric="euclidean"))
        fr = O1G.frame_ratio(z1.astype(np.float64), z2.astype(np.float64), l1, l2, m1, m2, q)
        rows.append({"view": name, "mix": mix, "sil": sil, "frame_ratio": fr})
        coords[name] = O1G.umap(rep)
    pd.DataFrame(rows).to_csv(OUT / "views.csv", index=False)
    fig, axes = plt.subplots(4, 2, figsize=(10, 17))
    colors = plt.get_cmap("tab10")
    for index, row in enumerate(rows):
        xy = coords[row["view"]]
        ax = axes[index, 0]
        for label, color in (("H1", "C0"), ("H2", "C1")):
            use = ds == label
            ax.scatter(xy[use, 0], xy[use, 1], s=4, alpha=.7, color=color, label=label, linewidths=0)
        ax.set_title(f"{row['view']}: half; mix={row['mix']:.3f}, sil={row['sil']:.3f}, frame={row['frame_ratio']:.3f}")
        ax.legend(markerscale=2, fontsize=8); ax.set_xticks([]); ax.set_yticks([])
        ax = axes[index, 1]
        for color_index, label in enumerate(sorted(set(types))):
            use = types == label
            ax.scatter(xy[use, 0], xy[use, 1], s=4, alpha=.7, color=colors(color_index % 10), label=label, linewidths=0)
        ax.set_title(f"{row['view']}: cell type; mix={row['mix']:.3f}, sil={row['sil']:.3f}, frame={row['frame_ratio']:.3f}")
        ax.legend(markerscale=2, fontsize=7, ncol=2); ax.set_xticks([]); ax.set_yticks([])
    fig.tight_layout(); fig.savefig(OUT / "fig_umap.png", dpi=180); plt.close(fig)
    return {row["view"]: row for row in rows}


def test_f_stage():
    if (OUT / "testF.csv").is_file() and (OUT / "Q_hat_F.npy").is_file() and (OUT / "views.csv").is_file() and (OUT / "fig_umap.png").is_file():
        return
    ensure_out()
    x, y, qload, check, h1, h2 = O1G.pca_halves()
    mx, my = load_models("F")
    base = CaseBase(mx, my, x, y, 0)
    problem, best, _ = solve(base, 5, "F", 0)
    qhat = np.asarray(best["Q"], dtype=np.float64)
    np.save(OUT / "Q_hat_F.npy", qhat)
    store_q("F_k5", qhat)
    q_e5b_path = O1G.OUT / "Q_hat_F.npy"
    if not q_e5b_path.is_file():
        raise RuntimeError(f"Missing O1g E5b Test F estimate: {q_e5b_path}")
    q_e5b = np.load(q_e5b_path).astype(np.float64)
    views = views_local(x, y, q_e5b, qhat, qload, h1, h2, load_f_operators())
    e = views["E5r_k5"]
    n = views["no_rotation"]
    b = views["E5b_estimate"]
    rec = {
        "mix_F": float(e["mix"]), "mix_F_noA": float(n["mix"]),
        "sil_F": float(e["sil"]), "sil_F_noA": float(n["sil"]),
        "frame_ratio_F": float(e["frame_ratio"]), "frame_ratio_F_noA": float(n["frame_ratio"]),
        "mix_E5b": float(b["mix"]), "sil_E5b": float(b["sil"]), "frame_ratio_E5b": float(b["frame_ratio"]),
        "err_F": float(np.linalg.norm(qhat - qload, ord="fro") / np.sqrt(60.0)),
        "best_start": int(best["candidate_order"]), "L_hat": float(best["L_hat"]),
        "latent_check_max": float(check), "gap_x": float(problem.gap_x), "gap_y": float(problem.gap_y),
    }
    pd.DataFrame([rec]).to_csv(OUT / "testF.csv", index=False)
    print(json.dumps({"case": "F_k5", **rec}, sort_keys=True), flush=True)


def symmetry_max():
    path = OUT / "fits_used.csv"
    if not path.is_file():
        return float("nan")
    return float(pd.read_csv(path)["median_asym"].max())


def fig_topk():
    if not O1H_TOPK.is_file():
        raise RuntimeError(f"Missing required O1h comparison: {O1H_TOPK}")
    tr = pd.read_csv(OUT / "testR.csv")
    current = tr.groupby("k")["err_k"].median()
    old = pd.read_csv(O1H_TOPK).set_index("k")
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.plot(list(KS), [float(current.loc[k]) for k in KS], marker="o", linewidth=2, label="E5r")
    ax.plot(list(KS), [float(old.loc[k, "err_hat"]) for k in KS], marker="s", linewidth=1.5, linestyle="--", label="O1h E5b")
    ax.axhline(.10, color="black", linestyle=":", linewidth=1, label="0.10")
    ax.set(xlabel="operator subspace dimension (k)", ylabel="median Test R top-k error", xticks=list(KS), ylim=(0, None))
    ax.legend(fontsize=8); fig.tight_layout(); fig.savefig(OUT / "fig_topk.png", dpi=180); plt.close(fig)


def entry(path):
    return {"path": path, "sha256": sha(OUT / path)}


def command():
    return 'cd "$HOME/Documents/scopalign-manus-O1i" && .venv/bin/python reports/O1i/scripts/o1i_run.py'


def report_stage(force_stopped=False):
    ensure_out()
    with (OUT / "requirements.lock").open("w") as f:
        subprocess.run([sys.executable, "-m", "pip", "freeze"], stdout=f, check=True)
    tr = pd.read_csv(OUT / "testR.csv")
    sym = symmetry_max()
    stop = force_stopped or sym > 1e-4
    gated = {
        "err_R_top5_2759": float(tr.loc[tr.k == 5, "err_k"].median()),
        "err_R_top3_2759": float(tr.loc[tr.k == 3, "err_k"].median()),
        "err_R_top8_2759": float(tr.loc[tr.k == 8, "err_k"].median()),
        "sub5_max": float(tr.loc[tr.k == 5, "sub_k"].max()),
        "frac_robust": float(np.mean(tr.loc[tr.k == 5, "n_good_starts"].to_numpy() >= 3)),
        "symmetry_max": float(sym),
    }
    if not stop:
        testf = pd.read_csv(OUT / "testF.csv").iloc[0]
        gated.update({
            "mix_F": float(testf.mix_F), "mix_F_noA": float(testf.mix_F_noA),
            "sil_F": float(testf.sil_F), "sil_F_noA": float(testf.sil_F_noA),
            "frame_ratio_F": float(testf.frame_ratio_F),
        })
        test_r_pass = gated["err_R_top5_2759"] <= .10 and gated["frac_robust"] >= 2 / 3
        test_f_pass = gated["mix_F"] >= .80 and gated["frame_ratio_F"] <= 1.5 and gated["sil_F"] >= .9 * gated["sil_F_noA"]
        gated["e5r_outcome"] = "works" if test_r_pass and test_f_pass else "fails"
        fig_topk()
        files = ["requirements.lock", "testR.csv", "testF.csv", "Q_hat_F.npy", "views.csv", "fig_umap.png", "fig_topk.png", "scripts/o1i_run.py"]
    else:
        files = ["requirements.lock", "testR.csv", "scripts/o1i_run.py"]
    base = command()
    regen = {path: base + " --stage all" for path in files}
    regen.update({
        "err_R_top5_2759": base + " --stage testR", "err_R_top3_2759": base + " --stage testR",
        "err_R_top8_2759": base + " --stage testR", "sub5_max": base + " --stage testR",
        "frac_robust": base + " --stage testR", "symmetry_max": base + " --stage testR",
    })
    if not stop:
        regen.update({"mix_F": base + " --stage testF", "mix_F_noA": base + " --stage testF", "sil_F": base + " --stage testF", "sil_F_noA": base + " --stage testF", "frame_ratio_F": base + " --stage testF", "e5r_outcome": base + " --stage report"})
    started = datetime.fromtimestamp((OUT / "testR.csv").stat().st_mtime, timezone.utc).isoformat()
    report = {
        "task_id": "O1i", "executor": "manus", "branch": "task/O1i",
        "commit": subprocess.check_output(["git", "-C", str(ROOT), "rev-parse", "HEAD"], text=True).strip(),
        "seed": 20261014, "started": started, "finished": datetime.now(timezone.utc).isoformat(),
        "status": "stopped" if stop else "complete", "gated": gated,
        "deliverables": [entry(path) for path in files],
        "artifacts": [entry(path) for path in ("Q_hats.npz", "screen.csv", "fits_used.csv")],
        "regeneration": regen, "deviations": [], "deviations_requested": [],
        "observations": [
            f"Authorized O1g import source SHA-256: {sha(O1G_SCRIPT)}.",
            f"O1h top-k comparison source SHA-256: {sha(O1H_TOPK)}.",
            "fits_used.csv records the absolute path and SHA-256 of every cached drift fit loaded by O1i.",
        ],
        "known_failures": KNOWN, "lockfile_sha256": sha(OUT / "requirements.lock"),
    }
    if stop:
        report["stop_rule"] = "Validity: symmetry_max > 1e-4 → REVISE."
    (OUT / "report.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps(gated, sort_keys=True), flush=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--stage", required=True, choices=("testR", "testF", "report", "all"))
    stage = parser.parse_args().stage
    if stage in ("testR", "all"):
        test_r_stage()
        if stage == "all" and symmetry_max() > 1e-4:
            report_stage(True)
            return
    if stage in ("testF", "all"):
        test_f_stage()
    if stage in ("report", "all"):
        report_stage()


if __name__ == "__main__":
    main()
