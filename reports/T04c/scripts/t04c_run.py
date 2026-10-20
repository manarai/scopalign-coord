#!/usr/bin/env python3
"""T04c: within-B operator contrasts using read-only cached T04 fits."""
import argparse
import hashlib
import importlib.util
import json
import os
import platform
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import torch
from scjdo2.compare import centroid_operator
from scjdo2.energy_drift import EnergyDriftAdapter, EnergyDriftConfig, _EnergyMLP
import threadpoolctl

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "reports" / "T04c"
CACHE = OUT / "cache"
T04_ROOT = ROOT / "reports" / "T04"
T04_SCRIPT = T04_ROOT / "scripts" / "t04_run.py"
T04_SYMMETRY = T04_ROOT / "fit_symmetry.csv"
SEEDS = (20261014, 20261016, 20261017, 20261018, 20261019)


def sha(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def utcnow() -> str:
    return datetime.now(timezone.utc).isoformat()


def import_t04():
    """Import the T04 runner by absolute path without creating T04 bytecode."""
    if not T04_SCRIPT.exists():
        raise RuntimeError(f"T04 script missing: {T04_SCRIPT}")
    sys.dont_write_bytecode = True
    spec = importlib.util.spec_from_file_location("t04_cache_readonly", T04_SCRIPT)
    if spec is None or spec.loader is None:
        raise RuntimeError("Could not load cached T04 script")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


T04 = import_t04()


def make_meta() -> dict:
    OUT.mkdir(parents=True, exist_ok=True)
    CACHE.mkdir(parents=True, exist_ok=True)
    (OUT / ".gitignore").write_text("cache/\n__pycache__/\ngates.json\nrun.log\nrun.pid\n")
    path = OUT / "run_meta.json"
    if not path.exists():
        path.write_text(json.dumps({
            "started": utcnow(),
            "platform": platform.platform(),
            "blas": threadpoolctl.threadpool_info(),
            "t04_script_sha256": sha(T04_SCRIPT),
            "t04_fit_symmetry_sha256": sha(T04_SYMMETRY),
        }, indent=2, sort_keys=True) + "\n")
    return json.loads(path.read_text())


def gates() -> dict:
    path = OUT / "gates.json"
    return json.loads(path.read_text()) if path.exists() else {}


def update_gates(**values) -> dict:
    value = gates()
    value.update(values)
    (OUT / "gates.json").write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    return value


def write_requirements() -> None:
    with (OUT / "requirements.lock").open("w") as handle:
        subprocess.run([sys.executable, "-m", "pip", "freeze"], check=True, stdout=handle)


def cache_manifest() -> dict:
    paths = [T04_SCRIPT, T04_SYMMETRY]
    for scenario in ("S0", "S1"):
        paths.append(T04.emb_cache(scenario))
        for seed in SEEDS:
            paths.extend((T04.drift_path(scenario, "B", seed), T04.operator_path(scenario, "B", seed)))
    records = []
    for path in paths:
        if not path.exists():
            records.append({"path": str(path.relative_to(ROOT)), "exists": False})
        else:
            records.append({"path": str(path.relative_to(ROOT)), "exists": True, "sha256": sha(path)})
    payload = {"inputs": records}
    (OUT / "input_manifest.json").write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    return payload


def restore_drift(saved: dict):
    config = EnergyDriftConfig(**saved["config"])
    drift = EnergyDriftAdapter(T04.D, config=config)
    model = _EnergyMLP(T04.D, int(saved["config"]["hidden"]), int(saved["config"]["depth"]))
    model.load_state_dict(saved["state_dict"])
    drift._model = model
    drift._device = "cpu"
    drift._fitted = True
    return drift


def check_stage() -> None:
    make_meta()
    cache_manifest()
    # T04's fit_worker set this before calling symmetry_check; retain that
    # setting so exact float32 replay is evaluated under the saved protocol.
    torch.set_num_threads(1)
    # The committed values were serialized from float32 autodiff results.
    # Pandas' default fast parser may round some decimals by one float64 ULP;
    # round_trip preserves the committed decimal values for exact comparison.
    expected = pd.read_csv(T04_SYMMETRY, float_precision="round_trip")
    expected = expected[expected["side"] == "B"].sort_values(["scenario", "seed"]).reset_index(drop=True)
    rows, matched = [], True
    for scenario in ("S0", "S1"):
        embedding = T04.load_embedding(scenario)
        z_b = np.asarray(embedding["z_b"], dtype=np.float32)
        for seed in SEEDS:
            saved = torch.load(T04.drift_path(scenario, "B", seed), map_location="cpu", weights_only=False)
            reproduced = restore_drift(saved).symmetry_check(z_b, n_cells=100)
            hit = expected[(expected["scenario"] == scenario) & (expected["seed"] == seed)]
            if len(hit) != 1:
                raise RuntimeError(f"Expected one T04 B symmetry row for {scenario}/{seed}, got {len(hit)}")
            target = hit.iloc[0]
            median_match = float(reproduced["median_asym"]) == float(target["median_asym"])
            max_match = float(reproduced["max_asym"]) == float(target["max_asym"])
            matched = matched and median_match and max_match
            rows.append({
                "scenario": scenario,
                "side": "B",
                "seed": int(seed),
                "expected_median_asym": float(target["median_asym"]),
                "recomputed_median_asym": float(reproduced["median_asym"]),
                "expected_max_asym": float(target["max_asym"]),
                "recomputed_max_asym": float(reproduced["max_asym"]),
                "median_exact_match": bool(median_match),
                "max_exact_match": bool(max_match),
            })
    if len(expected) != 10:
        matched = False
    table = pd.DataFrame(rows).sort_values(["scenario", "seed"])
    table.to_csv(OUT / "input_symmetry.csv", index=False)
    symmetry_max = float(table["recomputed_median_asym"].max())
    result = {
        "inputs_match": int(bool(matched)),
        "symmetry_max": symmetry_max,
        "n_b_fits": int(len(table)),
    }
    (OUT / "check.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    update_gates(**result)
    print(json.dumps(result, sort_keys=True), flush=True)


def require_inputs_match() -> None:
    value = gates()
    if "inputs_match" not in value:
        raise RuntimeError("Run --stage check first")
    if int(value["inputs_match"]) != 1:
        raise RuntimeError("inputs_mismatch: contrasts must not be run")


def cached_operator(scenario: str, seed: int, z_b: np.ndarray) -> np.ndarray:
    path = T04.operator_path(scenario, "B", seed)
    if path.exists():
        return np.asarray(np.load(path), dtype=np.float64)
    # T04 cache is read-only for this task. If an operator cache is absent,
    # recompute T04's individual operator in memory rather than writing T04.
    saved = torch.load(T04.drift_path(scenario, "B", seed), map_location="cpu", weights_only=False)
    return np.asarray(T04.individual_operator(T04.build_model(saved), z_b), dtype=np.float64)


def centroid(z: np.ndarray, matrices: np.ndarray, indices: np.ndarray) -> np.ndarray:
    return centroid_operator(z[indices], matrices[indices], n_cells_centroid=30).astype(np.float64)


def ratio(numerator: float, denominator: float) -> float:
    return float(numerator / denominator) if denominator else float("inf")


def relabel_plan(indices: np.ndarray, n_left: int) -> dict:
    rng = np.random.default_rng(20261200)
    plan = {}
    for seed in SEEDS:
        draws = []
        for _ in range(20):
            shuffled = rng.permutation(indices)
            draws.append((shuffled[:n_left], shuffled[n_left:]))
        plan[int(seed)] = draws
    return plan


def operator_contrast(z: np.ndarray, matrices: dict, left: np.ndarray, right: np.ndarray, plan: dict):
    rows, all_nulls, null_by_seed = [], [], {}
    for seed in SEEDS:
        observed = T04.residual(centroid(z, matrices[seed], left), centroid(z, matrices[seed], right))
        values = [T04.residual(centroid(z, matrices[seed], draw_left), centroid(z, matrices[seed], draw_right))
                  for draw_left, draw_right in plan[int(seed)]]
        null_by_seed[int(seed)] = values
        all_nulls.extend(values)
        null_median = float(np.median(values))
        rows.append({"seed": int(seed), "observed": float(observed), "null_median": null_median,
                     "ratio": ratio(float(observed), null_median)})
    contrast = float(np.median([row["observed"] for row in rows]))
    null = float(np.median(all_nulls))
    return {"rows": rows, "contrast": contrast, "null": null, "ratio": ratio(contrast, null),
            "per_seed_ratio": np.asarray([row["ratio"] for row in rows], dtype=np.float64),
            "null_by_seed": null_by_seed}


def position_contrast(z: np.ndarray, left: np.ndarray, right: np.ndarray, plan: dict):
    observed = float(np.linalg.norm(np.mean(z[left], axis=0) - np.mean(z[right], axis=0)))
    rows, all_nulls = [], []
    for seed in SEEDS:
        values = [float(np.linalg.norm(np.mean(z[draw_left], axis=0) - np.mean(z[draw_right], axis=0)))
                  for draw_left, draw_right in plan[int(seed)]]
        all_nulls.extend(values)
        null_median = float(np.median(values))
        rows.append({"seed": int(seed), "observed": observed, "null_median": null_median,
                     "ratio": ratio(observed, null_median)})
    return {"rows": rows, "observed": observed, "null": float(np.median(all_nulls)),
            "ratio": ratio(observed, float(np.median(all_nulls)))}


def contrasts_for_scenario(scenario: str, response_plan: dict):
    embedding = T04.load_embedding(scenario)
    z_b = np.asarray(embedding["z_b"], dtype=np.float64)
    n_a = int(len(embedding["z_a"]))
    types_b = np.asarray(embedding["types"].astype(str)[n_a:], dtype=str)
    responder_b = np.asarray(embedding["responder"][n_a:], dtype=bool)
    left = np.flatnonzero(responder_b)
    right = np.flatnonzero((types_b == "Monocytes") & ~responder_b)
    if len(left) == 0 or len(right) == 0:
        raise RuntimeError(f"Empty responder or nonresponder group for {scenario}")
    matrices = {int(seed): cached_operator(scenario, seed, z_b) for seed in SEEDS}
    contrast = operator_contrast(z_b, matrices, left, right, response_plan)
    position = position_contrast(z_b, left, right, response_plan)
    return {"z_b": z_b, "types_b": types_b, "responder_b": responder_b, "matrices": matrices,
            "left": left, "right": right, "contrast": contrast, "position": position}


def contrast_stage() -> None:
    require_inputs_match()
    base = T04.load_embedding("S1")
    n_a = int(len(base["z_a"]))
    types_b = np.asarray(base["types"].astype(str)[n_a:], dtype=str)
    responder_b = np.asarray(base["responder"][n_a:], dtype=bool)
    response_left = np.flatnonzero(responder_b)
    response_right = np.flatnonzero((types_b == "Monocytes") & ~responder_b)
    response_plan = relabel_plan(np.concatenate([response_left, response_right]), len(response_left))
    values = {scenario: contrasts_for_scenario(scenario, response_plan) for scenario in ("S0", "S1")}

    type_left = np.flatnonzero(values["S1"]["types_b"] == "Monocytes")
    type_right = np.flatnonzero(values["S1"]["types_b"] == "T cells")
    type_plan = relabel_plan(np.concatenate([type_left, type_right]), len(type_left))
    type_contrast = operator_contrast(values["S1"]["z_b"], values["S1"]["matrices"], type_left, type_right, type_plan)

    rows = []
    for scenario in ("S0", "S1"):
        for row in values[scenario]["contrast"]["rows"]:
            rows.append({"scenario": scenario, "comparison": "responders_vs_nonresponders", **row})
        for row in values[scenario]["position"]["rows"]:
            rows.append({"scenario": scenario, "comparison": "position_responders_vs_nonresponders", **row})
    for row in type_contrast["rows"]:
        rows.append({"scenario": "S1", "comparison": "monocytes_vs_t_cells", **row})
    pd.DataFrame(rows).sort_values(["scenario", "comparison", "seed"]).to_csv(OUT / "contrasts.csv", index=False)

    summary = {
        "contrast_ratio_S1": float(values["S1"]["contrast"]["ratio"]),
        "contrast_ratio_S0": float(values["S0"]["contrast"]["ratio"]),
        "seed_agreement_S1": float(np.mean(values["S1"]["contrast"]["per_seed_ratio"] >= 2.0)),
        "type_ratio_S1": float(type_contrast["ratio"]),
        "pos_ratio_S1": float(values["S1"]["position"]["ratio"]),
        "pos_ratio_S0": float(values["S0"]["position"]["ratio"]),
    }
    (OUT / "contrast_details.json").write_text(json.dumps({
        "summary": summary,
        "response_nulls": {scenario: {str(seed): [float(x) for x in values[scenario]["contrast"]["null_by_seed"][seed]] for seed in SEEDS}
                           for scenario in ("S0", "S1")},
    }, indent=2, sort_keys=True) + "\n")
    update_gates(**summary)
    print(json.dumps(summary, sort_keys=True), flush=True)


def contrast_figure() -> None:
    details = json.loads((OUT / "contrast_details.json").read_text())
    table = pd.read_csv(OUT / "contrasts.csv")
    summary = details["summary"]
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    x = np.arange(len(SEEDS))
    for scenario, color in (("S1", "C3"), ("S0", "C0")):
        rows = table[(table["scenario"] == scenario) & (table["comparison"] == "responders_vs_nonresponders")].sort_values("seed")
        normalized_null = []
        for seed, null_median in zip(rows["seed"], rows["null_median"]):
            draws = np.asarray(details["response_nulls"][scenario][str(int(seed))], dtype=np.float64)
            normalized_null.extend(draws / float(null_median))
        lower, upper = np.quantile(normalized_null, [0.025, 0.975])
        axes[0].fill_between(x, lower, upper, color=color, alpha=0.12)
        axes[0].plot(x, rows["ratio"].to_numpy(dtype=float), marker="o", color=color, label=scenario)
    axes[0].axhline(1.0, color="0.25", ls="--", lw=1, label="null median")
    axes[0].axhline(2.0, color="0.25", ls=":", lw=1, label="ratio 2")
    axes[0].set_xticks(x, [str(seed) for seed in SEEDS], rotation=25)
    axes[0].set_ylabel("operator contrast / relabelling null")
    axes[0].set_xlabel("drift-fit seed")
    axes[0].set_title("Responder contrast; shaded empirical null band")
    axes[0].legend()

    labels = ["type S1", "position S1", "position S0"]
    bars = [summary["type_ratio_S1"], summary["pos_ratio_S1"], summary["pos_ratio_S0"]]
    colors = ["C2", "C3", "C0"]
    axes[1].bar(labels, bars, color=colors)
    axes[1].axhline(1.0, color="0.25", ls="--", lw=1)
    axes[1].axhline(2.0, color="0.25", ls=":", lw=1)
    axes[1].set_ylabel("contrast / relabelling null")
    axes[1].set_title("Sanity and position baselines")
    for index, value in enumerate(bars):
        axes[1].text(index, value, f"{value:.2f}", ha="center", va="bottom", fontsize=9)
    fig.tight_layout()
    fig.savefig(OUT / "fig_contrast.png", dpi=180)
    plt.close(fig)


def report_stage() -> None:
    make_meta()
    write_requirements()
    value = gates()
    if "inputs_match" not in value:
        raise RuntimeError("Run --stage check before --stage report")
    stopped = int(value["inputs_match"]) == 0
    if stopped:
        gated = {"inputs_match": int(value["inputs_match"]), "symmetry_max": float(value["symmetry_max"])}
        outcome = None
        stop_rule = "inputs_mismatch"
        files = ["requirements.lock", "scripts/t04c_run.py"]
    else:
        required = ("contrast_ratio_S1", "contrast_ratio_S0", "seed_agreement_S1", "type_ratio_S1", "pos_ratio_S1", "pos_ratio_S0")
        if any(key not in value for key in required):
            raise RuntimeError("Run --stage contrast before --stage report")
        if value["contrast_ratio_S1"] >= 2.0 and value["contrast_ratio_S0"] <= 1.5 and value["seed_agreement_S1"] >= 0.8:
            outcome = "detects"
        elif value["contrast_ratio_S1"] <= 1.5:
            outcome = "no_detect"
        else:
            outcome = "weak"
        contrast_figure()
        gated = {
            "inputs_match": int(value["inputs_match"]),
            "contrast_ratio_S1": float(value["contrast_ratio_S1"]),
            "contrast_ratio_S0": float(value["contrast_ratio_S0"]),
            "seed_agreement_S1": float(value["seed_agreement_S1"]),
            "type_ratio_S1": float(value["type_ratio_S1"]),
            "pos_ratio_S1": float(value["pos_ratio_S1"]),
            "pos_ratio_S0": float(value["pos_ratio_S0"]),
            "symmetry_max": float(value["symmetry_max"]),
            "t04c_outcome": outcome,
        }
        stop_rule = None
        files = ["requirements.lock", "contrasts.csv", "fig_contrast.png", "scripts/t04c_run.py"]
    command = 'cd "$HOME/scopalign-cloud" && .venv/bin/python reports/T04c/scripts/t04c_run.py'
    regeneration = {path: command + " --stage all" for path in files}
    for key in gated:
        stage = "check" if key in {"inputs_match", "symmetry_max"} else "contrast" if key != "t04c_outcome" else "report"
        regeneration[key] = command + " --stage " + stage
    artifacts = []
    for path in ("input_manifest.json", "input_symmetry.csv", "check.json", "contrast_details.json", "run_meta.json"):
        target = OUT / path
        if target.exists():
            artifacts.append({"path": path, "sha256": sha(target)})
    meta = make_meta()
    report = {
        "task_id": "T04c",
        "executor": "manus",
        "branch": "task/T04c",
        "commit": subprocess.check_output(["git", "-C", str(ROOT), "rev-parse", "HEAD"], text=True).strip(),
        "seed": 20261014,
        "started": meta["started"],
        "finished": utcnow(),
        "status": "stopped" if stopped else "complete",
        "gated": gated,
        "deliverables": [{"path": path, "sha256": sha(OUT / path)} for path in files],
        "artifacts": artifacts,
        "regeneration": regeneration,
        "deviations": [],
        "deviations_requested": [],
        "observations": [
            f"Read-only T04 fit-symmetry SHA-256: {sha(T04_SYMMETRY)}.",
            f"Read-only T04 runner SHA-256: {sha(T04_SCRIPT)}.",
        ],
        "known_failures": T04.KNOWN,
        "lockfile_sha256": sha(OUT / "requirements.lock"),
    }
    if stopped:
        report["stop_rule"] = stop_rule
    (OUT / "report.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    payload = dict(gated)
    if outcome is not None:
        payload["t04c_outcome"] = outcome
    print(json.dumps(payload, sort_keys=True), flush=True)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--stage", required=True, choices=("check", "contrast", "report", "all"))
    stage = parser.parse_args().stage
    stages = ("check", "contrast", "report") if stage == "all" else (stage,)
    for current in stages:
        if current == "check":
            check_stage()
        elif current == "contrast":
            if int(gates().get("inputs_match", 0)) == 0:
                continue
            contrast_stage()
        elif current == "report":
            report_stage()


if __name__ == "__main__":
    main()
