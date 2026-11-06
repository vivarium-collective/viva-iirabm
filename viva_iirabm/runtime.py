"""Build and drive the vendored IIRABM serial binary.

The upstream simulator is a C++ agent-based model. ``ensure_binary()`` applies
the reproducible port (``vendor/patch_iirabm.py``) to the pristine upstream
sources and compiles a serial executable, caching the result. ``run_sim()``
invokes that executable for one (parameter-set x seed x intervention) and parses
the per-timestep trajectory it writes.

This is a *real bridge*: every number returned here is produced by the genuine
An & Cockrell IIRABM C++ code, not a reimplementation.
"""
from __future__ import annotations

import hashlib
import os
import subprocess
import sys
import tempfile
import threading
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

# ---- layout -----------------------------------------------------------------
_PKG = Path(__file__).resolve().parent
_REPO = _PKG.parent
_VENDOR = _REPO / "vendor"
_UPSTREAM_V1 = _VENDOR / "iirabm_upstream" / "V1"
_PATCH = _VENDOR / "patch_iirabm.py"

# Series order written by the serial driver (see patch_iirabm.py / sepsisOutput).
SERIES = [
    "pmn", "mono", "TH1", "TH2", "oxygen_deficit", "total_infection",
    "TNF", "sTNFr", "IL10", "GCSF", "pTH1", "pTH2",
    "PAF", "IFNg", "IL4", "IL8", "IL12", "sIL1r", "IL1ra", "IL1",
]
CYTOKINES = ["PAF", "TNF", "sTNFr", "IL1", "sIL1r", "IL1ra",
             "IFNg", "IL4", "IL8", "IL10", "IL12", "GCSF"]
CELL_POPS = ["pmn", "mono", "TH1", "TH2"]

# Aggregate endothelial-cell-damage death threshold used by the model
# (oxyDeficit > 8161 ~= 80% of the 101x101 grid -> system failure / death).
DEATH_THRESHOLD = 8161.0

_build_lock = threading.Lock()


def _build_dir() -> Path:
    override = os.environ.get("VIVA_IIRABM_BUILD_DIR")
    if override:
        d = Path(override)
    elif os.access(_VENDOR, os.W_OK):
        d = _VENDOR / "_build"
    else:  # installed wheel: use a user cache dir
        d = Path(os.environ.get("XDG_CACHE_HOME", Path.home() / ".cache")) / "viva_iirabm"
    d.mkdir(parents=True, exist_ok=True)
    return d


def _sources_hash() -> str:
    h = hashlib.sha256()
    for f in sorted(_UPSTREAM_V1.glob("*")) + [_PATCH]:
        h.update(f.read_bytes())
    return h.hexdigest()[:12]


def ensure_binary(force: bool = False) -> Path:
    """Return the path to the compiled serial IIRABM binary, building if needed."""
    with _build_lock:
        bdir = _build_dir()
        tag = _sources_hash()
        binary = bdir / f"iirabm-{tag}"
        if binary.exists() and not force:
            return binary
        patched = bdir / f"src-{tag}"
        subprocess.run([sys.executable, str(_PATCH), str(_UPSTREAM_V1), str(patched)],
                       check=True, capture_output=True, text=True)
        cxx = os.environ.get("CXX", "c++")
        cpps = sorted(str(p) for p in patched.glob("*.cpp"))
        proc = subprocess.run([cxx, "-O2", "-o", str(binary), *cpps, "-std=c++11"],
                              capture_output=True, text=True)
        if proc.returncode != 0 or not binary.exists():
            raise RuntimeError(f"IIRABM build failed:\n{proc.stderr[-4000:]}")
        return binary


# ---- intervention encoding --------------------------------------------------
# Order of the 12 targets within each intervention stage (matches applyIntervention
# in the vendored simulationFunctions.cpp).
INTERVENTION_TARGETS = ["PAF", "TNF", "sTNFr", "IL1", "sIL1r", "IL1ra",
                        "IFNg", "IL4", "IL8", "IL10", "IL12", "GCSF"]


def write_intervention_file(path: Path, stages: list[list[float]]) -> None:
    """Write an intervention genome (list of 12-value stages) to a .dat file."""
    flat = [str(float(v)) for stage in stages for v in stage]
    path.write_text(" ".join(flat) + "\n")


# ---- running ----------------------------------------------------------------
def run_sim(
    *,
    injury: int = 33,
    host_resilience: float = 0.1,
    toxigenesis: int = 5,
    env_toxicity: int = 2,
    invasiveness: int = 2,
    seed: int = 1,
    n_steps: int = 10000,
    intervention: list[list[float]] | None = None,
    intervention_start: int = 505,
    intervention_separation: int = 100,
    intervention_length: int = 99,
    workdir: str | os.PathLike | None = None,
) -> dict:
    """Run one IIRABM simulation; return trajectories + outcome.

    Returns a dict with:
      - every key in SERIES -> np.ndarray[n_steps]
      - 'cytokines' -> {name: array}, 'cell_populations' -> {name: array}
      - 'died' (bool), 'death_step' (int|None), 'healed' (bool)
      - 'final_oxygen_deficit', 'peak_oxygen_deficit', 'final_infection'
      - the input parameter set echoed back.
    """
    binary = ensure_binary()
    tmp = Path(workdir) if workdir else Path(tempfile.mkdtemp(prefix="iirabm_"))
    tmp.mkdir(parents=True, exist_ok=True)
    out = tmp / f"run_s{seed}.dat"

    do_int = 0
    n_stages = 0
    int_path = ""
    if intervention:
        do_int = 1
        n_stages = len(intervention)
        ipath = tmp / f"int_s{seed}.dat"
        write_intervention_file(ipath, intervention)
        int_path = str(ipath)

    args = [str(binary), str(injury), str(host_resilience), str(toxigenesis),
            str(env_toxicity), str(invasiveness), str(seed), str(n_steps),
            str(do_int), str(n_stages), str(intervention_start),
            str(intervention_separation), str(intervention_length),
            int_path, str(out)]
    subprocess.run(args, check=True, capture_output=True, text=True)

    lines = out.read_text().splitlines()
    result: dict = {}
    for i, name in enumerate(SERIES):
        result[name] = np.array([float(x) for x in lines[1 + i].split()])

    oxy = result["oxygen_deficit"]
    # Death is encoded as -1 fill once oxyDeficit crosses the threshold.
    dead_mask = oxy < 0
    died = bool(dead_mask.any()) or bool(oxy[-1] > DEATH_THRESHOLD)
    death_step = int(np.argmax(dead_mask)) if dead_mask.any() else None
    final_infection = float(result["total_infection"][-1])
    # Healed = survived to end, infection resolved, damage recovering.
    healed = (not died) and final_infection == 0.0

    result.update(
        cytokines={k: result[k] for k in CYTOKINES},
        cell_populations={k: result[k] for k in CELL_POPS},
        died=died, death_step=death_step, healed=healed,
        final_oxygen_deficit=float(oxy[oxy >= 0][-1]) if (oxy >= 0).any() else float("nan"),
        peak_oxygen_deficit=float(oxy[oxy >= 0].max()) if (oxy >= 0).any() else float("nan"),
        final_infection=final_infection,
        params=dict(injury=injury, host_resilience=host_resilience,
                    toxigenesis=toxigenesis, env_toxicity=env_toxicity,
                    invasiveness=invasiveness, seed=seed, n_steps=n_steps),
    )
    return result


def run_batch(seeds, *, max_workers: int | None = None, **kwargs) -> list[dict]:
    """Run ``run_sim`` across many seeds in parallel (one process each)."""
    max_workers = max_workers or os.cpu_count()
    with tempfile.TemporaryDirectory(prefix="iirabm_batch_") as td:
        with ThreadPoolExecutor(max_workers=max_workers) as ex:
            futs = [ex.submit(run_sim, seed=s, workdir=td, **kwargs) for s in seeds]
            return [f.result() for f in futs]


def mortality_rate(seeds, **kwargs) -> float:
    """Fraction of seeds that die (the paper's population-level mortality)."""
    runs = run_batch(seeds, **kwargs)
    return sum(r["died"] for r in runs) / len(runs)
