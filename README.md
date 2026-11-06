# viva-iirabm

<!-- BEGIN:dashboard -->
<!-- END:dashboard -->

A [process-bigraph](https://github.com/vivarium-collective/process-bigraph)
**real bridge** to the **Innate Immune Response Agent-Based Model (IIRABM)** of
An & Cockrell, used to reproduce:

> Cockrell RC, An G. *Examining the controllability of sepsis using genetic
> algorithms on an agent-based model of systemic inflammation.*
> PLoS Comput Biol 14(2): e1005876 (2018).
> <https://doi.org/10.1371/journal.pcbi.1005876>

The IIRABM is a 2-D agent-based model of the endothelial–blood interface that
simulates the human innate inflammatory response to injury/infection. Depending
on five parameters (injury size, microbial invasiveness, microbial toxigenesis,
environmental toxicity, host resilience) the *in silico* patient heals or dies
of overwhelming inflammation (sepsis). The paper searches for multi-stage
cytokine interventions — using a genetic algorithm — that steer the system from
death back to health.

This package wraps the **genuine An & Cockrell C++ simulator** (vendored under
`vendor/`, GPLv2+; see `vendor/PROVENANCE.md`). Nothing here reimplements the
immune model — every trajectory is produced by the original code, ported to a
portable **serial** build (the deprecated MPI C++ bindings are replaced with a
serial single-run driver; the physics is byte-for-byte upstream).

## What it provides

- **`IIRABMProcess`** — runs one real IIRABM simulation (parameter set × seed ×
  optional intervention) and replays the per-timestep trajectory (oxygen
  deficit, 12 cytokines, 4 cell populations, alive/dead) as a process-bigraph
  `Process`.
- **Composite generators** (dashboard Composites tab): `iirabm_baseline`
  (untreated sepsis) and `iirabm_treated` (under a multi-stage cytokine
  intervention).
- **`viva_iirabm.ga`** — a reduced-scale reproduction of the paper's genetic
  algorithm (tournament selection + uniform crossover over the allowed
  multiplier set), driving the real model to discover effective interventions.
- **A showcase investigation** reproducing the paper's figures (see
  `investigations/`).

## Installation

```bash
# development (editable)
uv venv .venv && source .venv/bin/activate
uv pip install -e ".[dev]"
```

Building the IIRABM binary needs a C++11 compiler (`c++` / `clang++` / `g++`);
it is compiled automatically on first use and cached. No MPI required.

> Once installed, the process registers via `viva_iirabm.core.build_core()`.

## Quick start

```python
from viva_iirabm.runtime import run_sim, mortality_rate

# one untreated trajectory (paper's 82%-death parameter set)
r = run_sim(seed=1, n_steps=8000)
print("died:", r["died"], "peak oxygen deficit:", r["peak_oxygen_deficit"])

# population mortality across 100 seeds (paper reports ~82%)
print("mortality:", mortality_rate(range(1, 101), n_steps=12500))
```

## Reproduction status

| Paper artifact | What it shows | Status here |
|---|---|---|
| **Fig 6 / baseline** | 100 stochastic trajectories, ~82% death | ✅ reproduced (~77–82% across 100 seeds) |
| **Fig 4** | Responder vs non-responder dynamics (oxy, GCSF, neutrophils) | ✅ reproduced |
| **Fig 2 / Table 1** | GA-identified multi-stage intervention lowers mortality | ◐ GA reproduced at reduced scale + Fig 2 digitized for replay |
| **Fig 3** | Generalizability across 8,291 parameter sets | ◐ subsampled (full sweep is supercomputer-scale) |

The genetic-algorithm *training* in the paper used up to ~1,000,000 CPU-hours on
Cray supercomputers; this repo reproduces the **method** at laptop scale, and is
explicit about what is scaled down (see the investigation's verdict).

## License

MIT for the Python package; GPLv2+ for the vendored IIRABM C++ (and the compiled
binary). See `LICENSE` and `vendor/PROVENANCE.md`.
