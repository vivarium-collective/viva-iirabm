"""Process-bigraph bridge for the IIRABM (Innate Immune Response ABM).

``IIRABMProcess`` is a **real bridge** to the genuine An & Cockrell C++
simulator (vendored under ``vendor/``). Because the IIRABM is a whole-trajectory
agent-based model — it advances every agent internally for the full run — the
Process runs the real binary once (lazily, on first ``update``), caches the
per-timestep trajectory, and then *replays* it one step per interval. Replaying
lets the bigraph engine, the dashboard's live view, and the study behavior-test
evaluator all read the real trajectory as an ordinary stream of observables.

Outputs are emitted as **deltas** against the previous step, so the bare
``float`` / ``map[string,float]`` ports accumulate to the correct absolute
value under the engine's additive apply rule (the v2ecoli mass/length pattern).
"""
from __future__ import annotations

from process_bigraph import Process

from .runtime import CELL_POPS, CYTOKINES, DEATH_THRESHOLD, run_sim


class IIRABMProcess(Process):
    """Replay one real IIRABM simulation as a process-bigraph Process.

    Config
    ------
    injury, host_resilience, toxigenesis, env_toxicity, invasiveness :
        The five IIRABM initialization parameters (injury radius in cells;
        resilience = per-step oxy heal rate; toxigenesis = infect spread;
        env_toxicity = recurring-injury count; invasiveness = infect repeat).
        Defaults are the paper's trained parameter set (82%-death case).
    seed : int
        RNG seed — selects the stochastic replicate.
    n_steps : int
        Simulation length in model timesteps (1 step ~= 7 min).
    intervention_flat : list[float]
        A flattened intervention genome (``intervention_stages`` blocks of 12
        per-cytokine multipliers, target order PAF, TNF, sTNFr, IL1, sIL1r,
        IL1ra, IFNg, IL4, IL8, IL10, IL12, GCSF). Empty = untreated.
    intervention_stages : int
        Number of sequential intervention stages encoded in ``intervention_flat``.
    intervention_start, intervention_separation, intervention_length : int
        Intervention schedule (paper defaults: start 505, 100 apart, length 99).

    Outputs
    -------
    oxygen_deficit, total_infection : float (delta)
        Systemic oxygen deficit (inverse health proxy) and total infection.
    cytokines, cell_populations : map[string,float] (per-key delta)
    alive : boolean
        False once aggregate damage crosses the death threshold.
    """

    config_schema = {
        "injury": {"_type": "integer", "_default": 33},
        "host_resilience": {"_type": "float", "_default": 0.1},
        "toxigenesis": {"_type": "integer", "_default": 5},
        "env_toxicity": {"_type": "integer", "_default": 2},
        "invasiveness": {"_type": "integer", "_default": 2},
        "seed": {"_type": "integer", "_default": 1},
        "n_steps": {"_type": "integer", "_default": 10000},
        "intervention_flat": {"_type": "list[float]", "_default": []},
        "intervention_stages": {"_type": "integer", "_default": 0},
        "intervention_start": {"_type": "integer", "_default": 505},
        "intervention_separation": {"_type": "integer", "_default": 100},
        "intervention_length": {"_type": "integer", "_default": 99},
    }

    def __init__(self, config=None, core=None):
        super().__init__(config, core)
        self._traj = None
        self._i = 0
        self._prev = None

    # No live inputs: the IIRABM computes an entire trajectory internally, so
    # nothing upstream can influence a step mid-replay. The scientific knobs are
    # construction-time config. (Documented rather than faked with a dummy port.)
    def inputs(self):
        return {}

    def outputs(self):
        return {
            "oxygen_deficit": "float",
            "total_infection": "float",
            "cytokines": "map[string,float]",
            "cell_populations": "map[string,float]",
            "alive": "boolean",
        }

    def initial_state(self):
        return {
            "oxygen_deficit": 0.0,
            "total_infection": 0.0,
            "cytokines": {k: 0.0 for k in CYTOKINES},
            "cell_populations": {k: 0.0 for k in CELL_POPS},
            "alive": True,
        }

    def _ensure_run(self):
        if self._traj is not None:
            return
        c = self.config
        stages = None
        if c["intervention_stages"] and c["intervention_flat"]:
            flat = list(c["intervention_flat"])
            n = c["intervention_stages"]
            stages = [flat[i * 12:(i + 1) * 12] for i in range(n)]
        self._traj = run_sim(
            injury=c["injury"], host_resilience=c["host_resilience"],
            toxigenesis=c["toxigenesis"], env_toxicity=c["env_toxicity"],
            invasiveness=c["invasiveness"], seed=c["seed"], n_steps=c["n_steps"],
            intervention=stages, intervention_start=c["intervention_start"],
            intervention_separation=c["intervention_separation"],
            intervention_length=c["intervention_length"],
        )
        self._prev = {
            "oxygen_deficit": 0.0, "total_infection": 0.0,
            "cytokines": {k: 0.0 for k in CYTOKINES},
            "cell_populations": {k: 0.0 for k in CELL_POPS},
        }

    def update(self, state, interval):
        self._ensure_run()
        step = min(self._i + max(1, int(interval)), len(self._traj["oxygen_deficit"]))
        j = step - 1
        self._i = step

        oxy = float(self._traj["oxygen_deficit"][j])
        # A dead patient's series is -1 from death on; clamp to the threshold so
        # the replayed absolute stays meaningful.
        oxy_abs = DEATH_THRESHOLD if oxy < 0 else oxy
        infe = float(self._traj["total_infection"][j])
        infe_abs = 0.0 if infe < 0 else infe

        def delta(scalar, key, sub=None):
            prev = self._prev[key] if sub is None else self._prev[key][sub]
            d = scalar - prev
            if sub is None:
                self._prev[key] = scalar
            else:
                self._prev[key][sub] = scalar
            return d

        upd = {
            "oxygen_deficit": delta(oxy_abs, "oxygen_deficit"),
            "total_infection": delta(infe_abs, "total_infection"),
            "cytokines": {k: delta(float(self._traj[k][j]), "cytokines", k)
                          for k in CYTOKINES},
            "cell_populations": {k: delta(float(self._traj[k][j]), "cell_populations", k)
                                 for k in CELL_POPS},
            "alive": not (oxy < 0 or oxy_abs >= DEATH_THRESHOLD),
        }
        return upd
