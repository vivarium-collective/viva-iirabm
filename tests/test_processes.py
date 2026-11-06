"""Tests for the IIRABM real bridge. Offline-safe (builds + runs the vendored
C++ locally; no network)."""
import numpy as np
import pytest

from process_bigraph import Composite, gather_emitter_results
from viva_iirabm.core import build_core
from viva_iirabm.processes import IIRABMProcess
from viva_iirabm.runtime import CYTOKINES, SERIES, ensure_binary, run_sim


def test_binary_builds():
    assert ensure_binary().exists()


def test_run_sim_shapes():
    r = run_sim(seed=1, n_steps=400)
    for name in SERIES:
        assert r[name].shape == (400,)
    assert set(r["cytokines"]) == set(CYTOKINES)
    assert isinstance(r["died"], bool)
    assert r["total_infection"][0] > 0  # injury seeds an infection


def test_run_sim_is_deterministic_per_seed():
    a = run_sim(seed=7, n_steps=300)["oxygen_deficit"]
    b = run_sim(seed=7, n_steps=300)["oxygen_deficit"]
    assert np.array_equal(a, b)


def test_process_update_ports():
    core = build_core()
    proc = IIRABMProcess(config={"seed": 1, "n_steps": 50}, core=core)
    assert proc.inputs() == {}
    out = proc.update({}, 1.0)
    assert set(out) == {"oxygen_deficit", "total_infection", "cytokines",
                        "cell_populations", "alive"}
    assert set(out["cytokines"]) == set(CYTOKINES)


def test_composite_runs_and_emits():
    from viva_iirabm.composites.sepsis import iirabm_baseline
    core = build_core()
    sim = Composite({"state": iirabm_baseline(core, seed=1, n_steps=120)}, core=core)
    sim.run(120)
    rows = gather_emitter_results(sim)[("emitter",)]
    assert len(rows) == 121
    # replayed deltas reconstruct a positive absolute oxygen deficit
    assert rows[-1]["oxygen_deficit"] > 0


def test_intervention_applies():
    """A treated run should differ from the untreated run (same seed)."""
    from viva_iirabm.interventions import FIG2_DIGITIZED
    base = run_sim(seed=1, n_steps=1500)["oxygen_deficit"]
    treated = run_sim(seed=1, n_steps=1500, intervention=FIG2_DIGITIZED)["oxygen_deficit"]
    assert not np.array_equal(base, treated)


@pytest.mark.slow
def test_mortality_regime():
    """The trained parameter set is a high-mortality sepsis case (paper ~82%)."""
    from viva_iirabm.runtime import mortality_rate
    mr = mortality_rate(range(1, 21), n_steps=9000)
    assert mr > 0.5
