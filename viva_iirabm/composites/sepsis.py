"""Composite generators for the IIRABM — the dashboard-visible entry points."""
from __future__ import annotations

from process_bigraph.composite_generator import composite_generator

from ..interventions import NAMED, flatten

_EMIT = {
    "oxygen_deficit": "float",
    "total_infection": "float",
    "cytokines": "map[string,float]",
    "cell_populations": "map[string,float]",
    "alive": "boolean",
}

_PARQUET_PATHS = [
    "stores/oxygen_deficit", "stores/total_infection",
    "stores/cytokines", "stores/cell_populations", "stores/alive",
]


def _document(core=None, *, intervention_flat=None, intervention_stages=0, **cfg):
    proc_cfg = dict(cfg)
    proc_cfg["intervention_flat"] = intervention_flat or []
    proc_cfg["intervention_stages"] = intervention_stages
    return {
        "iirabm": {
            "_type": "process",
            "address": "local:IIRABMProcess",
            "config": proc_cfg,
            "interval": 1.0,
            "outputs": {
                "oxygen_deficit": ["stores", "oxygen_deficit"],
                "total_infection": ["stores", "total_infection"],
                "cytokines": ["stores", "cytokines"],
                "cell_populations": ["stores", "cell_populations"],
                "alive": ["stores", "alive"],
            },
        },
        "stores": {
            "oxygen_deficit": 0.0,
            "total_infection": 0.0,
            "cytokines": {},
            "cell_populations": {},
            "alive": True,
        },
        "emitter": {
            "_type": "step",
            "address": "local:RAMEmitter",
            "config": {"emit": _EMIT},
            "inputs": {
                "oxygen_deficit": ["stores", "oxygen_deficit"],
                "total_infection": ["stores", "total_infection"],
                "cytokines": ["stores", "cytokines"],
                "cell_populations": ["stores", "cell_populations"],
                "alive": ["stores", "alive"],
            },
        },
    }


@composite_generator(
    name="iirabm_baseline",
    description="Untreated IIRABM sepsis trajectory (one seed). Defaults = the "
                "paper's trained 82%-death parameter set (I2/H0.1/T5/E2/inj33).",
    parameters={
        "injury": {"type": "integer", "default": 33, "description": "Injury radius (cells)"},
        "host_resilience": {"type": "float", "default": 0.1, "description": "Per-step oxy heal rate"},
        "toxigenesis": {"type": "integer", "default": 5, "description": "Microbial toxigenesis / infect spread"},
        "env_toxicity": {"type": "integer", "default": 2, "description": "Environmental toxicity (recurring injury)"},
        "invasiveness": {"type": "integer", "default": 2, "description": "Microbial invasiveness (infect repeat)"},
        "seed": {"type": "integer", "default": 1, "description": "RNG seed (stochastic replicate)"},
        "n_steps": {"type": "integer", "default": 8000, "description": "Simulation length (timesteps; 1 ~= 7 min)"},
    },
    emitters=[{"address": "local:ParquetEmitter", "paths": _PARQUET_PATHS}],
)
def iirabm_baseline(core=None, *, injury=33, host_resilience=0.1, toxigenesis=5,
                    env_toxicity=2, invasiveness=2, seed=1, n_steps=8000):
    return _document(core, injury=injury, host_resilience=host_resilience,
                     toxigenesis=toxigenesis, env_toxicity=env_toxicity,
                     invasiveness=invasiveness, seed=seed, n_steps=n_steps)


@composite_generator(
    name="iirabm_treated",
    description="IIRABM sepsis trajectory under a multi-stage cytokine intervention. "
                "'intervention' selects a named genome (e.g. the digitized Fig 2 "
                "8-stage GA solution).",
    parameters={
        "intervention": {"type": "string", "default": "fig2_digitized",
                         "description": "Named intervention genome: " + ", ".join(NAMED)},
        "seed": {"type": "integer", "default": 1, "description": "RNG seed"},
        "n_steps": {"type": "integer", "default": 8000, "description": "Simulation length (timesteps)"},
    },
    emitters=[{"address": "local:ParquetEmitter", "paths": _PARQUET_PATHS}],
)
def iirabm_treated(core=None, *, intervention="fig2_digitized", seed=1, n_steps=8000):
    genome = NAMED[intervention]
    return _document(core, seed=seed, n_steps=n_steps,
                     intervention_flat=flatten(genome),
                     intervention_stages=len(genome))
