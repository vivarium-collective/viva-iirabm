#!/usr/bin/env python3
"""Reproduce Fig 3 (subsampled): sweep a sample of clinically-relevant parameter
sets, compute base vs post-intervention mortality for each, and plot the
base-vs-post scatter (below-diagonal = intervention beneficial).

The paper sweeps 8,291 parameter sets x 100 replicates (supercomputer-scale);
this subsamples a grid at reduced replicate count and is explicit about it."""
import itertools

from viva_iirabm.reproduce import fig3_data, fig3_figure
from viva_iirabm.study_runner import execute

SEEDS_N = 15
N_STEPS = 9000

# a subsample across the clinically-relevant region (invasiveness x toxigenesis x
# resilience x env-toxicity), centered on the paper's trained set.
INVAS = [1, 2, 3]
TOXI = [4, 5, 6]
RESIL = [0.05, 0.1]
ENVTOX = [2]


def _param_sets():
    for inv, tox, res, env in itertools.product(INVAS, TOXI, RESIL, ENVTOX):
        yield {"invasiveness": inv, "toxigenesis": tox,
               "host_resilience": res, "env_toxicity": env, "injury": 33}


def build():
    param_sets = list(_param_sets())
    data = fig3_data(param_sets, seeds_n=SEEDS_N, n_steps=N_STEPS)
    fig = fig3_figure(data)
    pts = data["points"]
    improved = sum(1 for p in pts if p["post_mr"] < p["base_mr"])
    result = {
        "n_param_sets": len(pts), "seeds_n": SEEDS_N, "n_steps": N_STEPS,
        "improved": improved,
        "mean_base_mr": sum(p["base_mr"] for p in pts) / len(pts),
        "mean_post_mr": sum(p["post_mr"] for p in pts) / len(pts),
        "points": pts,
        "headline": f"{improved}/{len(pts)} parameter sets improved by the intervention",
    }
    return {"generalizability-scatter": fig}, result


if __name__ == "__main__":
    execute("fig3-generalizability", build)
