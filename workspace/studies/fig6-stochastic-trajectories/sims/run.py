#!/usr/bin/env python3
"""Reproduce Fig 6 / baseline mortality: 100 stochastic IIRABM trajectories for
the paper's trained parameter set (~82% death)."""
from viva_iirabm.reproduce import fig6_data, fig6_figure
from viva_iirabm.study_runner import execute

N_SEEDS = 100
N_STEPS = 12500


def build():
    data = fig6_data(n_seeds=N_SEEDS, n_steps=N_STEPS)
    fig = fig6_figure(data)
    result = {
        "n_seeds": data["n_seeds"], "n_steps": data["n_steps"],
        "deaths": data["deaths"], "mortality": data["mortality"],
        "paper_target_mortality": 0.82,
        "headline": f"mortality {100*data['mortality']:.0f}% (paper ~82%)",
    }
    return {"stochastic-trajectories": fig}, result


if __name__ == "__main__":
    execute("fig6-stochastic-trajectories", build)
