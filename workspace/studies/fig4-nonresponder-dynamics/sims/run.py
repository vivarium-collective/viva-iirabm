#!/usr/bin/env python3
"""Reproduce Fig 4: responder-average vs non-responder dynamics (oxygen deficit,
total GCSF, neutrophil population) under the GA-identified intervention."""
from viva_iirabm.reproduce import fig4_data, fig4_figure
from viva_iirabm.study_runner import execute

N_SEEDS = 40
N_STEPS = 12500


def build():
    data = fig4_data(n_seeds=N_SEEDS, n_steps=N_STEPS)
    fig = fig4_figure(data)
    result = {
        "n_seeds": N_SEEDS, "n_steps": N_STEPS,
        "base_mortality": data["base_mortality"],
        "treated_mortality": data["treated_mortality"],
        "n_responders": len(data["responders"]),
        "n_nonresponders": len(data["nonresponders"]),
        "headline": (f"treated mortality {100*data['treated_mortality']:.0f}% vs "
                     f"baseline {100*data['base_mortality']:.0f}%; "
                     f"{len(data['nonresponders'])} non-responders"),
    }
    return {"nonresponder-dynamics": fig}, result


if __name__ == "__main__":
    execute("fig4-nonresponder-dynamics", build)
