#!/usr/bin/env python3
"""Reproduce Fig 2 + Table 1: run the (reduced-scale) genetic algorithm against
the real IIRABM to discover a multi-stage cytokine intervention, show its
convergence and mortality reduction, and the mortality-vs-#interventions curve.
Also renders the digitized Fig 2 intervention for reference."""
from viva_iirabm.ga import run_ga
from viva_iirabm.interventions import FIG2_DIGITIZED
from viva_iirabm.reproduce import (fig2_figure, ga_convergence_figure,
                                   table1_data, table1_figure)
from viva_iirabm.runtime import mortality_rate
from viva_iirabm.study_runner import execute

# reduced-scale GA (paper used ~1e6 CPU-hours on supercomputers)
GA_KW = dict(n_stages=3, pop_size=10, generations=5,
             replicate_seeds=(1, 2, 3), n_steps=3500, ga_seed=0)
EVAL_SEEDS = 40
EVAL_STEPS = 9000


def build():
    best, best_f, history = run_ga(**GA_KW,
                                   progress=lambda g, f: print(f"  GA gen {g}: fitness {f:.0f}"))
    # mortality of the GA-found intervention vs untreated baseline
    base_mr = mortality_rate(range(1, EVAL_SEEDS + 1), n_steps=EVAL_STEPS)
    ga_mr = mortality_rate(range(1, EVAL_SEEDS + 1), n_steps=EVAL_STEPS, intervention=best)
    # Table-1-style: mortality vs number of sequential interventions (digitized genome)
    t1 = table1_data(seeds_n=EVAL_SEEDS, n_steps=EVAL_STEPS, genome=FIG2_DIGITIZED)

    figs = {
        "ga-convergence": ga_convergence_figure(history),
        "ga-intervention-bars": fig2_figure(best),
        "fig2-digitized-bars": fig2_figure(FIG2_DIGITIZED),
        "mortality-vs-interventions": table1_figure(t1),
    }
    result = {
        "ga": GA_KW, "best_fitness": best_f, "fitness_history": history,
        "best_intervention": best,
        "baseline_mortality": base_mr,
        "ga_intervention_mortality": ga_mr,
        "table1": t1["rows"],
        "headline": (f"GA lowered mortality {100*base_mr:.0f}% -> {100*ga_mr:.0f}% "
                     f"(reduced-scale; paper 82% -> 16% at full scale)"),
    }
    return figs, result


if __name__ == "__main__":
    execute("fig2-ga-intervention", build)
