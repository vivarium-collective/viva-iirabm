"""Reduced-scale reproduction of the paper's genetic-algorithm control search.

Cockrell & An (2018) train a GA over multi-stage cytokine interventions,
minimising fitness F = sum_i (OxygenDeficit_i + Infection_i) across stochastic
replicates, using tournament selection (size 2) + uniform crossover, with
per-gene values drawn from {0.05,0.25,0.5,0.75,1,2,5,10,20}. Their production
runs used up to ~1e6 CPU-hours on supercomputers. This module reproduces the
*method* faithfully against the real IIRABM at laptop scale (small population,
few generations, few replicate seeds, shorter sims) — enough to show the GA
discovers an effective multi-stage intervention, not the exact published genome.
"""
from __future__ import annotations

import random

from .interventions import ALLOWED, TARGETS, flatten
from .runtime import run_batch

N_TARGETS = len(TARGETS)  # 12


def random_stage(rng: random.Random) -> list[float]:
    return [rng.choice(ALLOWED) for _ in range(N_TARGETS)]


def random_genome(n_stages: int, rng: random.Random) -> list[list[float]]:
    return [random_stage(rng) for _ in range(n_stages)]


def evaluate(genome, *, seeds, n_steps, fitness_step=None, **sim_kw) -> float:
    """Mean fitness (oxy + infection) at ``fitness_step`` across replicate seeds.

    Lower is better. ``fitness_step`` defaults to the end of the run (the paper
    evaluates fitness 12 h after the last intervention)."""
    runs = run_batch(
        seeds, n_steps=n_steps,
        intervention=genome,
        **sim_kw,
    )
    vals = []
    for r in runs:
        oxy = r["oxygen_deficit"]
        infe = r["total_infection"]
        k = (fitness_step if fitness_step is not None else len(oxy)) - 1
        k = min(k, len(oxy) - 1)
        o = oxy[k]
        o = 8161.0 if o < 0 else o          # dead -> max damage
        i = infe[k]
        i = 0.0 if i < 0 else i
        vals.append(o + i)
    return sum(vals) / len(vals)


def tournament(pop, fits, rng, k=2):
    idx = rng.sample(range(len(pop)), k)
    return pop[min(idx, key=lambda i: fits[i])]


def uniform_crossover(a, b, rng):
    fa, fb = flatten(a), flatten(b)
    child = [fa[i] if rng.random() < 0.5 else fb[i] for i in range(len(fa))]
    n = len(a)
    return [child[i * N_TARGETS:(i + 1) * N_TARGETS] for i in range(n)]


def mutate(genome, rng, rate=0.1):
    for stage in genome:
        for i in range(N_TARGETS):
            if rng.random() < rate:
                stage[i] = rng.choice(ALLOWED)
    return genome


def run_ga(
    *,
    n_stages: int = 3,
    pop_size: int = 12,
    generations: int = 6,
    replicate_seeds=(1, 2, 3),
    n_steps: int = 4000,
    mutation_rate: float = 0.1,
    ga_seed: int = 0,
    progress=None,
    **sim_kw,
):
    """Evolve an intervention. Returns (best_genome, best_fitness, history)."""
    rng = random.Random(ga_seed)
    pop = [random_genome(n_stages, rng) for _ in range(pop_size)]
    fits = [evaluate(g, seeds=replicate_seeds, n_steps=n_steps, **sim_kw) for g in pop]
    history = [min(fits)]
    best_i = min(range(len(pop)), key=lambda i: fits[i])
    best, best_f = pop[best_i], fits[best_i]
    if progress:
        progress(0, best_f)

    for gen in range(1, generations + 1):
        new = [best]  # elitism
        while len(new) < pop_size:
            p1 = tournament(pop, fits, rng)
            p2 = tournament(pop, fits, rng)
            child = uniform_crossover(p1, p2, rng)
            child = mutate(child, rng, mutation_rate)
            new.append(child)
        pop = new
        fits = [evaluate(g, seeds=replicate_seeds, n_steps=n_steps, **sim_kw) for g in pop]
        gi = min(range(len(pop)), key=lambda i: fits[i])
        if fits[gi] < best_f:
            best, best_f = pop[gi], fits[gi]
        history.append(best_f)
        if progress:
            progress(gen, best_f)

    return best, best_f, history
