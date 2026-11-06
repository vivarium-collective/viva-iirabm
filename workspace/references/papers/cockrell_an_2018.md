# Cockrell & An (2018) — reference notes

**Citation.** Cockrell RC, An G. *Examining the controllability of sepsis using
genetic algorithms on an agent-based model of systemic inflammation.*
PLoS Comput Biol 14(2): e1005876 (2018).
https://doi.org/10.1371/journal.pcbi.1005876

**Model.** Innate Immune Response ABM (IIRABM) — a 2-D agent-based model of the
endothelial–blood interface (101×101 torus) simulating the innate inflammatory
response. Agents: endothelial cells, macrophages/monocytes, neutrophils (+marrow
precursors), TH0/TH1/TH2 (+germ precursors). Cytokines tracked: PAF, TNF, sTNFr,
IL1, sIL1r, IL1ra, IFNγ, IL4, IL8, IL10, IL12, GCSF. The system dies when
aggregate endothelial damage (oxygen deficit) exceeds 80% (`oxyDeficit > 8161`).

**Five initialization parameters.** injury size (radius, cells), microbial
invasiveness, microbial toxigenesis, environmental toxicity, host resilience.
Trained set: invasiveness=2, toxigenesis=5, resilience=0.1, env-toxicity=2,
injury=33 → ~82% mortality (general case), 68% for the specific patient at 12 h.

**Control.** An *intervention* augments/inhibits cytokine production over
sequential 12-h stages (each stage = 12 per-cytokine multipliers from
{0.05,0.25,0.5,0.75,1,2,5,10,20}; inhibition multiplicative, augmentation
additive). A genetic algorithm (tournament selection size 2, uniform crossover)
minimizes fitness F = Σ(oxygen deficit + infection) over stochastic replicates.

**Key figures.**
- **Fig 2 / Table 1** — the best GA-found 8-stage intervention; 82%→16% (general),
  68%→12% (specific). Table 1: mortality falls as #interventions rises.
- **Fig 3** — generalizability across 8,291 clinically-relevant parameter sets.
- **Fig 4** — responder vs non-responder dynamics (oxygen deficit, GCSF, neutrophils);
  non-responders over-produce GCSF → surplus neutrophils → sustained inflammation.
- **Fig 6** — 100 stochastic trajectories (red death / blue heal); 82% / 68%.

**Compute.** A single 90-day sim ≈ 4 min; each GA experiment up to ~1,000,000
CPU-hours (≈ $100k) on Cray XC30/XE6 supercomputers.

**Source.** IIRABM C++ (GPLv2+): https://bitbucket.org/cockrell/iirabm_public
(vendored commit `eaf7e08` under `vendor/`).
