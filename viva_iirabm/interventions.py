"""Named intervention genomes for the IIRABM.

``FIG2_DIGITIZED`` is the 8-stage, GA-identified intervention shown in Fig 2 of
Cockrell & An (2018), read off the figure's bar heights (base-2 log of the
per-cytokine augmentation/inhibition multiplier) and mapped to the nearest
value in the model's allowed set {0.05, 0.25, 0.5, 0.75, 1, 2, 5, 10, 20}.

This is a *digitization* of the published figure (the paper reports the
intervention only as bar heights, not as a numeric table), so treat it as an
approximate reconstruction, not the authors' exact genome. The GA module
(``viva_iirabm.ga``) reproduces the *method* that produced such an intervention.
"""
from __future__ import annotations

# target order within each 12-value stage
TARGETS = ["PAF", "TNF", "sTNFr", "IL1", "sIL1r", "IL1ra",
           "IFNg", "IL4", "IL8", "IL10", "IL12", "GCSF"]

# allowed multiplier set (IIRABM possibleInterventions)
ALLOWED = [0.05, 0.25, 0.5, 0.75, 1.0, 2.0, 5.0, 10.0, 20.0]
# log2 of each allowed value, for mapping bar heights
import math
_LOG2 = {round(math.log2(v), 3): v for v in ALLOWED}


def from_log2_bars(bars_by_stage: list[dict[str, float]]) -> list[list[float]]:
    """Map a list of {cytokine: log2(multiplier)} stages to multiplier genomes."""
    allowed_log2 = sorted(_LOG2)
    genome = []
    for stage in bars_by_stage:
        row = []
        for t in TARGETS:
            lv = stage.get(t, 0.0)
            nearest = min(allowed_log2, key=lambda a: abs(a - lv))
            row.append(_LOG2[nearest])
        genome.append(row)
    return genome


# Fig 2 bar heights (base-2 log multiplier) digitized per 12-hr stage.
# Bars at approx +4.3/-4.3 -> 20/0.05 ; +2.3/-2.3 -> 5/0.25 ; +1/-1 -> 2/0.5 ;
# ~0 -> 1 (no change). Values read from the published panels.
_FIG2_BARS = [
    # 12 hrs
    {"PAF": 4.3, "TNF": -2, "sTNFr": -4.3, "IL1": -4.3, "sIL1r": -1,
     "IL1ra": 2.3, "IFNg": 4.3, "IL4": 2.3, "IL8": -4.3, "IL10": 0, "IL12": 3.3, "GCSF": -4.3},
    # 24 hrs
    {"PAF": 4.3, "TNF": -0.4, "sTNFr": 2.3, "IL1": 0, "sIL1r": -4.3, "IL1ra": -0.4,
     "IFNg": -2, "IL4": 1, "IL8": -4.3, "IL10": -1, "IL12": 4.3, "GCSF": -4.3},
    # 36 hrs
    {"PAF": 4.3, "TNF": -2, "sTNFr": -4.3, "IL1": -4.3, "sIL1r": -4.3, "IL1ra": -1,
     "IFNg": 0, "IL4": -2, "IL8": 0, "IL10": -4.3, "IL12": -0.4, "GCSF": 4.3},
    # 48 hrs
    {"PAF": 4.3, "TNF": 1, "sTNFr": 1, "IL1": -0.4, "sIL1r": 1, "IL1ra": 3.3,
     "IFNg": 4.3, "IL4": -2, "IL8": -2, "IL10": 0, "IL12": 1, "GCSF": -2},
    # 60 hrs
    {"PAF": 4.3, "TNF": -4.3, "sTNFr": 4.3, "IL1": -4.3, "sIL1r": 1, "IL1ra": 0,
     "IFNg": 3.3, "IL4": 2.3, "IL8": -2, "IL10": 1, "IL12": -4.3, "GCSF": -4.3},
    # 72 hrs
    {"PAF": 4.3, "TNF": -0.4, "sTNFr": 4.3, "IL1": -4.3, "sIL1r": 3.3, "IL1ra": -1,
     "IFNg": 0, "IL4": 2.3, "IL8": -1, "IL10": -1, "IL12": -1, "GCSF": -2},
    # 84 hrs
    {"PAF": 4.3, "TNF": -2, "sTNFr": 1, "IL1": -4.3, "sIL1r": 4.3, "IL1ra": -0.4,
     "IFNg": -0.4, "IL4": 2.3, "IL8": -0.4, "IL10": -0.4, "IL12": -1, "GCSF": -4.3},
    # 96 hrs
    {"PAF": 4.3, "TNF": -2, "sTNFr": 3.3, "IL1": 4.3, "sIL1r": -4.3, "IL1ra": 0,
     "IFNg": 2.3, "IL4": -0.4, "IL8": -4.3, "IL10": 1, "IL12": -1, "GCSF": 1},
]

FIG2_DIGITIZED = from_log2_bars(_FIG2_BARS)


def flatten(genome: list[list[float]]) -> list[float]:
    return [v for stage in genome for v in stage]


NAMED = {
    "fig2_digitized": FIG2_DIGITIZED,
}
