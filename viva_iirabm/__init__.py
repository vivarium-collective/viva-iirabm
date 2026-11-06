"""viva-iirabm: process-bigraph wrapper for the Innate Immune Response ABM.

Real bridge to the An & Cockrell IIRABM C++ simulator, used to reproduce
Cockrell & An (2018), "Examining the controllability of sepsis using genetic
algorithms on an agent-based model of systemic inflammation" (PLoS Comput Biol).
"""
from .core import build_core
from .processes import IIRABMProcess
from .composites import iirabm_baseline, iirabm_treated

__all__ = ["IIRABMProcess", "build_core", "iirabm_baseline", "iirabm_treated"]
