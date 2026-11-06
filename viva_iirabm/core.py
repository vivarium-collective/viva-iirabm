"""Explicit core builder so the workspace's own process is always registered
regardless of how the editable install was made (see discovery convention)."""
from __future__ import annotations

from process_bigraph import allocate_core

from .processes import IIRABMProcess


def build_core(core=None):
    if core is None:
        core = allocate_core()
    core.register_link("IIRABMProcess", IIRABMProcess)
    return core
