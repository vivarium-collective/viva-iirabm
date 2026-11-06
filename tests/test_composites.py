"""Composite-generator registration + intervention-encoding tests."""
from viva_iirabm.interventions import ALLOWED, FIG2_DIGITIZED, TARGETS, flatten


def test_generators_registered():
    from process_bigraph.composite_generator import _REGISTRY
    for gen in ("iirabm_baseline", "iirabm_treated"):
        assert any(eid.endswith("." + gen) or eid == gen for eid in _REGISTRY), \
            f"{gen} not registered; have {list(_REGISTRY)[:8]}"


def test_fig2_genome_shape_and_values():
    assert len(FIG2_DIGITIZED) == 8                 # 8 stages
    for stage in FIG2_DIGITIZED:
        assert len(stage) == len(TARGETS) == 12     # 12 cytokine targets
        for v in stage:
            assert v in ALLOWED
    assert len(flatten(FIG2_DIGITIZED)) == 96


def test_generators_discoverable():
    from process_bigraph.composite_generator import discover_generators
    gens = discover_generators()
    assert any("viva_iirabm" in g for g in gens)
