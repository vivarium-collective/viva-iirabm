# Contributing to viva-iirabm

## Development setup

uv is required (`brew install uv` or `pip install uv`).

    uv venv .venv
    source .venv/bin/activate
    uv pip install -e ".[dev]"
    pytest                    # fast tests (slow sims deselected by default)
    pytest -m slow            # the full mortality-regime check (~minutes)

The IIRABM C++ binary is built automatically on first use (needs a C++11
compiler: `c++` / `clang++` / `g++`; no MPI). It is cached under
`vendor/_build/` (or `$VIVA_IIRABM_BUILD_DIR`).

## Reproducing the paper figures

    ./scripts/reproduce_all.sh          # runs all four study runners
    # or one at a time:
    python workspace/studies/<slug>/sims/run.py

Each runner writes an interactive figure to `workspace/studies/<slug>/viz/`
and a `data/result.json` with the headline metrics.

## Vendored upstream

`vendor/iirabm_upstream/` is the pristine An & Cockrell IIRABM (GPLv2+).
`vendor/patch_iirabm.py` applies the reproducible serial port. **Do not edit
the vendored sources directly** — change the patch script and rebuild, so the
provenance stays clean. See `vendor/PROVENANCE.md`.

## Releasing to PyPI

Tag a commit with `git tag v<VERSION>` and push the tag; the release workflow
publishes to PyPI via trusted publishing. One-time PyPI setup:
https://docs.pypi.org/trusted-publishers/
