# Vendored upstream: IIRABM (Innate Immune Response Agent-Based Model)

- **Source:** https://bitbucket.org/cockrell/iirabm_public
- **Vendored commit:** `eaf7e08b1774ba3860fed393e7484ab9cdc4076f` (branch `master`)
- **Authors:** Gary An & R. Chase Cockrell
- **License:** GNU General Public License v2 or later (GPLv2+). See
  `iirabm_upstream/UPSTREAM_README.md` for the upstream copyright notice.
- **Paper reproduced:** Cockrell RC, An G. *Examining the controllability of
  sepsis using genetic algorithms on an agent-based model of systemic
  inflammation.* PLoS Comput Biol 14(2): e1005876 (2018).
  https://doi.org/10.1371/journal.pcbi.1005876

## What `patch_iirabm.py` does

The upstream `V1/` targets a Cray/GCC toolchain with the deprecated **MPI C++
bindings** (`MPI::Init`, `MPI::COMM_WORLD.*`), which were removed in MPI-3 and
are absent from modern Open MPI / Apple clang. `patch_iirabm.py` applies a
small, reproducible transform to a *copy* of the pristine sources:

1. Adds `#include <random>` to the translation units that use `std::mt19937`
   but relied on transitive includes under the old toolchain; drops the
   `mpi.h` includes.
2. Converts a bounded set of **run-mode / timing constants**
   (`numTimeSteps`, `GA`, `testIntervention`, the intervention-schedule knobs,
   reseed knobs) from compile-time `const` in `Parameters.h` into **runtime
   globals**, so one binary covers baseline, intervention-replay and variable
   intervention schedules without recompilation. Array-sizing constants
   (`genomeLength`, `ga_popSize`, `numWriteSteps`, grid dims) stay compile-time.
3. Replaces the MPI orchestration `main()` / `postSimulationData()` with a
   **serial single-run driver** that takes the parameter set, seed and
   intervention from `argv` and writes the full per-timestep trajectory to a
   file.

**The physics is byte-for-byte upstream.** `mainSimulation()` (the agent loop)
and the cell / simulation / agent-construction functions are unchanged — only
the MPI *orchestration* layer is replaced. The paper's MPI build parallelised
by making each rank an independent simulation; this port parallelises the same
way, at the Python level (one OS process per replicate).

The combined compiled binary is a derivative work of GPLv2+ upstream code and
is therefore GPLv2+. The `viva_iirabm` Python package (which invokes the binary
as a separate process) is MIT-licensed; see the repo `LICENSE`.
