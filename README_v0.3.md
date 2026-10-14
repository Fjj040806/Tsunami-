# Tsunami Modeling — V0.3 foundation

Variable-depth model derivation and verified fully wet solvers, preparing for
controlled shoaling/reflection experiments. No tsunami observation validation
or coastal inundation claim is included in this release.

## Run

Python 3.10 or later, from this directory:

```bash
python -m pip install -r requirements.txt
python scripts/run_verification.py
python scripts/run_v03.py
```

Tests use standard-library unittest. No external data, compiler, notebook,
or display is needed. CSV/JSON/NPZ results and PNG figures are generated in
`results/` and `figures/`. On Windows, `py` can replace `python`.
Optional install: `python -m pip install -e .`.

## Implemented

- Independent variable-depth linear solver for `(eta,u)`, with an exact
  stationary-interface flux enforcing common surface and discharge.
- Nonlinear MUSCL hydrostatic reconstruction, HLL flux, centered intra-cell
  bed source and SSP-RK2. Pressure/source balance follows Audusse et al. (2004).
- Smooth-bed centered reconstruction; explicit constant-bed mode for aligned
  steps. No claim of second-order accuracy near discontinuous beds.
- Strict fully wet checks, stage CFL monitoring and weighted boundary/source
  mass accounting. No depth clipping or general positivity/entropy guarantee.
- Manufactured analytic fields and sources for algorithm verification only.

27 automatic tests pass. A 100-second static smooth-bed test gives maximum
nonlinear surface drift of 1.35e-14 m. Smooth variable-bed manufactured tests
show near-second-order mean L1 convergence. Forced test fields are not physical
sea-wave observations, and do not establish equivalence of the two models.

## Read and inspect

| Path | Contents |
|---|---|
| `docs/V03_VARIABLE_DEPTH_CN.md` | Derivation, balance proof, measured verification, limits |
| `paper/methods_v0.3.tex` | Standalone English methods/verification draft |
| `src/tsunami/variable.py` | Both new solvers and interface flux |
| `src/tsunami/manufactured.py` | Independent analytic test problem |
| `tests/test_variable.py` | Eight new numerical and contract tests |
| `configs/variable_verification.json` | Verification settings |
| `results/` | Tables, field data and source provenance |
| `figures/` | Static balance and manufactured convergence figures |

V0.1 and V0.2 code, derivations and historical results are retained. Their
readmes describe their original releases. `run_experiments.py` and `run_v02.py`
can regenerate those experiments in the current results/figures directories;
versioned filenames keep the new V0.3 results distinct. Historical manifests
refer to the source versions used for those runs.

Smooth-bed reconstruction can overshoot steep changes. Choose
`bed_reconstruction='constant'` for an aligned discontinuous step; always check
resolution. Outflow ghost copying is not an exact nonreflecting condition.
Optional forcing callbacks supply cell-average `(N,2)` arrays at each RK time;
use them only for explicit forced problems, not ordinary tsunami propagation.

Next: unforced pulse scattering across an aligned step, followed by a smooth
continental shelf and resolution-controlled linear/nonlinear comparisons.
External data validation and final paper assembly remain later milestones.

The methods are established, and implementation is original. Citations are in
`paper/references.bib` and the Chinese report. Original software is MIT licensed.
The complete project is ready for a GitHub repository; no GitHub push occurred.
