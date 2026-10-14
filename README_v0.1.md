# Tsunami Modeling — V0.1

**Model derivation and a verified minimal solver.** One-dimensional nonlinear
shallow-water equations, a static piecewise-constant bed, first-order finite
volumes, an HLL flux with endpoint/Roe speed bounds, hydrostatic reconstruction,
and adaptive explicit Euler time stepping.

This release establishes a computational foundation. It does not yet contain
a separate linear-model solver, parameter study, shoreline inundation, or
experimental/real-event validation. The final research report remains future work.

## Start here

Chinese derivation and walkthrough: `docs/MODEL_DERIVATION_CN.md`.
Measured verification results: `docs/VERIFICATION_RESULTS.md`.

Use Python 3.10 or later. From this extracted directory:

```bash
python -m pip install -r requirements.txt
python scripts/run_verification.py
python scripts/run_experiments.py
```

On Windows, `py` may be used in place of `python`. No compiler, external data,
Jupyter installation, or editable package install is required. Tests use the
Python standard-library unittest runner. The experiment script writes to
`results/` and `figures/`; plotting does not require a display. To skip the GIF:

```bash
python scripts/run_experiments.py --no-animation
```

For a short script with a plot window:

```bash
python examples/minimal_demo.py
```

Optional package install for importing in your own work: `python -m pip install -e .`.

## Included results

- Static free-surface preservation on smooth and stepped topography.
- Periodic, reflective-wall, and boundary-flux-adjusted open-domain mass checks.
- Small-amplitude phase-speed comparison with `sqrt(gD)`.
- Grid convergence against an independent smooth nonlinear characteristic solution.
- Wet-bed dam-break comparison with its exact entropy solution (test suite).
- A computed Gaussian pulse, snapshots, propagation plot, and GIF.

The first-order method damps wave amplitudes numerically. Do not interpret this
as calibrated physical dissipation. Zero-gradient outflow boundaries are not
exactly nonreflecting; the propagation demo keeps the pulse away from boundaries.

The solver stops if any physical cell reaches the minimum-depth threshold. It
never silently clips the evolved water depth. The nonnegative interface-depth
truncation belongs to hydrostatic reconstruction, not a wet/dry model.

## Files

| Directory | Contents |
|---|---|
| `src/tsunami/` | solver, initial conditions, independent exact references |
| `tests/` | numerical and input-contract tests |
| `scripts/` | portable tests and reproducible experiments |
| `examples/` | minimal editable pulse example |
| `docs/` | Chinese derivation and measured results |
| `paper/` | English methods draft and BibTeX for the future report |
| `figures/` | computed scientific plots and animation |
| `results/` | CSV errors, JSON diagnostics/manifest, NPZ snapshots |

The manifest records configurations, source hashes and observed dependency
versions. `requirements-tested.txt` pins the three direct dependencies from the
tested environment; it is not a complete cross-platform transitive lockfile.

## Attribution and next stage

Hydrostatic reconstruction follows Audusse et al. (2004),
DOI [10.1137/S1064827503431090](https://doi.org/10.1137/S1064827503431090).
The implementation is original; the numerical method and governing model are
established. See `paper/references.bib` and the derivation for scientific sources.
No third-party source code or observation data are copied into this release.

Next: establish a separately verified linear solver, quantify boundary
reflection, then introduce fully wet shelf topography and compare models.
High-order and wet/dry handling are separate later milestones.

Original software is MIT licensed. This package has not been pushed to GitHub.
