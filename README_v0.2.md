# Tsunami Modeling — V0.2

A reproducible comparison of linear and nonlinear one-dimensional shallow-water
models, with independent exact-solution verification and explicit mesh checks.
This is a research project in progress, not a calibrated tsunami prediction tool.

## Run

Python 3.10 or later; from this directory:

```bash
python -m pip install -r requirements.txt
python scripts/run_verification.py
python scripts/run_v02.py
```

Tests use standard-library unittest. Experiments run without a display or external
data and write CSV, JSON, NPZ, PNG and GIF outputs into `results/` and `figures/`.
On Windows, `py` can replace `python`. An optional editable install is
`python -m pip install -e .`. Tested direct dependency versions are recorded in
`requirements-tested.txt`; this is not a full transitive lockfile.

## Models and scope

- `flat.py`: independent exact-upwind linear flux and nonlinear HLL flux;
  componentwise MC-limited MUSCL reconstruction and SSP-RK2, constant depth only.
- `solver.py`: retained V0.1 first-order hydrostatic reconstruction for a static
  variable bed, including still-water preservation.
- Fully wet states only; no evolved-depth clipping or shoreline wet/dry method.
- Periodic smooth-wave comparison; no coast, source inversion, dispersion,
  observation fitting, or real-event validation.

Both models use the same continuous initial surface and velocity. Their respective
cell-averaged evolved variables differ: `(eta,u)` versus `(h,hu)`. Comparison uses
six incident amplitude/depth ratios, four nested meshes, and an additional 6400-cell
mesh for the smallest amplitude. The primary result is RMS model difference
normalized by incident amplitude, checked against last-two-grid sensitivity.
This sensitivity is empirical, not a rigorous error bound.

## Measured results

19 automatic tests pass. Smooth exact-reference tests show approximately second-order
mean L1 convergence (finest observed orders 1.996 and 1.998); this does not apply
to every norm or the retained variable-bed solver. At `t=0.2 L/sqrt(gD)`, waveform
RMS differences are 0.625%, 3.126%, 6.252% and 12.505% of the incident amplitude
for `A/D=0.01, 0.05, 0.1, 0.2`, respectively. These are model comparisons under
one prescribed setup, not measurement errors or universal validity thresholds.
Small-amplitude crest locations remain mesh sensitive.

## Read and inspect

| Path | Content |
|---|---|
| `docs/V02_LINEAR_NONLINEAR_CN.md` | Chinese derivation, experiment design, results, limitations |
| `paper/results_v0.2.tex` | English standalone methods/results draft |
| `configs/flat_comparison.json` | Fixed experiment parameters |
| `tests/` | 19 numerical and input-contract tests |
| `figures/` | Convergence, waveform comparison, amplitude effects, GIF |
| `results/` | Full tables, finest profiles, diagnostics and provenance |
| `docs/MODEL_DERIVATION_CN.md` | Historical V0.1 model derivation |
| `v0.1_results/`, `v0.1_figures/` | Historical first-order experiments |

To regenerate the historical V0.1 experiments, `python scripts/run_experiments.py`
writes its own outputs to the current results/figures directories; filenames differ
from V0.2. `README_v0.1.md` describes that historical release; its paths and test
count refer to the original package.

The numerical methods are established; the implementation and experiments are
original. Method citations are in `paper/references.bib` and the Chinese report.
Original software is MIT licensed. No third-party observation data are bundled.
The package is ready for a GitHub repository but has not been uploaded to GitHub.

Next milestone: a verified second-order well-balanced variable-bed method and a
variable-depth linear solver, followed by controlled shoaling/reflection studies.
