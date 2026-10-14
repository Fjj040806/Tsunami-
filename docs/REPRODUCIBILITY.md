# Reproducibility and metric definitions

All numerical parameters are fixed in JSON configurations. Units are SI unless explicitly normalized. Source hashes and measured runtimes for the final step, shelf and laboratory computations are in `results/final_*_run_manifest.json`. The release manifest separately records delivered-file hashes. Historical outputs are preserved under versioned directories.

## Study inventory

| Study | Parameters | Meshes | Purpose |
|---|---|---|---|
| Step | D=20/5 m, A=0.002 m, Gaussian center −1200 m, σ=120 m, t=220 s | 500, 1000, 2000, 4000 | Independent linear reference and coupling control |
| Shelf | Ddeep=20 m; ratios 0.5/0.25; widths 200/800/2400 m; A/Ddeep=0.001/0.01/0.05; σ=150 m; center −2000 m; t=450 s | 1000, 2000, 4000 | 18 cases, both models |
| Focused weak refinement | A/Ddeep=0.001, width=2400 m, both depth ratios | 4000, 8000, 16000 | Exploratory follow-up after standard screen failed |
| Laboratory | NOAA Case A, G4 input; G5–G10 comparison; original time 265.05–283 s | 400, 800, 1600 | Conditional fixed-parameter comparison |
| Offline animation | D=20/5 m, A=0.2 m, widths 200/2400 m, 61 frames over 450 s | 1000 | Illustration; not quantitative mesh evidence |

The final recorded step, shelf and laboratory stages took approximately 41, 391 and 22 seconds respectively in the build environment; these are observed runtimes, not hardware-independent promises. Fine refinement, earlier verification studies, plotting and animation add computation. Full reproduction is available through `scripts/reproduce_all.py`; its component commands were executed individually during project construction.

## Metrics

- Step `relative_L2_to_linear_reference`: discrete L2 discrepancy relative to the independent cell-averaged linear solution. For a nonlinear model, this includes physical and numerical differences; it is not a pure nonlinear numerical error.
- Shelf `relative_L2_model_difference`, E: L2 difference between nonlinear and linear surface, divided by the linear surface norm at the same final time. It is not an error normalized by initial amplitude and is not a local percentage height error.
- Shelf `combined_mesh_sensitivity`, S: sum of each model's fine-to-restricted-fine/coarse discrepancy, divided by the fine linear norm. See the exact restriction convention in the report and script. E/S≥5 is an empirical resolution screen, not an error bound or confidence interval.
- Shelf transmitted and reflected peak ratios: maxima in prescribed spatial masks, divided by A. Green's depth-ratio law is a slow-variation reference, not an exact solution for every transition.
- Laboratory `NRMSE_over_observed_peak`: RMS discrepancy divided by the observed maximum absolute baseline-corrected surface over the fixed window. Denominators are gauge-specific.
- Laboratory `normalized_two_grid_RMS_change`: temporal RMS change between the last two meshes divided by the same observed peak. It is numerical sensitivity, not measurement uncertainty.
- Laboratory `peak_time_difference_s`: time difference between global positive peaks. Multiple peaks can switch branches; this is not a robust first-arrival estimator.
- Mass residuals include integrated boundary flux where applicable. Roundoff-level residuals support conservation implementation; they do not certify all other errors.

## Verification and release checks

`python scripts/run_verification.py` executes all 35 tests without external data. Synthetic analytic, manufactured, still-water, mass, stage-CFL, incoming-boundary, wall and stationary-junction tests are included. The scientific studies additionally monitor minimum water depth and actual CFL. Smooth-order evidence does not imply second-order convergence at shocks. No universal positivity theorem is claimed.

The executed notebook is a lightweight companion: it reads supplied results and independently runs a lake-at-rest check. Its cells were executed in sequence in one shared Python namespace with captured figures and stdout; execution metadata records this method. It can be rerun using a conventional Jupyter kernel after installing the optional notebook requirements.

The delivered figures/tables are generated from CSV/JSON results, rather than manually transcribed. Raw laboratory observations must be fetched to regenerate comparison figures. The HTML explorer includes precomputed modest-grid data and needs no server or internet connection.

## Limits

The source implements a bounded computational study, not a general coastal tsunami package. The restricted Bernoulli interface is used for an aligned, piecewise-constant, fully wet subcritical bed step only. Smooth-bed runs use centered bed reconstruction and hydrostatic source balance. The linear and nonlinear models intentionally share the initial surface and velocity, rather than exactly identical nonlinear characteristic invariants. At very weak amplitude, differences may include this initialization effect and accumulated numerical phase differences.

No wetting/drying, dispersive correction, historical event source inversion, calibrated friction, physical breaking closure or statistical observation-error model is included. Four standard weak-wave scenarios remain below the empirical resolution threshold. Laboratory metadata discrepancies and boundary/geometry assumptions prevent an unconditional validation claim.
