"""First-order finite volumes: HLL flux + hydrostatic reconstruction.

The pressure corrections implement Audusse et al. (2004), DOI
10.1137/S1064827503431090. Implementation is original; the method is not.
Physical dry states are intentionally rejected in this release.
"""

from dataclasses import dataclass
from typing import Optional
import numpy as np


@dataclass(frozen=True)
class Grid:
    left: float
    right: float
    cells: int

    def __post_init__(self):
        if (not np.isfinite([self.left, self.right]).all()
                or self.right <= self.left or isinstance(self.cells, bool)
                or not isinstance(self.cells, (int, np.integer)) or self.cells < 3):
            raise ValueError("Grid needs finite left < right and integer cells >= 3")

    @property
    def dx(self):
        return (self.right - self.left) / self.cells

    @property
    def x(self):
        return self.left + (np.arange(self.cells) + 0.5) * self.dx


@dataclass(frozen=True)
class SolverConfig:
    g: float = 9.81
    cfl: float = 0.45
    boundary: str = "outflow"
    minimum_depth: float = 1e-10  # fail threshold, NOT a water-depth floor
    max_steps: int = 1_000_000

    def __post_init__(self):
        if not np.isfinite([self.g, self.cfl, self.minimum_depth]).all():
            raise ValueError("Solver parameters must be finite")
        if self.g <= 0 or not 0 < self.cfl <= 0.5 or self.minimum_depth <= 0:
            raise ValueError("Require g > 0, 0 < CFL <= 0.5, minimum_depth > 0")
        if self.boundary not in {"periodic", "outflow", "reflective"}:
            raise ValueError("Unknown boundary type")
        if isinstance(self.max_steps, bool) or not isinstance(self.max_steps, int) or self.max_steps < 1:
            raise ValueError("max_steps must be a positive integer")


@dataclass
class SimulationResult:
    x: np.ndarray
    bed: np.ndarray
    times: np.ndarray
    states: np.ndarray  # (outputs, cells, 2): depth h, discharge q=hu
    diagnostics: dict

    @property
    def eta(self):
        return self.states[:, :, 0] + self.bed[None, :]

    @property
    def velocity(self):
        return self.states[:, :, 1] / self.states[:, :, 0]


def _physical_flux(state, g):
    h, q = state[:, 0], state[:, 1]
    u = np.divide(q, h, out=np.zeros_like(q), where=h > 0)
    return np.column_stack((q, q * u + 0.5 * g * h * h))


def hll_flux(left, right, g):
    """HLL flux with endpoint and Roe/Einfeldt wave-speed bounds.

Inputs are reconstructed interface states, which may have zero depth.
Zero-depth reconstruction is not a physical wet/dry boundary capability.
"""
    hl, hr = left[:, 0], right[:, 0]
    ul = np.divide(left[:, 1], hl, out=np.zeros_like(hl), where=hl > 0)
    ur = np.divide(right[:, 1], hr, out=np.zeros_like(hr), where=hr > 0)
    cl, cr = np.sqrt(g * hl), np.sqrt(g * hr)
    roots = np.sqrt(hl) + np.sqrt(hr)
    uroe = np.divide(np.sqrt(hl) * ul + np.sqrt(hr) * ur, roots,
                     out=np.zeros_like(roots), where=roots > 0)
    croe = np.sqrt(0.5 * g * (hl + hr))
    sl = np.minimum.reduce((ul - cl, ur - cr, uroe - croe, np.zeros_like(hl)))
    sr = np.maximum.reduce((ul + cl, ur + cr, uroe + croe, np.zeros_like(hr)))
    fl, fr = _physical_flux(left, g), _physical_flux(right, g)
    numerator = sr[:, None] * fl - sl[:, None] * fr + (sl * sr)[:, None] * (right - left)
    denom = sr - sl
    flux = np.divide(numerator, denom[:, None], out=np.zeros_like(numerator),
                     where=denom[:, None] > 0)
    return flux, np.maximum(-sl, sr)


def interface_fluxes(state, bed, cfg):
    """Return flux used by left cell, right cell, and interface speed."""
    ext = np.vstack((state[-1] if cfg.boundary == "periodic" else state[0],
                     state,
                     state[0] if cfg.boundary == "periodic" else state[-1]))
    bext = np.concatenate(([bed[-1] if cfg.boundary == "periodic" else bed[0]],
                           bed,
                           [bed[0] if cfg.boundary == "periodic" else bed[-1]]))
    if cfg.boundary == "reflective":
        ext[0, 1] *= -1
        ext[-1, 1] *= -1
    left, right = ext[:-1], ext[1:]
    bstar = np.maximum(bext[:-1], bext[1:])
    hl = np.maximum(0.0, left[:, 0] + bext[:-1] - bstar)
    hr = np.maximum(0.0, right[:, 0] + bext[1:] - bstar)
    ul, ur = left[:, 1] / left[:, 0], right[:, 1] / right[:, 0]
    ls, rs = np.column_stack((hl, hl * ul)), np.column_stack((hr, hr * ur))
    base, speeds = hll_flux(ls, rs, cfg.g)
    for_left, for_right = base.copy(), base.copy()
    for_left[:, 1] += 0.5 * cfg.g * (left[:, 0] ** 2 - hl ** 2)
    for_right[:, 1] += 0.5 * cfg.g * (right[:, 0] ** 2 - hr ** 2)
    return for_left, for_right, speeds


def simulate(grid: Grid, bed, h0, q0, output_times,
             config: Optional[SolverConfig] = None) -> SimulationResult:
    """Integrate cell averages to prescribed output times.

    Output times must start at zero and be strictly increasing. Adaptive
    steps stop exactly at outputs. No clipping or mass correction is applied.
    An outflow boundary copies ghost cells: it is NOT exactly nonreflecting.
    Diagnostics integrate the actual numerical boundary mass flux.
    """
    cfg = config or SolverConfig()
    bed, h, q = (np.asarray(a, dtype=float).copy() for a in (bed, h0, q0))
    if any(a.shape != (grid.cells,) for a in (bed, h, q)):
        raise ValueError("bed, h0, q0 must each have shape (grid.cells,)")
    if not all(np.isfinite(a).all() for a in (bed, h, q)):
        raise ValueError("Initial arrays must be finite")
    if np.any(h <= cfg.minimum_depth):
        raise ValueError("V0.1 is fully wet: all depths must exceed minimum_depth")
    times = np.asarray(output_times, dtype=float)
    if (times.ndim != 1 or len(times) < 1 or not np.isfinite(times).all()
            or times[0] != 0 or np.any(np.diff(times) <= 0)):
        raise ValueError("output_times must start at zero and be strictly increasing")
    state = np.column_stack((h, q))
    states = [state.copy()]
    mass0 = float(h.sum() * grid.dx)
    boundary_integral, time, steps = 0.0, 0.0, 0
    minimum_h = float(h.min())
    max_cfl, max_residual, minimum_dt, maximum_dt = 0.0, 0.0, np.inf, 0.0
    masses, integrals, residuals = [mass0], [0.0], [0.0]
    for target in times[1:]:
        while time < target:
            if steps >= cfg.max_steps:
                raise RuntimeError("Maximum step count reached")
            fl, fr, interface_speeds = interface_fluxes(state, bed, cfg)
            cell_speed = np.max(np.abs(state[:, 1] / state[:, 0]) + np.sqrt(cfg.g * state[:, 0]))
            speed = max(float(cell_speed), float(interface_speeds.max()))
            dt = min(cfg.cfl * grid.dx / speed, float(target - time))
            if dt <= 0 or time + dt == time:
                raise RuntimeError("Time step underflow")
            updated = state - (dt / grid.dx) * (fl[1:] - fr[:-1])
            if not np.isfinite(updated).all() or np.any(updated[:, 0] <= cfg.minimum_depth):
                raise RuntimeError("Nonfinite or near-dry state: V0.1 fully-wet assumption failed; no clipping applied")
            boundary_integral += dt * float(fl[-1, 0] - fr[0, 0])
            state = updated
            time = float(target) if dt == target - time else time + dt
            steps += 1
            minimum_h = min(minimum_h, float(state[:, 0].min()))
            max_cfl = max(max_cfl, dt * speed / grid.dx)
            minimum_dt, maximum_dt = min(minimum_dt, dt), max(maximum_dt, dt)
            mass = float(state[:, 0].sum() * grid.dx)
            residual = mass - mass0 + boundary_integral
            max_residual = max(max_residual, abs(residual))
        states.append(state.copy())
        masses.append(float(state[:, 0].sum() * grid.dx))
        integrals.append(boundary_integral)
        residuals.append(masses[-1] - mass0 + boundary_integral)
    diagnostics = {
        "steps": steps, "minimum_depth_m": minimum_h, "maximum_actual_cfl": max_cfl,
        "minimum_dt_s": None if steps == 0 else minimum_dt,
        "maximum_dt_s": None if steps == 0 else maximum_dt,
        "initial_mass_per_unit_width_m2": mass0,
        "final_mass_per_unit_width_m2": masses[-1],
        "integrated_outward_boundary_mass_flux_m2": boundary_integral,
        "maximum_normalized_mass_balance_residual": max_residual / abs(mass0),
        "mass_per_output_m2": masses, "boundary_integral_per_output_m2": integrals,
        "mass_balance_residual_per_output_m2": residuals,
        "boundary": cfg.boundary, "depth_clipping_count": 0,
    }
    return SimulationResult(grid.x, bed, times.copy(), np.asarray(states), diagnostics)
