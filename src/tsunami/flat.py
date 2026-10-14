"""Flat-bed linear/nonlinear models with independently defined fluxes.

The models share reconstruction and RK2 infrastructure, not their equations.
V0.1's well-balanced variable-bed solver remains in solver.py unchanged.
"""
from dataclasses import dataclass
import numpy as np
from .solver import SolverConfig, hll_flux


@dataclass
class FlatResult:
    x: np.ndarray
    times: np.ndarray
    states: np.ndarray
    depth: float
    model: str
    diagnostics: dict

    @property
    def eta(self):
        return self.states[:, :, 0] if self.model == "linear" else self.states[:, :, 0]-self.depth

    @property
    def velocity(self):
        return self.states[:, :, 1] if self.model == "linear" else self.states[:, :, 1]/self.states[:, :, 0]


def _interfaces(state, boundary, order):
    if boundary == "periodic":
        ext = np.concatenate((state[-2:], state, state[:2]))
    elif boundary == "outflow":
        ext = np.pad(state, ((2, 2), (0, 0)), mode="edge")
    else:
        ext = np.concatenate((state[:2][::-1], state, state[-2:][::-1]))
        ext[:2, 1] *= -1
        ext[-2:, 1] *= -1
    center = ext[1:-1]
    if order == 1:
        slopes = np.zeros_like(center)
    else:
        dm, dp = center-ext[:-2], ext[2:]-center
        middle = 0.5*(dm+dp)
        # Monotonized-central limiter in each model's evolved variables.
        same = (dm*dp > 0)
        slopes = np.where(same, np.sign(middle)*np.minimum.reduce(
            (2*np.abs(dm), np.abs(middle), 2*np.abs(dp))), 0.0)
    return center[:-1]+0.5*slopes[:-1], center[1:]-0.5*slopes[1:]


def linear_flux(left, right, depth, g):
    """Exact upwind flux for W=(eta,u), A=[[0,D],[g,0]].

    Since |A|=sqrt(gD) I, Fhat=.5 A(WL+WR)-.5 c(WR-WL).
    This never calls the nonlinear HLL flux or nonlinear depth variables.
    """
    mean = 0.5*(left+right)
    return np.column_stack((depth*mean[:, 1], g*mean[:, 0])) - 0.5*np.sqrt(g*depth)*(right-left)


def _integrate(grid, initial, depth, output_times, model, cfg, order):
    if (np.ndim(depth) != 0 or not np.isfinite(depth) or depth <= 0):
        raise ValueError("V0.2 requires one positive scalar constant depth")
    if isinstance(order, bool) or order not in (1, 2):
        raise ValueError("order must be 1 or 2")
    state = np.asarray(initial, dtype=float).copy()
    if state.shape != (grid.cells, 2) or not np.isfinite(state).all():
        raise ValueError("Invalid initial state shape or values")
    times = np.asarray(output_times, dtype=float)
    if (times.ndim != 1 or len(times)<1 or not np.isfinite(times).all()
            or times[0] != 0 or np.any(np.diff(times)<=0)):
        raise ValueError("Output times must start at zero and increase strictly")
    def check(w):
        h = depth+w[:, 0] if model == "linear" else w[:, 0]
        if not np.isfinite(w).all() or np.any(h <= cfg.minimum_depth):
            raise ValueError("Fully-wet state violated; no depth clipping is applied")
    check(state)
    def fluxes(w):
        check(w)
        left, right = _interfaces(w, cfg.boundary, order)
        check(left); check(right)
        if model == "linear":
            f = linear_flux(left, right, depth, cfg.g)
            speed = float(np.sqrt(cfg.g*depth))
        else:
            f, speeds = hll_flux(left, right, cfg.g)
            speed = max(float(speeds.max()), float(np.max(np.abs(w[:,1]/w[:,0])+np.sqrt(cfg.g*w[:,0]))))
        return f, speed
    def mass(w):
        return float(np.sum((depth+w[:,0]) if model=="linear" else w[:,0])*grid.dx)
    mass0 = mass(state)
    integral, t, steps, rejected = 0.0, 0.0, 0, 0
    states=[state.copy()]; residuals=[0.0]; integrals=[0.0]
    max_residual, max_cfl = 0.0, 0.0
    min_h=float(np.min(depth+state[:,0] if model=="linear" else state[:,0]))
    for target in times[1:]:
        while t < target:
            if steps >= cfg.max_steps:
                raise RuntimeError("Maximum step count exceeded")
            f0,speed0=fluxes(state)
            dt=min(cfg.cfl*grid.dx/speed0,float(target-t))
            for attempt in range(20):
                if dt<=0 or t+dt==t:
                    raise RuntimeError("Time-step underflow")
                stage=state-dt/grid.dx*(f0[1:]-f0[:-1])
                if order==1:
                    updated=stage; net=float(f0[-1,0]-f0[0,0]); speed1=speed0
                    check(updated)
                else:
                    f1,speed1=fluxes(stage)
                    if dt*speed1/grid.dx > cfg.cfl*(1+1e-12):
                        dt=min(0.9*dt,cfg.cfl*grid.dx/speed1)
                        rejected+=1;continue
                    updated=.5*state+.5*(stage-dt/grid.dx*(f1[1:]-f1[:-1]))
                    check(updated)
                    net=.5*float(f0[-1,0]-f0[0,0]+f1[-1,0]-f1[0,0])
                break
            else:
                raise RuntimeError("RK-stage CFL adjustment failed")
            state=updated
            integral+=dt*net
            t=float(target) if dt==target-t else t+dt
            steps+=1
            max_cfl=max(max_cfl,dt*max(speed0,speed1)/grid.dx)
            min_h=min(min_h,float(np.min(depth+state[:,0] if model=="linear" else state[:,0])))
            max_residual=max(max_residual,abs(mass(state)-mass0+integral))
        states.append(state.copy());residuals.append(mass(state)-mass0+integral);integrals.append(integral)
    diagnostics={"model":model,"order":order,"boundary":cfg.boundary,"g":cfg.g,"cfl":cfg.cfl,
        "steps":steps,"rejected_cfl_trials":rejected,"max_actual_stage_cfl":max_cfl,
        "minimum_depth_m":min_h,"depth_clipping_count":0,
        "maximum_normalized_mass_balance_residual":max_residual/abs(mass0),
        "boundary_integral_per_output_m2":integrals,"mass_balance_residual_per_output_m2":residuals}
    return FlatResult(grid.x,times.copy(),np.asarray(states),float(depth),model,diagnostics)


def simulate_linear(grid, depth, eta0, u0, output_times, config=None, order=2):
    """Independent constant-depth linear model, cell averages of (eta,u)."""
    cfg=config or SolverConfig(boundary="periodic",cfl=0.4)
    return _integrate(grid,np.column_stack((eta0,u0)),depth,output_times,"linear",cfg,order)


def simulate_nonlinear_flat(grid, depth, h0, q0, output_times, config=None, order=2):
    """Flat-bed nonlinear model, cell averages of (h,q). No bed source."""
    cfg=config or SolverConfig(boundary="periodic",cfl=0.4)
    return _integrate(grid,np.column_stack((h0,q0)),depth,output_times,"nonlinear",cfg,order)
