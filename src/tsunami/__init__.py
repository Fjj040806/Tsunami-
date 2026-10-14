"""Verified fully wet shallow-water modeling for research education."""

from .solver import Grid, SolverConfig, SimulationResult, simulate

__version__ = "1.0.0"
from .flat import FlatResult, simulate_linear, simulate_nonlinear_flat
__all__ = ["Grid", "SolverConfig", "SimulationResult", "simulate",
           "FlatResult", "simulate_linear", "simulate_nonlinear_flat"]

from .variable import VariableResult, simulate_linear_variable, simulate_nonlinear_variable
__all__ += ["VariableResult", "simulate_linear_variable", "simulate_nonlinear_variable"]

from .variable import IncidentBoundary
__all__ += ["IncidentBoundary"]
