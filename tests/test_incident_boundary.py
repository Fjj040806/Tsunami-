import unittest
import numpy as np
from tsunami import Grid,SolverConfig,IncidentBoundary
from tsunami.variable import simulate_linear_variable,simulate_nonlinear_variable
from tsunami.initial import cell_average
class IncidentTests(unittest.TestCase):
    def test_linear_incoming_pulse_and_mass(self):
        grid=Grid(0,100,200);D=np.full(200,10.);c=np.sqrt(98.1)
        f=lambda t:.02*np.exp(-((t-2)/.5)**2)
        r=simulate_linear_variable(grid,D,np.zeros(200),np.zeros(200),[0,5],SolverConfig(boundary='outflow',cfl=.4),boundary=IncidentBoundary(f))
        exact=cell_average(grid,lambda x:f(5-x/c))
        self.assertLess(np.mean(np.abs(r.eta[-1]-exact)),1e-4)
        self.assertLess(r.diagnostics['maximum_normalized_mass_balance_residual'],1e-12)
        self.assertTrue(r.diagnostics['incident_boundary_enabled'])
    def test_linear_right_wall_reflection(self):
        grid=Grid(0,100,300);D=np.full(300,10.);c=np.sqrt(98.1)
        f=lambda t:.02*np.exp(-((t-2)/.5)**2)
        r=simulate_linear_variable(grid,D,np.zeros(300),np.zeros(300),[0,15],SolverConfig(boundary='outflow',cfl=.4),boundary=IncidentBoundary(f))
        exact=cell_average(grid,lambda x:f(15-x/c)+f(15-(200-x)/c))
        self.assertLess(np.mean(np.abs(r.eta[-1]-exact)),1e-4)
        self.assertLess(r.diagnostics['maximum_normalized_mass_balance_residual'],1e-12)
    def test_nonlinear_zero_input_variable_bed_balance(self):
        grid=Grid(0,100,100);D=10+2*np.cos(2*np.pi*grid.x/100)
        r=simulate_nonlinear_variable(grid,-D,D,np.zeros(100),[0,10],SolverConfig(boundary='outflow',cfl=.4),boundary=IncidentBoundary(lambda t:0.))
        self.assertLess(np.max(np.abs(r.eta)),2e-12)
        self.assertLess(r.diagnostics['maximum_normalized_mass_balance_residual'],1e-12)
    def test_boundary_contract(self):
        grid=Grid(0,100,20);D=np.full(20,10.);z=np.zeros(20)
        with self.assertRaises(ValueError):simulate_linear_variable(grid,D,z,z,[0,1],boundary=IncidentBoundary(lambda t:0.))
        with self.assertRaises(ValueError):simulate_linear_variable(grid,D,z,z,[0,1],SolverConfig(boundary='outflow'),boundary=IncidentBoundary(lambda t:np.nan))
