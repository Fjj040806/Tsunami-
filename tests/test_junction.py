import unittest
import numpy as np
from tsunami.variable import stationary_step_flux,simulate_nonlinear_variable
from tsunami import Grid,SolverConfig
from tsunami.initial import cell_average
from tsunami.studies import step_reference,gaussian
class JunctionTests(unittest.TestCase):
    def test_rest_pressure_and_common_mass(self):
        l=np.array([[.3,0.,-20.]]);r=np.array([[.3,0.,-5.]])
        fl,fr,_=stationary_step_flux(l,r,9.81)
        self.assertLess(abs(fl[0,0]),1e-12)
        self.assertEqual(fl[0,0],fr[0,0])
        self.assertAlmostEqual(fl[0,1],.5*9.81*20.3**2)
        self.assertAlmostEqual(fr[0,1],.5*9.81*5.3**2)
    def test_low_amplitude_junction_linear_limit(self):
        A=1e-4;cl=np.sqrt(9.81*20);cr=np.sqrt(9.81*5)
        l=np.array([[A,cl/20*A,-20.]]);r=np.array([[0.,0.,-5.]])
        fl,fr,_=stationary_step_flux(l,r,9.81)
        self.assertLess(abs(fl[0,0]/(cr*(4/3)*A)-1),1e-4)
        self.assertAlmostEqual(fl[0,0],fr[0,0])
    def test_complete_low_amplitude_step_wave(self):
        grid=Grid(-1000,1000,800);D=np.where(grid.x<0,20.,5.);A=.002;c=np.sqrt(9.81*20)
        def init(x):
            d=np.where(x<0,20.,5.);e=gaussian(x,A,-300,40);u=c/20*e
            return np.stack((d+e,(d+e)*u),axis=-1)
        z=cell_average(grid,init)
        r=simulate_nonlinear_variable(grid,-D,z[:,0],z[:,1],[0,60],SolverConfig(boundary='outflow',cfl=.4),bed_reconstruction='constant',step_coupling='bernoulli')
        exact=cell_average(grid,lambda x:step_reference(x,60,20,5,A,-300,40))[:,0]
        self.assertLess(np.linalg.norm(r.eta[-1]-exact)/np.linalg.norm(exact),.04)
        self.assertLess(r.diagnostics['maximum_normalized_mass_balance_residual'],1e-12)
    def test_reject_supercritical_junction(self):
        with self.assertRaises(ValueError):stationary_step_flux(np.array([[0.,50.,-20.]]),np.array([[0.,0.,-5.]]),9.81)
