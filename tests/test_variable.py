import unittest
import numpy as np
from tsunami import Grid,SolverConfig
from tsunami.initial import cell_average,sine_simple_wave
from tsunami.flat import simulate_linear
from tsunami.variable import simulate_linear_variable,simulate_nonlinear_variable,linear_interface_flux
from tsunami.manufactured import setup,exact_state,analytic_forcing,fields,G
class VariableTests(unittest.TestCase):
    def test_lake_at_rest_smooth_and_step(self):
        grid=Grid(0,1000,160)
        for kind in ['smooth','step']:
            d=10+2*np.cos(2*np.pi*grid.x/1000) if kind=='smooth' else np.where(grid.x<500,10.,4.)
            for boundary in ['periodic','outflow','reflective']:
                cfg=SolverConfig(boundary=boundary,cfl=.4)
                mode='centered' if kind=='smooth' else 'constant'
                for model in ['linear','nonlinear']:
                    if model=='linear':r=simulate_linear_variable(grid,d,np.full(160,.3),np.zeros(160),[0,20],cfg,mode)
                    else:r=simulate_nonlinear_variable(grid,-d,d+.3,np.zeros(160),[0,20],cfg,mode)
                    self.assertLess(np.max(np.abs(r.eta-.3)),2e-12)
                    self.assertLess(np.max(np.abs(r.velocity)),2e-12)
    def test_smooth_variable_bed_manufactured_convergence(self):
        for model in ['linear','nonlinear']:
            errors=[]
            for n in [100,200,400]:
                grid=Grid(0,1000,n);d,z,f=setup(grid,model)
                fn=simulate_linear_variable if model=='linear' else simulate_nonlinear_variable
                r=fn(grid,d if model=='linear' else -d,z[:,0],z[:,1],[0,20],forcing=f)
                exact=cell_average(grid,lambda x:exact_state(x,20,model))[:,0]
                errors.append(float(np.mean(np.abs(r.states[-1,:,0]-exact))))
                self.assertLess(r.diagnostics['maximum_normalized_mass_balance_residual'],1e-12)
            self.assertGreater(np.log2(errors[-2]/errors[-1]),1.7)
    def test_forcing_independent_complex_step_residual(self):
        x=np.array([110.,280.,710.]);t=3.;delta=1e-25
        for model in ['linear','nonlinear']:
            wt=np.imag(exact_state(x,t+1j*delta,model))/delta
            d,e,u=fields(x+1j*delta,t)
            physical=np.zeros((3,2))
            if model=='linear':flux=np.stack((d*u,G*e),axis=-1)
            else:
                h=d+e;flux=np.stack((h*u,h*u*u+.5*G*h*h),axis=-1)
                dr=np.imag(d)/delta;dreal,ereal,_=fields(x,t)
                physical[:,1]=G*(dreal+ereal)*dr
            residual=wt+np.imag(flux)/delta-physical
            np.testing.assert_allclose(residual,analytic_forcing(x,t,model),rtol=1e-11,atol=1e-12)
    def test_variable_linear_reduces_to_flat_solver(self):
        grid=Grid(0,1000,100);eta=.2*np.sin(2*np.pi*grid.x/1000);u=np.sqrt(9.81/10)*eta
        cfg=SolverConfig(boundary='periodic',cfl=.4)
        a=simulate_linear(grid,10,eta,u,[0,10],cfg)
        b=simulate_linear_variable(grid,np.full(100,10.),eta,u,[0,10],cfg)
        np.testing.assert_allclose(a.states,b.states,atol=1e-13,rtol=1e-12)
    def test_step_interface_transmission_identity(self):
        dl,dr=10.,4.;cl,cr=np.sqrt(9.81*dl),np.sqrt(9.81*dr)
        left=np.array([[1.,cl/dl,-dl]]);right=np.array([[0.,0.,-dr]])
        flux,_=linear_interface_flux(left,right,9.81)
        transmission=2*cl/(cl+cr);reflection=(cl-cr)/(cl+cr)
        self.assertAlmostEqual(flux[0,1]/9.81,transmission)
        self.assertAlmostEqual(flux[0,0],cr*transmission)
        self.assertAlmostEqual(1+reflection,transmission)
        self.assertAlmostEqual(reflection**2+(cr/cl)*transmission**2,1.)
    def test_periodic_and_outflow_mass_cfl(self):
        grid=Grid(0,1000,120);d=10+2*np.cos(2*np.pi*grid.x/1000)
        eta=.15*np.exp(-((grid.x-700)/80)**2);u=np.sqrt(9.81/10)*eta
        for boundary in ['periodic','outflow','reflective']:
            cfg=SolverConfig(boundary=boundary,cfl=.4)
            for model in ['linear','nonlinear']:
                if model=='linear':r=simulate_linear_variable(grid,d,eta,u,[0,10,30],cfg)
                else:r=simulate_nonlinear_variable(grid,-d,d+eta,(d+eta)*u,[0,10,30],cfg)
                self.assertLess(r.diagnostics['maximum_normalized_mass_balance_residual'],1e-12)
                self.assertLessEqual(r.diagnostics['maximum_actual_cfl'],.4*(1+2e-12))
                self.assertEqual(r.diagnostics['depth_clipping_count'],0)
    def test_flat_nonlinear_exact_convergence(self):
        errors=[]
        for n in [100,200,400]:
            grid=Grid(0,1000,n)
            z=cell_average(grid,lambda x:sine_simple_wave(x,0,10,.2,1000))
            r=simulate_nonlinear_variable(grid,np.full(n,-10.),z[:,0],z[:,1],[0,20])
            exact=cell_average(grid,lambda x:sine_simple_wave(x,20,10,.2,1000))[:,0]-10
            errors.append(np.mean(np.abs(r.eta[-1]-exact)))
        self.assertGreater(np.log2(errors[-2]/errors[-1]),1.7)
    def test_invalid_inputs_and_no_mutation(self):
        grid=Grid(0,100,20);d=np.full(20,10.);eta=np.zeros(20);u=np.zeros(20)
        originals=[a.copy() for a in [d,eta,u]]
        simulate_linear_variable(grid,d,eta,u,[0,1])
        for a,b in zip([d,eta,u],originals):np.testing.assert_array_equal(a,b)
        for bad in [np.zeros(20),np.full(20,np.nan),np.ones(19)]:
            with self.assertRaises(ValueError):simulate_linear_variable(grid,bad,eta,u,[0,1])
        with self.assertRaises(ValueError):simulate_linear_variable(grid,d,eta,u,[0,1],bed_reconstruction='unknown')
        with self.assertRaises(ValueError):simulate_nonlinear_variable(grid,-d,np.zeros(20),u,[0,1])
        with self.assertRaises(ValueError):simulate_linear_variable(grid,d,eta,u,[1,2])
        with self.assertRaises(ValueError):simulate_linear_variable(grid,d,eta,u,[0,1],forcing=lambda t:np.zeros(20))
