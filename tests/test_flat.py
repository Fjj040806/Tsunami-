import unittest
import numpy as np
from tsunami import Grid, SolverConfig
from tsunami.flat import simulate_linear, simulate_nonlinear_flat
from tsunami.initial import cell_average, sine_simple_wave


class FlatTests(unittest.TestCase):
    def test_linear_exact_both_directions(self):
        grid=Grid(0,1000,400);D=10;g=9.81;c=np.sqrt(g*D);t=20
        def exact(x,t):
            right=.1*np.sin(2*np.pi*(x-c*t)/1000)
            left=.04*np.cos(4*np.pi*(x+c*t)/1000)
            return np.stack((right+left,c/D*(right-left)),axis=-1)
        initial=cell_average(grid,lambda x:exact(x,0))
        reference=cell_average(grid,lambda x:exact(x,t))
        run=simulate_linear(grid,D,initial[:,0],initial[:,1],[0,t])
        self.assertLess(np.mean(np.abs(run.states[-1,:,0]-reference[:,0]))/.1,0.001)
        self.assertLess(np.mean(np.abs(run.states[-1,:,1]-reference[:,1]))/(.1*c/D),0.001)

    def test_linear_second_order_convergence(self):
        errors=[];D=10;c=np.sqrt(9.81*D)
        for n in [100,200,400]:
            grid=Grid(0,1000,n)
            initial=cell_average(grid,lambda x:np.stack((.1*np.sin(2*np.pi*x/1000),
                c/D*.1*np.sin(2*np.pi*x/1000)),axis=-1))
            reference=cell_average(grid,lambda x:.1*np.sin(2*np.pi*(x-c*20)/1000))
            run=simulate_linear(grid,D,initial[:,0],initial[:,1],[0,20])
            errors.append(np.mean(np.abs(run.eta[-1]-reference)))
        self.assertTrue(np.all(np.log2(np.array(errors[:-1])/errors[1:])>1.7),errors)

    def test_nonlinear_second_order_convergence(self):
        errors=[]
        for n in [100,200,400]:
            grid=Grid(0,1000,n)
            initial=cell_average(grid,lambda x:sine_simple_wave(x,0,10,.2,1000))
            reference=cell_average(grid,lambda x:sine_simple_wave(x,20,10,.2,1000))
            run=simulate_nonlinear_flat(grid,10,initial[:,0],initial[:,1],[0,20])
            errors.append(np.mean(np.abs(run.states[-1,:,0]-reference[:,0])))
        self.assertTrue(np.all(np.log2(np.array(errors[:-1])/errors[1:])>1.7),errors)

    def test_conservation_cfl_and_positive_water(self):
        grid=Grid(0,1000,200);D=10;c=np.sqrt(9.81*D)
        eta=.5*np.sin(2*np.pi*grid.x/1000);u=c/D*eta
        for run in [simulate_linear(grid,D,eta,u,[0,20]),
                    simulate_nonlinear_flat(grid,D,D+eta,(D+eta)*u,[0,20])]:
            self.assertLess(run.diagnostics['maximum_normalized_mass_balance_residual'],1e-12)
            self.assertLessEqual(run.diagnostics['max_actual_stage_cfl'],.4+1e-12)
            self.assertGreater(run.diagnostics['minimum_depth_m'],0)

    def test_linear_amplitude_scaling(self):
        grid=Grid(0,1000,100);eta=.1*np.sin(2*np.pi*grid.x/1000)
        a=simulate_linear(grid,10,eta,eta,[0,10])
        b=simulate_linear(grid,10,3*eta,3*eta,[0,10])
        np.testing.assert_allclose(b.states,3*a.states,atol=1e-12)

    def test_uniform_nonlinear_state_and_reflective_mass(self):
        grid=Grid(0,1000,100)
        uniform=simulate_nonlinear_flat(grid,10,np.full(100,10),np.full(100,20),[0,10])
        np.testing.assert_allclose(uniform.states[-1],uniform.states[0],atol=1e-13)
        eta=.1*np.exp(-((grid.x-700)/100)**2)
        for model in ['linear','nonlinear']:
            cfg=SolverConfig(boundary='reflective',cfl=.4)
            run=(simulate_linear(grid,10,eta,np.zeros(100),[0,50],cfg) if model=='linear'
                 else simulate_nonlinear_flat(grid,10,10+eta,np.zeros(100),[0,50],cfg))
            self.assertLess(run.diagnostics['maximum_normalized_mass_balance_residual'],1e-12)
            self.assertEqual(run.diagnostics['boundary_integral_per_output_m2'][-1],0)

    def test_invalid_depth_and_near_dry_input(self):
        grid=Grid(0,100,10)
        with self.assertRaises(ValueError):simulate_linear(grid,np.ones(10),np.zeros(10),np.zeros(10),[0,1])
        with self.assertRaises(ValueError):simulate_linear(grid,10,np.full(10,-10),np.zeros(10),[0,1])
        with self.assertRaises(ValueError):simulate_nonlinear_flat(grid,10,np.ones(10),np.zeros(10),[0,1],order=3)

    def test_rk2_open_boundary_flux_balance(self):
        grid=Grid(0,1000,200);D=10;c=np.sqrt(9.81*D)
        eta=.1*np.exp(-.5*((grid.x-750)/50)**2);u=c/D*eta
        cfg=SolverConfig(boundary='outflow',cfl=.4)
        for run in [simulate_linear(grid,D,eta,u,[0,50],cfg),
                    simulate_nonlinear_flat(grid,D,D+eta,(D+eta)*u,[0,50],cfg)]:
            self.assertGreater(abs(run.diagnostics['boundary_integral_per_output_m2'][-1]),.1)
            self.assertLess(run.diagnostics['maximum_normalized_mass_balance_residual'],1e-12)


if __name__=='__main__':unittest.main()
