"""Physics/analysis-based checks; run with python -m unittest discover."""
import unittest
import numpy as np
from tsunami import Grid, SolverConfig, simulate
from tsunami.initial import cell_average, right_going_state, sine_simple_wave, wet_dam_break


class SolverTests(unittest.TestCase):
    def test_lake_at_rest_smooth_and_step(self):
        grid = Grid(0, 1000, 200)
        for bed in (-10 + 3*np.sin(2*np.pi*grid.x/1000),
                    np.where(grid.x < 500, -10.0, -4.0)):
            run = simulate(grid, bed, -bed, np.zeros(grid.cells), [0, 100],
                           SolverConfig(boundary="reflective"))
            self.assertLess(np.max(np.abs(run.eta)), 1e-10)
            self.assertLess(np.max(np.abs(run.velocity)), 1e-10)

    def test_periodic_mass_conservation(self):
        grid = Grid(0, 1000, 200)
        initial = cell_average(grid, lambda x: sine_simple_wave(x, 0, 10, 0.2, 1000))
        run = simulate(grid, np.full(grid.cells, -10), initial[:, 0], initial[:, 1],
                       [0, 20], SolverConfig(boundary="periodic"))
        self.assertLess(run.diagnostics["maximum_normalized_mass_balance_residual"], 1e-12)
        self.assertEqual(run.diagnostics["integrated_outward_boundary_mass_flux_m2"], 0)

    def test_outflow_mass_flux_balance(self):
        grid = Grid(0, 1000, 200)
        initial = cell_average(grid, lambda x: right_going_state(
            10 + 0.1*np.exp(-0.5*((x-750)/50)**2), 10))
        run = simulate(grid, np.full(grid.cells, -10), initial[:, 0], initial[:, 1],
                       [0, 50], SolverConfig(boundary="outflow"))
        self.assertGreater(abs(run.diagnostics["integrated_outward_boundary_mass_flux_m2"]), 0.1)
        self.assertLess(run.diagnostics["maximum_normalized_mass_balance_residual"], 1e-12)

    def test_reflective_walls_conserve_mass(self):
        grid = Grid(0, 1000, 150)
        h = 10 + 0.2*np.exp(-((grid.x-800)/40)**2)
        run = simulate(grid, np.full(grid.cells, -10), h, np.zeros(grid.cells),
                       [0, 50], SolverConfig(boundary="reflective"))
        self.assertLess(run.diagnostics["maximum_normalized_mass_balance_residual"], 1e-12)
        self.assertEqual(run.diagnostics["integrated_outward_boundary_mass_flux_m2"], 0)

    def test_uniform_moving_water(self):
        grid = Grid(0, 100, 50)
        run = simulate(grid, np.full(50, -10), np.full(50, 10), np.full(50, 20), [0, 10])
        np.testing.assert_array_equal(run.states[-1], run.states[0])

    def test_nonlinear_exact_wave_convergence(self):
        errors = []
        for n in [100, 200, 400]:
            grid = Grid(0, 1000, n)
            initial = cell_average(grid, lambda x: sine_simple_wave(x, 0, 10, 0.2, 1000))
            reference = cell_average(grid, lambda x: sine_simple_wave(x, 20, 10, 0.2, 1000))
            run = simulate(grid, np.full(n, -10), initial[:, 0], initial[:, 1], [0, 20],
                           SolverConfig(boundary="periodic"))
            errors.append(float(np.mean(np.abs(run.states[-1, :, 0]-reference[:, 0]))))
        orders = np.log2(np.array(errors[:-1])/errors[1:])
        self.assertTrue(np.all(orders > 0.8), (errors, orders))
        self.assertTrue(np.all(orders < 1.2), (errors, orders))

    def test_small_wave_speed(self):
        grid = Grid(0, 1000, 400)
        c = np.sqrt(9.81*10)
        amp = 1e-4
        init = cell_average(grid, lambda x: np.stack((10+amp*np.sin(2*np.pi*x/1000),
                                                     c*amp*np.sin(2*np.pi*x/1000)), axis=-1))
        time = 0.25*1000/c
        run = simulate(grid, np.full(grid.cells, -10), init[:, 0], init[:, 1], [0, time],
                       SolverConfig(boundary="periodic"))
        mode = np.exp(-2j*np.pi*grid.x/1000)
        initial_mode = np.sum((init[:, 0]-10)*mode)
        final_mode = np.sum(run.eta[-1]*mode)
        speed = -np.angle(final_mode/initial_mode)*1000/(2*np.pi*time)
        self.assertLess(abs(speed-c)/c, 0.005)

    def test_fully_wet_and_cfl(self):
        grid = Grid(0, 100, 100)
        h = 10+0.2*np.cos(2*np.pi*grid.x/100)
        run = simulate(grid, np.full(100, -10), h, np.zeros(100), [0, 5],
                       SolverConfig(boundary="periodic", cfl=0.45))
        self.assertGreater(run.diagnostics["minimum_depth_m"], 0)
        self.assertLessEqual(run.diagnostics["maximum_actual_cfl"], 0.45+1e-14)
        self.assertEqual(run.diagnostics["depth_clipping_count"], 0)

    def test_input_validation_and_no_mutation(self):
        grid = Grid(0, 100, 10)
        bed, h, q = np.full(10, -10.0), np.full(10, 10.0), np.zeros(10)
        simulate(grid, bed, h, q, [0, 1])
        np.testing.assert_array_equal(h, np.full(10, 10.0))
        np.testing.assert_array_equal(q, np.zeros(10))
        with self.assertRaises(ValueError):
            simulate(grid, bed, np.zeros(10), q, [0, 1])
        with self.assertRaises(ValueError):
            simulate(grid, bed, h, q, [0, 1, 1])
        with self.assertRaises(ValueError):
            SolverConfig(cfl=0.8)

    def test_zero_duration(self):
        grid = Grid(0, 100, 10)
        run = simulate(grid, np.full(10, -10), np.full(10, 10), np.zeros(10), [0])
        self.assertEqual(run.states.shape, (1, 10, 2))
        self.assertEqual(run.diagnostics["steps"], 0)
        self.assertEqual(run.diagnostics["maximum_normalized_mass_balance_residual"], 0)

    def test_wet_dam_break_entropy_solution(self):
        errors = []
        for n in [400, 800]:
            grid = Grid(-100, 100, n)
            initial = wet_dam_break(grid.x, 0)
            reference = cell_average(grid, lambda x: wet_dam_break(x, 8), order=16)
            run = simulate(grid, np.zeros(n), initial[:, 0], initial[:, 1], [0, 8])
            errors.append(float(np.mean(np.abs(run.states[-1, :, 0]-reference[:, 0]))))
        self.assertLess(errors[-1], errors[0])
        self.assertLess(errors[-1], 0.02)


if __name__ == "__main__":
    unittest.main()
