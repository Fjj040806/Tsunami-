"""A small experiment you can edit without installing the package."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np
import matplotlib.pyplot as plt
from tsunami import Grid, SolverConfig, simulate
from tsunami.initial import cell_average, right_going_state

grid = Grid(0, 20000, 800)  # metres, 800 finite-volume cells
depth = 50.0
amplitude = 1.0
initial = cell_average(grid, lambda x: right_going_state(
    depth + amplitude * np.exp(-0.5*((x-4000)/800)**2), depth))
run = simulate(grid, np.full(grid.cells, -depth), initial[:, 0], initial[:, 1],
               [0, 100, 200, 300], SolverConfig(boundary="outflow"))
for t, eta in zip(run.times, run.eta):
    plt.plot(run.x/1000, eta, label=f"{t:.0f} s")
plt.xlabel("Distance (km)")
plt.ylabel("Surface elevation (m)")
plt.title("Fully-wet long-wave propagation, V0.1")
plt.legend()
plt.show()
print(run.diagnostics)
