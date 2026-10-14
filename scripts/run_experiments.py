"""Reproduce V0.1 verification metrics, figures and animation."""
from pathlib import Path
import sys
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import argparse
import csv
import hashlib
import json
import platform
import time
from datetime import datetime, timezone
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
from tsunami import Grid, SolverConfig, simulate, __version__
from tsunami.initial import cell_average, right_going_state, sine_simple_wave


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False), encoding="utf-8")


def lake_case(figures):
    grid = Grid(0, 1000, 300)
    bed = -10+3*np.sin(2*np.pi*grid.x/1000) + np.where(grid.x>650, 1.5, 0)
    run = simulate(grid, bed, -bed, np.zeros(grid.cells), [0, 10, 50, 100],
                   SolverConfig(boundary="reflective"))
    fig, axes = plt.subplots(2, 1, figsize=(9, 6), constrained_layout=True)
    axes[0].fill_between(grid.x, bed, 0, alpha=0.25, color="#168aad")
    axes[0].plot(grid.x, bed, color="#6c584c", label="Bed elevation")
    axes[0].plot(grid.x, run.eta[-1], color="#005f73", label="Surface after 100 s")
    axes[0].set(ylabel="Elevation (m)", title="Lake at rest over smooth and stepped bed")
    axes[0].legend()
    axes[1].plot(grid.x, run.eta[-1], label="Surface drift (m)", color="#005f73")
    axes[1].set(xlabel="Distance (m)", ylabel="Surface drift (m)")
    axes[1].ticklabel_format(axis="y", style="sci", scilimits=(0,0))
    axes[1].grid(alpha=0.2)
    fig.savefig(figures/"lake_at_rest.png", dpi=180)
    plt.close(fig)
    return {"max_surface_drift_m": float(np.max(np.abs(run.eta))),
            "max_velocity_m_s": float(np.max(np.abs(run.velocity))),
            "configuration": {"cells":300,"domain_m":[0,1000],"output_times_s":[0,10,50,100],"g":9.81,"cfl":0.45,"boundary":"reflective"},
            "diagnostics": run.diagnostics}


def convergence_case(figures, results):
    rows, runs = [], []
    for n in [100, 200, 400, 800]:
        grid = Grid(0, 1000, n)
        start = cell_average(grid, lambda x: sine_simple_wave(x, 0, 10, 0.2, 1000))
        exact = cell_average(grid, lambda x: sine_simple_wave(x, 20, 10, 0.2, 1000))
        run = simulate(grid, np.full(n,-10), start[:,0], start[:,1], [0,20],
                       SolverConfig(boundary="periodic"))
        error = run.states[-1,:,0]-exact[:,0]
        rows.append({"cells":n,"dx_m":grid.dx,"depth_L1_mean_error_m":float(np.mean(np.abs(error))),
                     "depth_Linf_error_m":float(np.max(np.abs(error))),
                     "q_L1_mean_error_m2_s":float(np.mean(np.abs(run.states[-1,:,1]-exact[:,1]))),
                     "mass_balance_residual":run.diagnostics["maximum_normalized_mass_balance_residual"],
                     "observed_L1_order":None})
        runs.append((grid,run,exact))
    for i in range(1,len(rows)):
        rows[i]["observed_L1_order"] = float(np.log2(rows[i-1]["depth_L1_mean_error_m"]/rows[i]["depth_L1_mean_error_m"]))
    with (results/"convergence.csv").open("w", newline="") as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    fig,axes=plt.subplots(1,2,figsize=(11,4.5),constrained_layout=True)
    xfine=np.linspace(0,1000,2000)
    axes[0].plot(xfine,sine_simple_wave(xfine,20,10,0.2,1000)[:,0]-10,"k--",label="Exact characteristic solution")
    for grid,run,exact in runs[::2]:
        axes[0].plot(grid.x,run.eta[-1],label=f"Finite volumes, N={grid.cells}")
    axes[0].set(xlabel="Distance (m)",ylabel="Surface elevation (m)",title="Smooth nonlinear wave, t = 20 s")
    axes[0].legend(fontsize=8)
    dx=np.array([r["dx_m"] for r in rows]);err=np.array([r["depth_L1_mean_error_m"] for r in rows])
    axes[1].loglog(dx,err,"o-",color="#005f73",label="Measured mean L1 error")
    axes[1].loglog(dx,err[-1]*dx/dx[-1],"k--",label="First-order slope")
    axes[1].set(xlabel="Cell width (m)",ylabel="Mean absolute depth error (m)",title="Grid convergence")
    axes[1].legend(fontsize=8);axes[1].grid(alpha=.2,which="both")
    fig.savefig(figures/"convergence.png",dpi=180);plt.close(fig)
    return {"configuration":{"depth_m":10,"amplitude_m":0.2,"domain_length_m":1000,"final_time_s":20,"cfl":0.45,"boundary":"periodic","quadrature_order":8},"rows":rows}


def speed_case():
    grid=Grid(0,1000,400);D=10;A=1e-4;c=np.sqrt(9.81*D);t=.25*1000/c
    start=cell_average(grid,lambda x: np.stack((D+A*np.sin(2*np.pi*x/1000),c*A*np.sin(2*np.pi*x/1000)),axis=-1))
    run=simulate(grid,np.full(grid.cells,-D),start[:,0],start[:,1],[0,t],SolverConfig(boundary="periodic"))
    mode=np.exp(-2j*np.pi*grid.x/1000)
    z0=np.sum((start[:,0]-D)*mode);z1=np.sum(run.eta[-1]*mode)
    numerical_speed=float(-np.angle(z1/z0)*1000/(2*np.pi*t))
    return {"reference_speed_m_s":float(c),"measured_phase_speed_m_s":numerical_speed,
            "relative_speed_error":float(abs(numerical_speed-c)/c),
            "amplitude_ratio":float(abs(z1/z0)),
            "note":"Nonlinear solver versus a linear small-amplitude wave; not an exact nonlinear solution.",
            "configuration":{"cells":400,"depth_m":10,"amplitude_m":A,"length_m":1000,"final_time_s":float(t),"cfl":0.45,"boundary":"periodic"}}


def pulse_case(figures, results, animate):
    grid=Grid(0,20000,800);D=50;A=1;x0=4000;sigma=800;times=np.linspace(0,300,61)
    start=cell_average(grid,lambda x:right_going_state(D+A*np.exp(-.5*((x-x0)/sigma)**2),D))
    run=simulate(grid,np.full(grid.cells,-D),start[:,0],start[:,1],times,SolverConfig(boundary="outflow"))
    np.savez_compressed(results/"pulse_snapshots.npz",x_m=run.x,time_s=run.times,bed_m=run.bed,h_m=run.states[:,:,0],q_m2_s=run.states[:,:,1])
    fig,ax=plt.subplots(figsize=(10,4.5),constrained_layout=True)
    for index in [0,20,40,60]:ax.plot(grid.x/1000,run.eta[index],label=f"t = {times[index]:.0f} s")
    ax.set(xlabel="Distance (km)",ylabel="Surface elevation (m)",title="Right-going pulse over a flat bed | first-order solver")
    ax.legend();ax.grid(alpha=.2);fig.savefig(figures/"pulse_propagation.png",dpi=180);plt.close(fig)
    if animate:
        fig,ax=plt.subplots(figsize=(9,4),constrained_layout=True)
        line,=ax.plot(grid.x/1000,run.eta[0],color="#005f73",lw=2)
        label=ax.text(.02,.92,"",transform=ax.transAxes)
        ax.set(xlim=(0,20),ylim=(-.08,1.1),xlabel="Distance (km)",ylabel="Surface elevation (m)",title="Computed fully-wet long-wave pulse (not event data)")
        ax.grid(alpha=.2)
        def update(i):line.set_ydata(run.eta[i]);label.set_text(f"t = {times[i]:.0f} s");return line,label
        anim=FuncAnimation(fig,update,frames=len(times),interval=80,blit=False)
        anim.save(figures/"pulse_propagation.gif",writer=PillowWriter(fps=12),dpi=100);plt.close(fig)
    return {"configuration":{"cells":800,"domain_m":[0,20000],"depth_m":D,"amplitude_m":A,"gaussian_center_m":x0,"gaussian_sigma_m":sigma,"final_time_s":300,"output_interval_s":5,"boundary":"outflow","cfl":0.45},
            "initial_peak_m":float(run.eta[0].max()),"final_peak_m":float(run.eta[-1].max()),
            "note":"Peak reduction includes first-order numerical diffusion; it is not evidence of physical dissipation.","diagnostics":run.diagnostics}


def outflow_case():
    grid=Grid(0,1000,200)
    start=cell_average(grid,lambda x:right_going_state(10+.1*np.exp(-.5*((x-750)/50)**2),10))
    run=simulate(grid,np.full(200,-10),start[:,0],start[:,1],[0,50],SolverConfig(boundary="outflow"))
    return {"configuration":{"cells":200,"domain_m":[0,1000],"depth_m":10,"amplitude_m":.1,"pulse_center_m":750,"pulse_sigma_m":50,"final_time_s":50,"cfl":.45,"boundary":"outflow"},"diagnostics":run.diagnostics}


def main():
    parser=argparse.ArgumentParser();parser.add_argument("--no-animation",action="store_true");args=parser.parse_args()
    figures,results=ROOT/"figures",ROOT/"results"
    figures.mkdir(exist_ok=True);results.mkdir(exist_ok=True)
    plt.rcParams.update({"font.family":"DejaVu Sans","axes.spines.top":False,"axes.spines.right":False})
    start=time.perf_counter()
    metrics={"lake_at_rest":lake_case(figures),"nonlinear_convergence":convergence_case(figures,results),
             "small_amplitude_speed":speed_case(),"outflow_mass_balance":outflow_case(),
             "pulse_demo":pulse_case(figures,results,not args.no_animation)}
    write_json(results/"verification_metrics.json",metrics)
    source_hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()
                   for pattern in ["src/**/*.py","scripts/*.py","tests/*.py"] for p in ROOT.glob(pattern)}
    write_json(results/"run_manifest.json",{"created_utc":datetime.now(timezone.utc).isoformat(),"version":__version__,
        "python":platform.python_version(),"platform":platform.platform(),"numpy":np.__version__,"matplotlib":matplotlib.__version__,
        "source_sha256":source_hashes,"elapsed_seconds":time.perf_counter()-start,
        "command":"python scripts/run_experiments.py"+(" --no-animation" if args.no_animation else ""),
        "scope":"Numerical verification only; no laboratory or field validation"})
    print(json.dumps({"lake_drift_m":metrics["lake_at_rest"]["max_surface_drift_m"],
        "convergence_orders":[r["observed_L1_order"] for r in metrics["nonlinear_convergence"]["rows"]],
        "relative_speed_error":metrics["small_amplitude_speed"]["relative_speed_error"],
        "outflow_mass_balance_residual":metrics["outflow_mass_balance"]["diagnostics"]["maximum_normalized_mass_balance_residual"]},indent=2))


if __name__=="__main__":main()
