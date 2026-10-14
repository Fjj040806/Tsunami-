"""Run independent verification and amplitude/model comparisons for V0.2."""
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
import csv
import json
import hashlib
import platform
import time
from datetime import datetime,timezone
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation,PillowWriter
from tsunami import Grid,SolverConfig,__version__
from tsunami.flat import simulate_linear,simulate_nonlinear_flat
from tsunami.initial import cell_average,sine_simple_wave


def rms(a):return float(np.sqrt(np.mean(np.asarray(a)**2)))


def write_json(path,value):path.write_text(json.dumps(value,indent=2,allow_nan=False),encoding='utf-8')


def write_csv(path,rows):
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)


def linear_initial(grid,A,D,L,g):
    c=np.sqrt(g*D)
    return cell_average(grid,lambda x:np.stack((A*np.sin(2*np.pi*x/L),c/D*A*np.sin(2*np.pi*x/L)),axis=-1))


def matched_initial(grid,A,D,L,g):
    c=np.sqrt(g*D)
    def f(x):
        eta=A*np.sin(2*np.pi*x/L);u=c/D*eta
        return np.stack((D+eta,(D+eta)*u),axis=-1)
    return cell_average(grid,f)


def linear_reference(grid,t,A,D,L,g):
    c=np.sqrt(g*D)
    return cell_average(grid,lambda x:A*np.sin(2*np.pi*(x-c*t)/L))


def peak_x(x,y,length):
    i=int(np.argmax(y));ym,yp=y[(i-1)%len(y)],y[(i+1)%len(y)]
    denom=ym-2*y[i]+yp
    shift=0 if abs(denom)<1e-30 else float(np.clip(.5*(ym-yp)/denom,-.5,.5))
    return float((x[i]+shift*(x[1]-x[0]))%length)


def verification(figures,results):
    rows=[]
    for model in ['linear','nonlinear']:
        errors=[]
        for n in [100,200,400,800]:
            grid=Grid(0,1000,n);D=10;A=.2;t=20;cfg=SolverConfig(boundary='periodic',cfl=.4)
            if model=='linear':
                initial=linear_initial(grid,A,D,1000,9.81)
                exact=linear_reference(grid,t,A,D,1000,9.81)
                run=simulate_linear(grid,D,initial[:,0],initial[:,1],[0,t],cfg)
            else:
                initial=cell_average(grid,lambda x:sine_simple_wave(x,0,D,A,1000))
                exact=cell_average(grid,lambda x:sine_simple_wave(x,t,D,A,1000))[:,0]-D
                run=simulate_nonlinear_flat(grid,D,initial[:,0],initial[:,1],[0,t],cfg)
            error=float(np.mean(np.abs(run.eta[-1]-exact)));errors.append(error)
            order=None if len(errors)==1 else float(np.log2(errors[-2]/error))
            rows.append({'model':model,'cells':n,'dx_m':grid.dx,'mean_L1_surface_error_m':error,
                'Linf_surface_error_m':float(np.max(np.abs(run.eta[-1]-exact))),
                'observed_L1_order':order,'mass_balance_residual':run.diagnostics['maximum_normalized_mass_balance_residual']})
    write_csv(results/'v02_convergence.csv',rows)
    fig,ax=plt.subplots(figsize=(8,4.5),constrained_layout=True)
    for model,color in [('linear','#0077b6'),('nonlinear','#d1495b')]:
        r=[row for row in rows if row['model']==model]
        dx=np.array([row['dx_m'] for row in r]);err=np.array([row['mean_L1_surface_error_m'] for row in r])
        ax.loglog(dx,err,'o-',color=color,label=f'{model.capitalize()} versus independent exact solution')
    ax.loglog(dx,err[-1]*(dx/dx[-1])**2,'k--',label='Second-order reference slope')
    ax.set(xlabel='Cell width (m)',ylabel='Mean absolute surface error (m)',title='MUSCL + SSP-RK2 | smooth periodic verification')
    ax.legend(fontsize=8);ax.grid(alpha=.2,which='both')
    fig.savefig(figures/'v02_convergence.png',dpi=180);plt.close(fig)
    return rows


def comparison(cfg,figures,results):
    D,L,g=cfg['depth_m'],cfg['periodic_length_m'],cfg['gravity_m_s2']
    c=np.sqrt(g*D);t=cfg['final_time_periods']*L/c
    solver_cfg=SolverConfig(g=g,cfl=cfg['cfl'],boundary='periodic')
    rows,detail,finest={},[],{}
    snapshots={}
    for eps in cfg['epsilon_values']:
        A=eps*D;results_by_mesh={}
        mesh_cells=list(cfg['grid_cells'])
        if eps==min(cfg['epsilon_values']):mesh_cells.append(cfg['extra_grid_for_smallest_epsilon'])
        for n in mesh_cells:
            grid=Grid(0,L,n)
            wi=linear_initial(grid,A,D,L,g);zi=matched_initial(grid,A,D,L,g)
            linear=simulate_linear(grid,D,wi[:,0],wi[:,1],[0,t],solver_cfg,order=cfg['order'])
            nonlinear=simulate_nonlinear_flat(grid,D,zi[:,0],zi[:,1],[0,t],solver_cfg,order=cfg['order'])
            le,ne=linear.eta[-1],nonlinear.eta[-1]
            reference=linear_reference(grid,t,A,D,L,g)
            results_by_mesh[n]=(linear,nonlinear)
            d=rms(ne-le)/A
            fine_sensitivity=None
            current_shift=(peak_x(grid.x,ne,L)-peak_x(grid.x,le,L)+L/2)%L-L/2
            crest_sensitivity=None
            if n!=mesh_cells[0]:
                prev=mesh_cells[mesh_cells.index(n)-1]
                ratio=n//prev
                if ratio*prev!=n:raise ValueError('Nested meshes required')
                old_l,old_n=results_by_mesh[prev]
                # Restrict fine cell averages conservatively to coarse grid.
                fine_sensitivity=(rms(le.reshape(prev,ratio).mean(axis=1)-old_l.eta[-1])+
                    rms(ne.reshape(prev,ratio).mean(axis=1)-old_n.eta[-1]))/A
                previous_shift=(peak_x(old_n.x,old_n.eta[-1],L)-peak_x(old_l.x,old_l.eta[-1],L)+L/2)%L-L/2
                crest_sensitivity=abs(current_shift-previous_shift)
            detail.append({'epsilon':eps,'cells':n,'model_difference_NRMSE':d,
                'linear_exact_NRMSE':rms(le-reference)/A,
                'combined_mesh_sensitivity':fine_sensitivity,
                'crest_shift_m':current_shift,'crest_mesh_sensitivity_m':crest_sensitivity,
                'linear_mass_residual':linear.diagnostics['maximum_normalized_mass_balance_residual'],
                'nonlinear_mass_residual':nonlinear.diagnostics['maximum_normalized_mass_balance_residual'],
                'linear_max_cfl':linear.diagnostics['max_actual_stage_cfl'],
                'nonlinear_max_cfl':nonlinear.diagnostics['max_actual_stage_cfl']})
        n=mesh_cells[-1];lr,nr=results_by_mesh[n];le,ne=lr.eta[-1],nr.eta[-1]
        x=lr.x
        linpeak,nonpeak=peak_x(x,le,L),peak_x(x,ne,L)
        shift=(nonpeak-linpeak+L/2)%L-L/2
        sensitivity=detail[-1]['combined_mesh_sensitivity']
        difference=detail[-1]['model_difference_NRMSE']
        dominance=difference/sensitivity if sensitivity>0 else None
        row={'epsilon':eps,'incident_amplitude_m':A,'finest_cells':n,
            'model_difference_NRMSE':difference,'linear_exact_NRMSE':detail[-1]['linear_exact_NRMSE'],
            'combined_mesh_sensitivity':sensitivity,'difference_to_sensitivity_ratio':dominance,
            'resolved_by_mesh_check':bool(dominance is not None and dominance>=cfg['resolved_sensitivity_ratio_threshold']),
            'crest_mesh_sensitivity_m':detail[-1]['crest_mesh_sensitivity_m'],
            'crest_shift_resolved_by_mesh_check':bool(abs(shift)>=cfg['resolved_sensitivity_ratio_threshold']*max(detail[-1]['crest_mesh_sensitivity_m'],1e-14)),
            'linear_peak_m':float(le.max()),'nonlinear_peak_m':float(ne.max()),
            'signed_peak_difference_over_A':float((ne.max()-le.max())/A),
            'nonlinear_peak_ahead_of_linear_m':float(shift),
            'mass_balance_residual_max':max(lr.diagnostics['maximum_normalized_mass_balance_residual'],nr.diagnostics['maximum_normalized_mass_balance_residual'])}
        rows[eps]=row;finest[eps]=(lr,nr)
        snapshots[f'epsilon_{eps:g}_linear_eta_m']=le
        snapshots[f'epsilon_{eps:g}_nonlinear_eta_m']=ne
        snapshots[f'epsilon_{eps:g}_x_m']=x
    summary=list(rows.values())
    write_csv(results/'v02_amplitude_comparison.csv',summary)
    write_csv(results/'v02_mesh_sensitivity.csv',detail)
    np.savez_compressed(results/'v02_finest_profiles.npz',time_s=np.array(t),**snapshots)
    fig,axes=plt.subplots(1,3,figsize=(13,4),constrained_layout=True)
    for ax,eps in zip(axes,[.001,.05,.2]):
        lr,nr=finest[eps];A=eps*D
        ax.plot(lr.x/L,lr.eta[-1]/A,color='#0077b6',label='Linear')
        ax.plot(nr.x/L,nr.eta[-1]/A,color='#d1495b',label='Nonlinear')
        ax.set(xlabel='x / L',ylabel='Surface elevation / A',title=f'A / D = {eps:g}')
        ax.grid(alpha=.2);ax.legend(fontsize=8)
    fig.suptitle(f'Matched physical initial data | t = {cfg["final_time_periods"]} L / c | N = {n}')
    fig.savefig(figures/'v02_model_comparison.png',dpi=180);plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(11,4.3),constrained_layout=True)
    eps=np.array([r['epsilon'] for r in summary]);diff=np.array([r['model_difference_NRMSE'] for r in summary]);sens=np.array([r['combined_mesh_sensitivity'] for r in summary])
    axes[0].loglog(eps,diff,'o-',color='#d1495b',label='Model waveform difference / A')
    axes[0].loglog(eps,sens,'s--',color='#555555',label='Combined last-two-grid sensitivity / A')
    axes[0].set(xlabel='Incident relative amplitude A / D',ylabel='Normalized RMS',title='Model difference versus mesh sensitivity')
    axes[0].legend(fontsize=8);axes[0].grid(alpha=.2,which='both')
    axes[1].errorbar(eps,[r['nonlinear_peak_ahead_of_linear_m'] for r in summary],
        yerr=[r['crest_mesh_sensitivity_m'] for r in summary],fmt='o-',color='#005f73',capsize=3,
        label='Bars: two-grid shift change (not confidence intervals)')
    axes[1].set_xscale('log');axes[1].legend(fontsize=7)
    axes[1].set(xlabel='Incident relative amplitude A / D',ylabel='Nonlinear crest ahead of linear crest (m)',title='Interpolated crest-position difference')
    axes[1].grid(alpha=.2,which='both')
    fig.savefig(figures/'v02_amplitude_effect.png',dpi=180);plt.close(fig)
    return {'configuration':cfg,'final_time_s':float(t),'summary':summary,'mesh_detail':detail,
        'interpretation':'Mesh sensitivity is empirical, not a rigorous numerical error bound. Both models omit dispersion; no experimental validation.'}


def animation(cfg,figures):
    D,L,g=cfg['depth_m'],cfg['periodic_length_m'],cfg['gravity_m_s2'];A=.2*D
    grid=Grid(0,L,800);c=np.sqrt(g*D);times=np.linspace(0,cfg['final_time_periods']*L/c,51)
    wi=linear_initial(grid,A,D,L,g);zi=matched_initial(grid,A,D,L,g)
    linear=simulate_linear(grid,D,wi[:,0],wi[:,1],times)
    nonlinear=simulate_nonlinear_flat(grid,D,zi[:,0],zi[:,1],times)
    fig,ax=plt.subplots(figsize=(9,4),constrained_layout=True)
    ll,=ax.plot(grid.x/L,linear.eta[0]/A,label='Linear',color='#0077b6',lw=2)
    nl,=ax.plot(grid.x/L,nonlinear.eta[0]/A,label='Nonlinear',color='#d1495b',lw=2)
    label=ax.text(.03,.9,'',transform=ax.transAxes)
    ax.set(xlim=(0,1),ylim=(-1.15,1.15),xlabel='x / L',ylabel='Surface elevation / A',title='Same initial surface and velocity | A / D = 0.2 | N = 800')
    ax.legend(loc='lower right');ax.grid(alpha=.2)
    def update(i):
        ll.set_ydata(linear.eta[i]/A);nl.set_ydata(nonlinear.eta[i]/A);label.set_text(f't = {times[i]:.2f} s');return ll,nl,label
    anim=FuncAnimation(fig,update,frames=len(times),blit=False)
    anim.save(figures/'v02_linear_nonlinear.gif',writer=PillowWriter(fps=12),dpi=100);plt.close(fig)


def main():
    start=time.perf_counter();results=ROOT/'results';figures=ROOT/'figures'
    results.mkdir(exist_ok=True);figures.mkdir(exist_ok=True)
    cfg=json.loads((ROOT/'configs/flat_comparison.json').read_text())
    plt.rcParams.update({'font.family':'DejaVu Sans','axes.spines.top':False,'axes.spines.right':False})
    convergence=verification(figures,results)
    comparisons=comparison(cfg,figures,results)
    animation(cfg,figures)
    write_json(results/'v02_metrics.json',{'convergence':convergence,'comparison':comparisons})
    source_hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()
        for pattern in ['src/**/*.py','scripts/*.py','tests/*.py','configs/*.json'] for p in ROOT.glob(pattern)}
    write_json(results/'v02_run_manifest.json',{'version':__version__,'created_utc':datetime.now(timezone.utc).isoformat(),
        'python':platform.python_version(),'numpy':np.__version__,'matplotlib':matplotlib.__version__,
        'source_sha256':source_hashes,'elapsed_seconds':time.perf_counter()-start,'command':'python scripts/run_v02.py',
        'parameter_combinations':len(cfg['epsilon_values'])*len(cfg['grid_cells'])+1,'model_runs_in_sweep':2*(len(cfg['epsilon_values'])*len(cfg['grid_cells'])+1),
        'scope':'flat-bed periodic numerical model comparison; no coastal shoaling or observations'})
    print(json.dumps({'convergence':convergence,'amplitude_summary':comparisons['summary']},indent=2))


if __name__=='__main__':main()
