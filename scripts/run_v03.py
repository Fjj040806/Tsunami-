"""Reproduce variable-bed verification tables and figures."""
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
import json,csv,hashlib,platform,time
from datetime import datetime,timezone
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from tsunami import Grid,SolverConfig,__version__
from tsunami.initial import cell_average
from tsunami.variable import simulate_linear_variable,simulate_nonlinear_variable,_reconstruct
from tsunami.manufactured import setup,exact_state,L,D0,B,A,G,C

def dump(path,data):path.write_text(json.dumps(data,indent=2,allow_nan=False),encoding='utf-8')
def table(path,rows):
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

def main():
    start=time.perf_counter();figures=ROOT/'figures';results=ROOT/'results'
    figures.mkdir(exist_ok=True);results.mkdir(exist_ok=True)
    cfg=json.loads((ROOT/'configs/variable_verification.json').read_text())
    # Fail rather than silently label a hard-coded analytic field with other parameters.
    assert cfg['gravity_m_s2']==G and cfg['manufactured_length_m']==L and cfg['manufactured_amplitude_m']==A
    solver=SolverConfig(g=cfg['gravity_m_s2'],cfl=cfg['cfl'],boundary='periodic')
    grid=Grid(0,cfg['lake_length_m'],cfg['lake_cells'])
    times=np.linspace(0,cfg['lake_final_time_s'],21);surface=cfg['lake_surface_m']
    lake=[];history={};beds={}
    for kind in ['smooth','step']:
        if kind=='smooth':d=cell_average(grid,lambda x:10+2*np.cos(2*np.pi*x/grid.right));mode='centered'
        else:d=np.where(grid.x<500,10.,4.);mode='constant'
        beds[kind]=-d
        for model in ['linear','nonlinear']:
            if model=='linear':r=simulate_linear_variable(grid,d,np.full(grid.cells,surface),np.zeros(grid.cells),times,solver,mode)
            else:r=simulate_nonlinear_variable(grid,-d,d+surface,np.zeros(grid.cells),times,solver,mode)
            drift=np.max(np.abs(r.eta-surface),axis=1);vel=np.max(np.abs(r.velocity),axis=1)
            history[f'{kind}_{model}_surface_drift_m']=drift
            history[f'{kind}_{model}_velocity_m_s']=vel
            lake.append({'bed':kind,'model':model,'cells':grid.cells,'duration_s':float(times[-1]),
                'maximum_surface_drift_m':float(drift.max()),'maximum_velocity_m_s':float(vel.max()),
                'mass_balance_residual':r.diagnostics['maximum_normalized_mass_balance_residual'],
                'maximum_actual_cfl':r.diagnostics['maximum_actual_cfl'],
                'minimum_depth_m':r.diagnostics['minimum_depth_m']})
    table(results/'v03_lake_at_rest.csv',lake)
    np.savez_compressed(results/'v03_lake_history.npz',time_s=times,x_m=grid.x,**history,**{f'{k}_bed_m':v for k,v in beds.items()})
    # Instantaneous momentum RHS if the centered source were omitted.
    w=np.column_stack((np.full(grid.cells,surface),np.zeros(grid.cells),beds['smooth']))
    lc,rc=_reconstruct(w,'periodic','centered')
    hl=lc[1:-1,0]-lc[1:-1,2];hr=rc[1:-1,0]-rc[1:-1,2]
    centered=.5*solver.g*(hl+hr)*(lc[1:-1,2]-rc[1:-1,2])/grid.dx
    omission=float(np.max(np.abs(centered)))
    rows=[];profiles={}
    for model in ['linear','nonlinear']:
        old=None
        for n in cfg['manufactured_cells']:
            grid=Grid(0,L,n);d,z,f=setup(grid,model)
            fn=simulate_linear_variable if model=='linear' else simulate_nonlinear_variable
            final=cfg['manufactured_final_time_s']
            r=fn(grid,d if model=='linear' else -d,z[:,0],z[:,1],[0,final],solver,forcing=f)
            exact=cell_average(grid,lambda x:exact_state(x,final,model))
            err=r.states[-1]-exact;e=float(np.mean(np.abs(err[:,0])))
            rows.append({'model':model,'cells':n,'dx_m':grid.dx,'surface_mean_L1_error_m':e,
                'surface_Linf_error_m':float(np.max(np.abs(err[:,0]))),
                'second_state_mean_L1_error':float(np.mean(np.abs(err[:,1]))),
                'observed_surface_L1_order':None if old is None else float(np.log2(old/e)),
                'mass_balance_residual':r.diagnostics['maximum_normalized_mass_balance_residual'],
                'maximum_actual_cfl':r.diagnostics['maximum_actual_cfl'],'minimum_depth_m':r.diagnostics['minimum_depth_m']})
            old=e
            if n==cfg['manufactured_cells'][-1]:
                profiles[model+'_eta_m']=r.eta[-1]
                profiles[model+'_exact_eta_m']=exact[:,0] if model=='linear' else exact[:,0]-d
                profiles[model+'_velocity_m_s']=r.velocity[-1]
    table(results/'v03_manufactured_convergence.csv',rows)
    np.savez_compressed(results/'v03_manufactured_profiles.npz',x_m=grid.x,**profiles)
    plt.rcParams.update({'axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(1,2,figsize=(11,4.3),constrained_layout=True)
    for kind,c in [('smooth','#0077b6'),('step','#d1495b')]:
        axes[0].plot(np.linspace(0,1000,cfg['lake_cells'],endpoint=False)+500/cfg['lake_cells'],beds[kind],color=c,label=kind.capitalize()+' bed')
        for model,style in [('linear','-'),('nonlinear','--')]:
            axes[1].plot(times,history[f'{kind}_{model}_surface_drift_m'],style,color=c,label=f'{kind}, {model}')
    axes[0].axhline(surface,color='k',label='Still water surface')
    axes[0].set(xlabel='x (m)',ylabel='Elevation (m)',title='Fully wet static beds')
    axes[1].set(xlabel='Time (s)',ylabel='Maximum surface drift (m)',title='Static water: residual at floating-point scale')
    for ax in axes:ax.grid(alpha=.2);ax.legend(fontsize=8)
    fig.savefig(figures/'v03_lake_at_rest.png',dpi=180);plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(11,4.3),constrained_layout=True)
    for model,c in [('linear','#0077b6'),('nonlinear','#d1495b')]:
        rr=[r for r in rows if r['model']==model];dx=np.array([r['dx_m'] for r in rr]);err=np.array([r['surface_mean_L1_error_m'] for r in rr])
        axes[0].loglog(dx,err,'o-',color=c,label=model.capitalize())
        axes[1].plot(grid.x,profiles[model+'_eta_m'],color=c,label=model.capitalize()+' numerical')
    axes[0].loglog(dx,err[-1]*(dx/dx[-1])**2,'k--',label='Second-order reference')
    axes[1].plot(grid.x,profiles['linear_exact_eta_m'],'k--',label='Manufactured exact surface')
    axes[0].set(xlabel='Cell width (m)',ylabel='Mean L1 surface error (m)',title='Smooth variable-bed manufactured solution')
    axes[1].set(xlabel='x (m)',ylabel='Surface elevation (m)',title=f'Forced verification only | N={grid.cells}')
    for ax in axes:ax.grid(alpha=.2);ax.legend(fontsize=8)
    fig.savefig(figures/'v03_manufactured_convergence.png',dpi=180);plt.close(fig)
    metrics={'configuration':cfg,'lake_at_rest':lake,'manufactured_convergence':rows,
        'omitting_centered_source_instantaneous_max_q_rhs_m2_s2':omission,
        'forcing_parameters':{'D0_m':D0,'B_m':B,'A_m':A,'L_m':L,'c_m_s':float(C)},
        'scope':'Numerical verification only. Manufactured forcing is artificial. No observational validation, shoaling experiment or wet/dry capability.'}
    dump(results/'v03_metrics.json',metrics)
    dump(results/'v03_run_manifest.json',{'version':__version__,'created_utc':datetime.now(timezone.utc).isoformat(),
        'python':platform.python_version(),'numpy':np.__version__,'matplotlib':matplotlib.__version__,
        'command':'python scripts/run_v03.py','elapsed_seconds':time.perf_counter()-start,
        'source_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()
          for pattern in ['src/**/*.py','scripts/*.py','tests/*.py','configs/*.json'] for p in ROOT.glob(pattern)}})
    print(json.dumps(metrics,indent=2))
if __name__=='__main__':main()
