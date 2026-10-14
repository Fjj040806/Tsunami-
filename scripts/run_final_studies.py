"""Run reproducible step, shelf and conditional laboratory comparisons."""
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
import argparse,csv,json,hashlib,zipfile,io,time
from datetime import datetime,timezone
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from tsunami import Grid,SolverConfig,__version__,IncidentBoundary
from tsunami.variable import simulate_linear_variable,simulate_nonlinear_variable
from tsunami.initial import cell_average
from tsunami.studies import gaussian,shelf_depth,step_reference,composite_depth,G

def dump(p,v):p.write_text(json.dumps(v,indent=2,allow_nan=False),encoding='utf-8')
def table(p,rows):
    with p.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def norm(x):return float(np.linalg.norm(x))
def solve(grid,d,amp,center,sigma,t,model,mode='centered',depth_function=None,step_coupling='hydrostatic'):
    df=depth_function if depth_function is not None else lambda x:np.interp(x,grid.x,d)
    # The incoming pulse is wholly on the deep plateau to numerical precision.
    deep=float(d[0]);c=np.sqrt(G*deep)
    def init(x):
        e=gaussian(x,amp,center,sigma);u=c/deep*e
        if model=='linear':return np.stack((e,u),axis=-1)
        h=df(x)+e;return np.stack((h,h*u),axis=-1)
    z=cell_average(grid,init)
    fn=simulate_linear_variable if model=='linear' else simulate_nonlinear_variable
    return fn(grid,d if model=='linear' else -d,z[:,0],z[:,1],[0,t],SolverConfig(cfl=.4,boundary='outflow'),mode,step_coupling=step_coupling)

def step(cfg,results,figures):
    rows=[];last=None
    for n in cfg['cells']:
        grid=Grid(cfg['left_m'],cfg['right_m'],n)
        df=lambda x:np.where(x<0,cfg['deep_m'],cfg['shallow_m'])
        d=df(grid.x)
        ex=cell_average(grid,lambda x:step_reference(x,cfg['final_time_s'],cfg['deep_m'],cfg['shallow_m'],cfg['amplitude_m'],cfg['center_m'],cfg['sigma_m']))[:,0]
        for model in ['linear','nonlinear','nonlinear_hydrostatic']:
            r=solve(grid,d,cfg['amplitude_m'],cfg['center_m'],cfg['sigma_m'],cfg['final_time_s'],model,'constant',df,'bernoulli' if model=='nonlinear' else 'hydrostatic')
            eta=r.eta[-1];R=float(eta[grid.x<0].max()/cfg['amplitude_m']);T=float(eta[grid.x>0].max()/cfg['amplitude_m'])
            cl,cr=np.sqrt(G*cfg['deep_m']),np.sqrt(G*cfg['shallow_m'])
            rows.append({'model':model,'cells':n,'R_measured':R,'R_linear_reference':float((cl-cr)/(cl+cr)),
                'T_measured':T,'T_linear_reference':float(2*cl/(cl+cr)),
                'relative_L2_to_linear_reference':norm(eta-ex)/norm(ex),
                'normalized_mass_residual':r.diagnostics['maximum_normalized_mass_balance_residual'],
                'minimum_depth_m':r.diagnostics['minimum_depth_m']})
            if n==cfg['cells'][-1]:
                if last is None:last={'x_m':grid.x,'exact_eta_m':ex,'bed_m':-d}
                last[model+'_eta_m']=eta
        print('step',n,'complete',flush=True)
    table(results/'step_scattering.csv',rows);np.savez_compressed(results/'step_profiles.npz',**last)
    fig,axes=plt.subplots(1,2,figsize=(11,4),constrained_layout=True)
    A=cfg['amplitude_m']
    for name,c in [('linear','#0077b6'),('nonlinear','#d1495b'),('nonlinear_hydrostatic','#888888')]:axes[0].plot(last['x_m'],last[name+'_eta_m']/A,label=name,color=c)
    axes[0].plot(last['x_m'],last['exact_eta_m']/A,'k--',label='Linear exact scattering')
    axes[0].set(xlabel='x (m)',ylabel='Surface / incident amplitude',title='Low-amplitude step scattering');axes[0].legend(fontsize=8)
    for model,c in [('linear','#0077b6'),('nonlinear','#d1495b'),('nonlinear_hydrostatic','#888888')]:
        rr=[r for r in rows if r['model']==model]
        axes[1].loglog([r['cells'] for r in rr],[r['relative_L2_to_linear_reference'] for r in rr],'o-',color=c,label=model)
    axes[1].set(xlabel='Cells',ylabel='Relative L2 difference from linear exact solution',title='Grid refinement; nonlinear retains model difference')
    axes[1].legend(fontsize=8)
    for ax in axes:ax.grid(alpha=.2)
    fig.savefig(figures/'step_scattering.png',dpi=180);plt.close(fig)
    dump(results/'step_metrics.json',{'configuration':cfg,'rows':rows,'boundary_note':'Pulse remains far from domain ends; Gaussian tails are not compactly supported.'})

def shelf(cfg,results,figures,threshold):
    rows=[];detail=[];profiles={};case=0
    for ratio in cfg['depth_ratio']:
      for length in cfg['transition_length_m']:
       for eps in cfg['epsilon']:
        case+=1;previous=None;mesh=[]
        for n in cfg['cells']:
            grid=Grid(cfg['left_m'],cfg['right_m'],n);deep=cfg['deep_m'];shallow=deep*ratio;A=eps*deep
            df=lambda x:shelf_depth(x,deep,shallow,length)
            d=cell_average(grid,df)
            lr=solve(grid,d,A,cfg['center_m'],cfg['sigma_m'],cfg['final_time_s'],'linear','centered',df)
            nr=solve(grid,d,A,cfg['center_m'],cfg['sigma_m'],cfg['final_time_s'],'nonlinear','centered',df)
            le,ne=lr.eta[-1],nr.eta[-1];scale=norm(le)
            E=norm(ne-le)/scale;S=None
            if previous is not None:
                old_l,old_n=previous
                factor=n//len(old_l);restricted_l=le.reshape(len(old_l),factor).mean(axis=1)
                restricted_n=ne.reshape(len(old_n),factor).mean(axis=1)
                S=(norm(restricted_l-old_l)+norm(restricted_n-old_n))/norm(restricted_l)
            mask=grid.x>length/2+3*cfg['sigma_m']
            maskleft=grid.x<-length/2-3*cfg['sigma_m']
            record={'depth_ratio':ratio,'transition_length_m':length,'epsilon':eps,'cells':n,
                'relative_L2_model_difference':E,'combined_mesh_sensitivity':S,
                'linear_transmitted_peak_over_A':float(le[mask].max()/A),
                'nonlinear_transmitted_peak_over_A':float(ne[mask].max()/A),
                'linear_reflected_peak_over_A':float(le[maskleft].max()/A),
                'nonlinear_reflected_peak_over_A':float(ne[maskleft].max()/A),
                'linear_mass_residual':lr.diagnostics['maximum_normalized_mass_balance_residual'],
                'nonlinear_mass_residual':nr.diagnostics['maximum_normalized_mass_balance_residual'],
                'minimum_nonlinear_depth_m':nr.diagnostics['minimum_depth_m'],
                'maximum_actual_cfl':max(lr.diagnostics['maximum_actual_cfl'],nr.diagnostics['maximum_actual_cfl'])}
            detail.append(record);mesh.append(record);previous=(le,ne)
        row=mesh[-1].copy();row['difference_to_mesh_ratio']=row['relative_L2_model_difference']/row['combined_mesh_sensitivity']
        row['passes_empirical_mesh_screen']=row['difference_to_mesh_ratio']>=threshold
        row['green_slow_variation_reference']=float(ratio**(-.25))
        row['linear_peak_last_two_mesh_change']=abs(mesh[-1]['linear_transmitted_peak_over_A']-mesh[-2]['linear_transmitted_peak_over_A'])
        row['nonlinear_peak_last_two_mesh_change']=abs(mesh[-1]['nonlinear_transmitted_peak_over_A']-mesh[-2]['nonlinear_transmitted_peak_over_A'])
        rows.append(row)
        key=f'r{ratio:g}_l{length:g}_e{eps:g}'
        if eps==max(cfg['epsilon']):profiles[key+'_linear']=le;profiles[key+'_nonlinear']=ne
        table(results/'shelf_summary.csv',rows);table(results/'shelf_mesh_detail.csv',detail)
        print('shelf case',case,'ratio',ratio,'length',length,'epsilon',eps,'E',E,'screen',row['passes_empirical_mesh_screen'],flush=True)
    np.savez_compressed(results/'shelf_profiles.npz',x_m=grid.x,**profiles)
    fig,axes=plt.subplots(1,2,figsize=(11,4.3),constrained_layout=True)
    for ratio,c in [(.5,'#0077b6'),(.25,'#d1495b')]:
        for length,marker in [(200,'o'),(800,'s'),(2400,'^')]:
            rr=[r for r in rows if r['depth_ratio']==ratio and r['transition_length_m']==length]
            axes[0].loglog([r['epsilon'] for r in rr],[r['relative_L2_model_difference'] for r in rr],marker+'-',color=c,label=f'Ds/Dd={ratio:g}, length={length} m')
        rr=[r for r in rows if r['depth_ratio']==ratio and r['epsilon']==min(cfg['epsilon'])]
        axes[1].errorbar([r['transition_length_m']/cfg['sigma_m'] for r in rr],[r['linear_transmitted_peak_over_A'] for r in rr],yerr=[r['linear_peak_last_two_mesh_change'] for r in rr],fmt='o-',color=c,label=f'Linear, Ds/Dd={ratio:g}')
        axes[1].axhline(ratio**(-.25),color=c,ls='--',label=f'Green limit, Ds/Dd={ratio:g}')
    axes[0].set(xlabel='Incident amplitude / deep depth',ylabel='Relative L2 model difference',title='Shelf: nonlinearity and geometry')
    axes[1].set(xlabel='Transition length / initial Gaussian sigma',ylabel='Transmitted peak / incident amplitude',title='Low-amplitude linear shoaling (bars: mesh changes)')
    for ax in axes:ax.legend(fontsize=7);ax.grid(alpha=.2,which='both')
    fig.savefig(figures/'shelf_parameter_study.png',dpi=180);plt.close(fig)
    fig,axes=plt.subplots(1,3,figsize=(13,4),constrained_layout=True)
    for ax,length in zip(axes,cfg['transition_length_m']):
        key=f'r0.25_l{length:g}_e0.05';A=.05*cfg['deep_m']
        ax.plot(grid.x,profiles[key+'_linear']/A,label='Linear',color='#0077b6');ax.plot(grid.x,profiles[key+'_nonlinear']/A,label='Nonlinear',color='#d1495b')
        ax.set(xlim=(-5000,5000),xlabel='x (m)',ylabel='Surface / A',title=f'Transition {length} m');ax.legend(fontsize=8);ax.grid(alpha=.2)
    fig.suptitle('Depth ratio 0.25; incident amplitude/depth 0.05; t=450 s')
    fig.savefig(figures/'shelf_waveforms.png',dpi=180);plt.close(fig)
    dump(results/'shelf_metrics.json',{'configuration':cfg,'rows':rows,'mesh_screen_threshold':threshold,
        'interpretation':'Model differences are not observational errors. Two-grid sensitivity is empirical, not an upper bound. Green reference assumes slow variation and weak reflection.'})

def load_archive(path):
    with zipfile.ZipFile(path) as z:
        lines=z.read('ts3a.txt').decode('latin1').splitlines();rows=[]
        for line in lines:
            try:r=[float(x) for x in line.split()]
            except ValueError:continue
            if len(r)==8:rows.append(r)
    a=np.array(rows)
    if a.shape!=(600,8) or not np.isfinite(a).all() or not np.allclose(np.diff(a[:,0]),.05):raise ValueError('Unexpected NOAA data layout')
    return a

def lab(cfg,results,figures):
    archive=ROOT/'data/raw/noaa_composite.zip'
    if not archive.exists():
        import subprocess
        subprocess.run([sys.executable,str(ROOT/'scripts/fetch_noaa_data.py')],check=True)
    a=load_archive(archive);baseline=a[:cfg['baseline_samples'],1:].mean(axis=0);eta=a[:,1:]-baseline
    absolute=a[:,0];times=absolute-absolute[0];keep=absolute<=cfg['final_absolute_time_s'];times=times[keep]
    coords=[12.64,15.04,17.22,19.40,20.86,22.33,22.80]
    incident=lambda t:float(np.interp(t,absolute-absolute[0],eta[:,0],left=0,right=0))
    rows=[];history={};previous={}
    for n in cfg['cells']:
        grid=Grid(cfg['upstream_x_m'],cfg['wall_x_m'],n);d=cell_average(grid,composite_depth)
        for model in ['linear','nonlinear']:
            fn=simulate_linear_variable if model=='linear' else simulate_nonlinear_variable
            r=fn(grid,d if model=='linear' else -d,np.zeros(n) if model=='linear' else d,np.zeros(n),times,
                SolverConfig(cfl=.4,boundary='outflow'),boundary=IncidentBoundary(incident))
            pred=np.column_stack([np.array([np.interp(x,grid.x,e) for e in r.eta]) for x in coords[1:]])
            old=previous.get(model)
            for j,x in enumerate(coords[1:]):
                obs=eta[keep,j+1];diff=pred[:,j]-obs
                scale=float(np.max(np.abs(obs)))
                rows.append({'model':model,'cells':n,'gauge':j+5,'x_m':x,
                    'RMSE_m':float(np.sqrt(np.mean(diff**2))),'NRMSE_over_observed_peak':float(np.sqrt(np.mean(diff**2))/scale),
                    'predicted_peak_m':float(pred[:,j].max()),'observed_peak_m':float(obs.max()),
                    'peak_time_difference_s':float(times[np.argmax(pred[:,j])]-times[np.argmax(obs)]),
                    'normalized_two_grid_RMS_change':None if old is None else float(np.sqrt(np.mean((pred[:,j]-old[:,j])**2))/scale),
                    'mass_residual':r.diagnostics['maximum_normalized_mass_balance_residual'],
                    'minimum_depth_m':r.diagnostics['minimum_depth_m']})
            previous[model]=pred
            if n==cfg['cells'][-1]:history[model]=pred
        print('laboratory',n,'complete',flush=True)
    table(results/'lab_comparison.csv',rows)
    # Bundle predictions and timestamps, not third-party measured traces.
    np.savez_compressed(results/'lab_predictions.npz',absolute_time_s=absolute[keep],x_m=np.array(coords[1:]),**history)
    fig,axes=plt.subplots(2,3,figsize=(12,6.5),constrained_layout=True)
    for j,ax in enumerate(axes.flat):
        ax.plot(absolute[keep],eta[keep,j+1]*1000,'k.',ms=2,label='USACE observed (baseline removed)')
        for model,c in [('linear','#0077b6'),('nonlinear','#d1495b')]:ax.plot(absolute[keep],history[model][:,j]*1000,color=c,label=model)
        ax.set(xlabel='Original experiment time (s)',ylabel='Surface (mm)',title=f'G{j+5}, x={coords[j+1]:.2f} m');ax.grid(alpha=.2)
        if j==0:ax.legend(fontsize=6.5)
    fig.suptitle('Conditional Case A comparison: measured G4 input; no fitted time shift')
    fig.savefig(figures/'lab_case_a_comparison.png',dpi=180);plt.close(fig)
    dump(results/'lab_metrics.json',{'configuration':cfg,'rows':rows,'baseline_m':baseline.tolist(),
        'archive_sha256':hashlib.sha256(archive.read_bytes()).hexdigest(),
        'raw_rows':len(a),'comparison_samples':len(times),'gauge_coordinates_m':coords,
        'archive_url':'https://nctr.pmel.noaa.gov/benchmark/Solitary_wave/tsunami3_runup.zip',
        'geometry_junctions_m':[15.039,19.399,22.329],
        'wall_still_water_depth_m':float(composite_depth(23.22)),
        'metadata_conflicts':'Web G4 Case A 12.22 m versus archive 12.64 m; webpage amplitudes differ from archive actual. Archive README chosen, figure geometry converted to wavemaker origin; conflict remains disclosed.',
        'conditional_scope':'G4 provides incoming characteristic only; downstream gauge predictions are compared without parameter/time fitting. Not independent source reconstruction, complete official benchmark certification, or real-event validation.',
        'window_note':'Ends at 283 s; estimated peak reflection returns to G4 around 289 s; earlier tails and decomposition uncertainty are limitations.'})

def main():
    p=argparse.ArgumentParser();p.add_argument('--section',choices=['all','step','shelf','lab'],default='all');args=p.parse_args()
    cfg=json.loads((ROOT/'configs/final_studies.json').read_text());results=ROOT/'results';figures=ROOT/'figures'
    results.mkdir(exist_ok=True);figures.mkdir(exist_ok=True);start=time.perf_counter()
    for name in ['step','shelf','lab']:
        if args.section in ('all',name):
            if name=='shelf':shelf(cfg[name],results,figures,cfg['mesh_screen_ratio'])
            else:globals()[name](cfg[name],results,figures)
    dump(results/f'final_{args.section}_run_manifest.json',{'version':__version__,'created_utc':datetime.now(timezone.utc).isoformat(),
        'elapsed_seconds':time.perf_counter()-start,'command':f'python scripts/run_final_studies.py --section {args.section}',
        'source_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for pattern in ['src/**/*.py','scripts/*.py','configs/*.json'] for p in ROOT.glob(pattern)}})
if __name__=='__main__':main()
