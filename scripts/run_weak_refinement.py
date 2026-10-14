"""Focused fine-mesh checks of two weak, slowly varying shelf cases."""
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'));sys.path.insert(0,str(ROOT/'src'))
import argparse,json,numpy as np
from run_final_studies import solve,table,dump,norm
from tsunami import Grid
from tsunami.initial import cell_average
from tsunami.studies import shelf_depth

def main():
    p=argparse.ArgumentParser();p.add_argument('--ratio',type=float,choices=[.5,.25],required=True);args=p.parse_args()
    cfg=json.loads((ROOT/'configs/final_studies.json').read_text())['shelf'];rows=[];previous=None
    for n in [4000,8000,16000]:
        grid=Grid(cfg['left_m'],cfg['right_m'],n);deep=cfg['deep_m'];shallow=deep*args.ratio;A=.001*deep;length=2400
        df=lambda x:shelf_depth(x,deep,shallow,length);d=cell_average(grid,df)
        lr=solve(grid,d,A,cfg['center_m'],cfg['sigma_m'],cfg['final_time_s'],'linear','centered',df)
        nr=solve(grid,d,A,cfg['center_m'],cfg['sigma_m'],cfg['final_time_s'],'nonlinear','centered',df)
        le,ne=lr.eta[-1],nr.eta[-1];E=norm(ne-le)/norm(le);S=None
        if previous is not None:
            ol,on=previous;factor=n//len(ol);rl=le.reshape(len(ol),factor).mean(axis=1);rn=ne.reshape(len(on),factor).mean(axis=1)
            S=(norm(rl-ol)+norm(rn-on))/norm(rl)
        mask=grid.x>length/2+3*cfg['sigma_m']
        rows.append({'depth_ratio':args.ratio,'transition_length_m':length,'epsilon':.001,'cells':n,
            'relative_L2_model_difference':E,'combined_mesh_sensitivity':S,
            'difference_to_mesh_ratio':None if S is None else E/S,
            'passes_empirical_mesh_screen':False if S is None else E/S>=5,
            'linear_transmitted_peak_over_A':float(le[mask].max()/A),'nonlinear_transmitted_peak_over_A':float(ne[mask].max()/A),
            'linear_mass_residual':lr.diagnostics['maximum_normalized_mass_balance_residual'],
            'nonlinear_mass_residual':nr.diagnostics['maximum_normalized_mass_balance_residual']})
        previous=(le,ne)
        table(ROOT/f'results/weak_refinement_r{args.ratio:g}.csv',rows)
        print(args.ratio,n,'E',E,'S',S,flush=True)
    np.savez_compressed(ROOT/f'results/weak_profiles_r{args.ratio:g}.npz',x_m=grid.x,linear_eta_m=le,nonlinear_eta_m=ne)
    dump(ROOT/f'results/weak_metrics_r{args.ratio:g}.json',{'rows':rows,'selection':'Two slowly varying cases chosen for focused follow-up after coarse low-amplitude mesh checks failed; exploratory, not a pre-registered selection.'})
if __name__=='__main__':main()
