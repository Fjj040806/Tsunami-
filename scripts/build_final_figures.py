"""Build publication-sized figures and LaTeX tables from measured result files."""
from pathlib import Path
import json,csv,sys
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from run_final_studies import load_archive

BLUE='#0077b6';PINK='#c83f5c';GRAY='#777777'
def read(name):return json.loads((ROOT/'results'/name).read_text())
def save(fig,name):fig.savefig(ROOT/'figures'/name,dpi=220,bbox_inches='tight');plt.close(fig)
def label(ax):ax.grid(alpha=.2,which='both')
def tab(path,header,rows,caption):
    columns='l'+'r'*(len(header)-1)
    text='\\begin{table}[htbp]\\centering\\small\n\\begin{tabular}{'+columns+'}\\toprule\n'
    text+=' & '.join(header)+' \\\\\n\\midrule\n'
    text+=''.join(' & '.join(row)+' \\\\\n' for row in rows)
    text+='\\bottomrule\\end{tabular}\n\\caption{'+caption+'}\\end{table}\n'
    path.write_text(text)
def main():
    plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False,'legend.fontsize':9})
    f=ROOT/'figures';p=ROOT/'paper';f.mkdir(exist_ok=True)
    verpath=ROOT/'results/v03_metrics.json'
    ver=json.loads((verpath if verpath.exists() else ROOT/'v0.3_results/v03_metrics.json').read_text())
    fig,ax=plt.subplots(figsize=(7,4.3),constrained_layout=True)
    for model,c,m in [('linear',BLUE,'o'),('nonlinear',PINK,'s')]:
        rr=[r for r in ver['manufactured_convergence'] if r['model']==model]
        dx=np.array([r['dx_m'] for r in rr]);err=np.array([r['surface_mean_L1_error_m'] for r in rr])
        ax.loglog(dx,err,m+'-',color=c,label=model.capitalize())
    ax.loglog(dx,err[-1]*(dx/dx[-1])**2,'k--',label='Second-order reference')
    ax.set(xlabel='Cell width (m)',ylabel='Mean L1 surface error (m)',title='Smooth variable-bed manufactured verification')
    ax.legend();label(ax);save(fig,'final_verification.png')
    st=read('step_metrics.json');z=np.load(ROOT/'results/step_profiles.npz')
    fig,ax=plt.subplots(figsize=(7,4.3),constrained_layout=True)
    for name,c,ls in [('linear',BLUE,'-'),('nonlinear',PINK,'-'),('nonlinear_hydrostatic',GRAY,':')]:
        ax.plot(z['x_m'],z[name+'_eta_m']/st['configuration']['amplitude_m'],color=c,ls=ls,label={'nonlinear':'Nonlinear, stationary coupling','nonlinear_hydrostatic':'Nonlinear, hydrostatic control'}.get(name,'Linear'))
    ax.plot(z['x_m'],z['exact_eta_m']/st['configuration']['amplitude_m'],'k--',label='Linear exact scattering')
    ax.set(xlim=(-3000,2500),xlabel='x (m)',ylabel='Surface / incident amplitude',title='20-to-5 m step; incident A/D = 0.0001')
    ax.legend(fontsize=8.5);label(ax);save(fig,'final_step_waveform.png')
    finest=max(r['cells'] for r in st['rows']);rr=[r for r in st['rows'] if r['cells']==finest]
    names={'linear':'Linear','nonlinear':'Nonlinear, junction','nonlinear_hydrostatic':'Nonlinear, HR control'}
    tab(p/'generated_step_table.tex',['Model','$R$','$T$','Relative $L^2$ (\\%)'],[[names[r['model']],f"{r['R_measured']:.5f}",f"{r['T_measured']:.5f}",f"{100*r['relative_L2_to_linear_reference']:.3f}"] for r in rr],f'Step results on {finest} cells. Exact linear coefficients: $R=1/3$, $T=4/3$.')
    sh=read('shelf_metrics.json');rows=sh['rows']
    fig,ax=plt.subplots(figsize=(7,4.6),constrained_layout=True)
    for ratio,c in [(.5,BLUE),(.25,PINK)]:
        for length,m in [(200,'o'),(800,'s'),(2400,'^')]:
            rr=[r for r in rows if r['depth_ratio']==ratio and r['transition_length_m']==length]
            ax.loglog([r['epsilon'] for r in rr],[100*r['relative_L2_model_difference'] for r in rr],color=c,label=f'Ds/Dd={ratio}, length={length} m')
            for r in rr:ax.plot(r['epsilon'],100*r['relative_L2_model_difference'],marker=m,color=c,mfc=c if r['passes_empirical_mesh_screen'] else 'white',ms=7)
    ax.set(xlabel='Incident amplitude / deep depth',ylabel='Relative L2 model difference (%)',title='Shelf comparison at 450 s; N=4000')
    ax.legend(fontsize=8.2,ncol=2,loc='lower right');label(ax)
    ax.text(.02,.98,'Open markers: empirical mesh screen fails',va='top',transform=ax.transAxes,fontsize=9)
    save(fig,'final_shelf_difference.png')
    tab(p/'generated_shelf_table.tex',['$D_s/D_d$','$\\ell$ (m)','$A/D_d$','$100E$','$100S$','Screen'],[[f"{r['depth_ratio']:.2f}",str(r['transition_length_m']),f"{r['epsilon']:.3f}",f"{100*r['relative_L2_model_difference']:.2f}",f"{100*r['combined_mesh_sensitivity']:.2f}",'Pass' if r['passes_empirical_mesh_screen'] else 'Unresolved'] for r in rows],'Standard-grid shelf study at 4000 cells. $E$ and $S$ have the denominators in equations (\\ref{eq:E}) and (\\ref{eq:S}); their percentages are not percentages of incident amplitude.')
    weak=[]
    for ratio in [.5,.25]:weak+=read(f'weak_metrics_r{ratio:g}.json')['rows']
    tab(p/'generated_weak_table.tex',['$D_s/D_d$','Cells','$100E$','$100S$','$E/S$'],[[f"{r['depth_ratio']:.2f}",str(r['cells']),f"{100*r['relative_L2_model_difference']:.3f}",'--' if r['combined_mesh_sensitivity'] is None else f"{100*r['combined_mesh_sensitivity']:.3f}",'--' if r['difference_to_mesh_ratio'] is None else f"{r['difference_to_mesh_ratio']:.2f}"] for r in weak],'Exploratory focused refinement at $\\ell=2400$ m and $A/D_d=0.001$. The two finest rows pass the fivefold screen; the other four weak physical cases have not received this follow-up.')
    fig,ax=plt.subplots(figsize=(7,4.4),constrained_layout=True)
    for ratio,c in [(.5,BLUE),(.25,PINK)]:
        rr=[r for r in rows if r['depth_ratio']==ratio and r['epsilon']==.001]
        ax.errorbar([r['transition_length_m']/150 for r in rr],[r['linear_transmitted_peak_over_A'] for r in rr],yerr=[r['linear_peak_last_two_mesh_change'] for r in rr],fmt='o-',color=c,capsize=3,label=f'N=4000; Ds/Dd={ratio}')
        ax.axhline(ratio**(-.25),ls='--',color=c,label=f'Green reference; Ds/Dd={ratio}')
        fine=[r for r in weak if r['depth_ratio']==ratio][-1];prev=[r for r in weak if r['depth_ratio']==ratio][-2]
        ax.errorbar([16],[fine['linear_transmitted_peak_over_A']],yerr=[abs(fine['linear_transmitted_peak_over_A']-prev['linear_transmitted_peak_over_A'])],fmt='D',color=c,mfc='white',capsize=4,label=f'N=16000 follow-up; ratio={ratio}')
    ax.set(xlabel='Transition length / incident Gaussian sigma',ylabel='Transmitted peak / incident amplitude',title='Weak-wave shoaling; bars are mesh changes')
    ax.legend(fontsize=8.5,ncol=2);label(ax);save(fig,'final_shoaling.png')
    lab=read('lab_metrics.json');pred=np.load(ROOT/'results/lab_predictions.npz')
    a=load_archive(ROOT/'data/raw/noaa_composite.zip');eta=a[:,1:]-np.array(lab['baseline_m']);keep=a[:,0]<=lab['configuration']['final_absolute_time_s']
    for indices,name in [([0,1,2],'final_lab_upstream.png'),([3,4,5],'final_lab_downstream.png')]:
        fig,axes=plt.subplots(3,1,figsize=(7,8.4),constrained_layout=True,sharex=True,sharey=True)
        for j,ax in zip(indices,axes):
            ax.plot(a[keep,0],eta[keep,j+1]*1000,'k.',ms=2.5,label='USACE observation, baseline removed')
            ax.plot(pred['absolute_time_s'],pred['linear'][:,j]*1000,color=BLUE,label='Linear')
            ax.plot(pred['absolute_time_s'],pred['nonlinear'][:,j]*1000,color=PINK,ls='--',label='Nonlinear')
            ax.set(ylabel='Surface (mm)',ylim=(-5,20),title=f'G{j+5}; wavemaker-origin x={pred["x_m"][j]:.2f} m');label(ax)
        axes[0].legend(fontsize=8);axes[-1].set(xlabel='Original experiment time (s)',xlim=(a[keep,0][0],a[keep,0][-1]))
        axes[-1].set_xticks([265.05,270,275,280,283]);save(fig,name)
    last=max(r['cells'] for r in lab['rows']);labrows=[r for r in lab['rows'] if r['cells']==last]
    paired=[]
    for gauge in range(5,11):
        l=next(r for r in labrows if r['gauge']==gauge and r['model']=='linear');n=next(r for r in labrows if r['gauge']==gauge and r['model']=='nonlinear')
        paired.append([f'G{gauge}',f"{l['observed_peak_m']*1000:.2f}",f"{l['NRMSE_over_observed_peak']*100:.2f}",f"{n['NRMSE_over_observed_peak']*100:.2f}",f"{l['normalized_two_grid_RMS_change']*100:.2f}",f"{n['normalized_two_grid_RMS_change']*100:.2f}"])
    tab(p/'generated_lab_table.tex',['Gauge','Obs. peak (mm)','L RMSE (\\%)','NL RMSE (\\%)','L mesh (\\%)','NL mesh (\\%)'],paired,'Conditional Case A comparison at 1600 cells. RMSE and last-two-grid RMS change are normalized by each observed absolute peak in the fixed window. These are not confidence intervals.')
    def rng(model):
        a=[100*r['NRMSE_over_observed_peak'] for r in labrows if r['model']==model]
        return f'{min(a):.2f}--{max(a):.2f}\\%'
    tests=read('test_summary.json')
    (p/'generated_macros.tex').write_text(f"\\newcommand{{\\TestCount}}{{{tests['tests_run']}}}\n\\newcommand{{\\ShelfPass}}{{{sum(r['passes_empirical_mesh_screen'] for r in rows)}}}\n\\newcommand{{\\LabSamples}}{{{lab['comparison_samples']}}}\n\\newcommand{{\\LabLinearRange}}{{{rng('linear')}}}\n\\newcommand{{\\LabNonlinearRange}}{{{rng('nonlinear')}}}\n")
    summary={'tests':tests['tests_run'],'step_finest':rr if False else [r for r in st['rows'] if r['cells']==finest],
        'shelf_standard_mesh_passes':sum(r['passes_empirical_mesh_screen'] for r in rows),'shelf_cases':len(rows),
        'weak_followup_finest':[r for r in weak if r['cells']==16000],
        'lab_finest':labrows,'laboratory_scope':lab['conditional_scope']}
    (ROOT/'results/final_summary.json').write_text(json.dumps(summary,indent=2))
    print('Final figures, LaTeX tables and summary generated from results.')
if __name__=='__main__':main()
