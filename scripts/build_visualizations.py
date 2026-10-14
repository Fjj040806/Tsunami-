"""Build offline wave explorer and an illustrative GIF from modest-grid runs."""
from pathlib import Path
import sys,json,io
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image
from tsunami import Grid,SolverConfig,simulate_linear_variable,simulate_nonlinear_variable
from tsunami.initial import cell_average
from tsunami.studies import gaussian,shelf_depth,G

def main():
    out=ROOT/'visualizations';out.mkdir(exist_ok=True)
    cases=[];images=[]
    for name,length in [('Short shelf transition',200.),('Broad shelf transition',2400.)]:
        grid=Grid(-8000,8000,1000);times=np.linspace(0,450,61);amp=.2
        df=lambda x:shelf_depth(x,20.,5.,length)
        d=df(grid.x)
        e=lambda x:gaussian(x,amp,-2000.,150.)
        u=lambda x:np.sqrt(G*20)/20*e(x)
        li=cell_average(grid,lambda x:np.stack((e(x),u(x)),axis=-1))
        ni=cell_average(grid,lambda x:np.stack((df(x)+e(x),(df(x)+e(x))*u(x)),axis=-1))
        cfg=SolverConfig(cfl=.4,boundary='outflow')
        linear=simulate_linear_variable(grid,d,li[:,0],li[:,1],times,cfg)
        nonlinear=simulate_nonlinear_variable(grid,-d,ni[:,0],ni[:,1],times,cfg)
        assert np.isfinite(linear.eta).all() and np.isfinite(nonlinear.eta).all()
        cases.append({'name':name,'x':(grid.x[::2]/1000).round(4).tolist(),'depth':d[::2].round(5).tolist(),
            'time':times.tolist(),'linear':(linear.eta[:,::2]/amp).round(5).tolist(),
            'nonlinear':(nonlinear.eta[:,::2]/amp).round(5).tolist()})
        if length==2400:
            for i,t in enumerate(times):
                fig,(ax,bx)=plt.subplots(2,1,figsize=(8,5),gridspec_kw={'height_ratios':[3,1]},sharex=True)
                ax.plot(grid.x/1000,linear.eta[i]/amp,label='Linear',color='#176bb4')
                ax.plot(grid.x/1000,nonlinear.eta[i]/amp,label='Nonlinear',color='#c93d62')
                ax.set(xlim=(-8,8),ylim=(-.55,1.9),ylabel='Surface / initial amplitude',title=f'Broad shelf: t = {t:.1f} s')
                ax.legend(loc='upper right');ax.grid(alpha=.2)
                bx.fill_between(grid.x/1000,-d,-23,color='#a0a8b0');bx.plot(grid.x/1000,-d,color='#525c65')
                bx.set(ylim=(-23,0),ylabel='Bed (m)',xlabel='Position (km)');bx.grid(alpha=.2)
                fig.text(.5,.01,'Illustrative N=1000; D: 20 to 5 m; A=0.2 m; transition=2.4 km',ha='center',fontsize=9)
                fig.tight_layout(rect=(0,.035,1,1));buf=io.BytesIO();fig.savefig(buf,format='png',dpi=100);plt.close(fig)
                images.append(Image.open(buf).convert('RGB').copy())
    images[0].save(out/'shelf_propagation.gif',save_all=True,append_images=images[1:],duration=85,loop=0)
    template=(ROOT/'visualizations/explorer_template.html').read_text()
    (out/'wave_explorer.html').write_text(template.replace('__DATA__',json.dumps(cases,allow_nan=False,separators=(',',':'))))
    (out/'demo_metadata.json').write_text(json.dumps({'cells':1000,'amplitude_m':.2,'depths_m':[20,5],'sigma_m':150,'times_s':[0,450],'frames':61,'purpose':'Illustration only; numerical conclusions use results/*.csv'},indent=2))
    print('Offline explorer and 61-frame GIF built.')
if __name__=='__main__':main()
