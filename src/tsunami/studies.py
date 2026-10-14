"""Physical study definitions and independent analytical references."""
import numpy as np
G=9.81

def gaussian(x,amplitude,center,sigma):
    return amplitude*np.exp(-.5*((x-center)/sigma)**2)

def shelf_depth(x,deep,shallow,length):
    # C2 transition: cubic derivative vanishes at both plateau endpoints.
    s=np.clip(x/length+.5,0,1)
    smooth=s*s*s*(10+s*(-15+6*s))
    return deep+(shallow-deep)*smooth

def step_reference(x,t,deep,shallow,amplitude,center,sigma):
    cl,cr=np.sqrt(G*deep),np.sqrt(G*shallow)
    refl=(cl-cr)/(cl+cr);trans=2*cl/(cl+cr)
    inc=gaussian(x-cl*t,amplitude,center,sigma)
    reflected=refl*gaussian(-x-cl*t,amplitude,center,sigma)
    transmitted=trans*gaussian(cl/cr*x-cl*t,amplitude,center,sigma)
    eta=np.where(x<0,inc+reflected,transmitted)
    u=np.where(x<0,cl/deep*(inc-reflected),cr/shallow*transmitted)
    return np.stack((eta,u),axis=-1)

def composite_depth(x):
    # NOAA figure: shore-origin junctions 8.181, 3.821, .891 m;
    # convert to wavemaker-origin coordinates using 23.22 m wall position.
    x=np.asarray(x)
    return (.218-np.clip(x-15.039,0,19.399-15.039)/53
            -np.clip(x-19.399,0,22.329-19.399)/150
            -np.maximum(x-22.329,0)/13)
