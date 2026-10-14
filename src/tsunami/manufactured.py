"""Analytic forced test fields; not physical tsunami observations."""
import numpy as np
from .initial import cell_average
L=1000.;D0=10.;B=2.;A=.2;G=9.81;C=np.sqrt(G*D0)
def fields(x,t):
    k=2*np.pi/L;phase=k*(x-C*t)
    d=D0+B*np.cos(k*x);e=A*np.sin(phase);u=C/D0*e
    return d,e,u
def exact_state(x,t,model):
    d,e,u=fields(x,t)
    return np.stack((e,u),axis=-1) if model=='linear' else np.stack((d+e,(d+e)*u),axis=-1)
def analytic_forcing(x,t,model):
    d,e,u=fields(x,t);h=d+e;k=2*np.pi/L
    dx=-B*k*np.sin(k*x);ex=A*k*np.cos(k*(x-C*t));et=-C*ex
    ux=C/D0*ex;ut=C/D0*et
    if model=='linear':return np.stack((et+dx*u+d*ux,ut+G*ex),axis=-1)
    return np.stack((et+(dx+ex)*u+h*ux,et*u+h*ut+(dx+ex)*u*u+2*h*u*ux+G*h*ex),axis=-1)
def setup(grid,model):
    d=cell_average(grid,lambda x:fields(x,0)[0])
    z=cell_average(grid,lambda x:exact_state(x,0,model))
    return d,z,lambda t:cell_average(grid,lambda x:analytic_forcing(x,t,model))
