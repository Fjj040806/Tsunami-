"""Fully wet variable-bed MUSCL/SSP-RK2 models.

Nonlinear hydrostatic reconstruction and centered source: Audusse et al.
(2004), equations 3.8--3.9. Bed slope is centered for smooth-bed accuracy;
select bed_reconstruction='constant' for aligned discontinuous steps.
No general positivity/entropy guarantee or wet/dry capability is claimed.
"""
from dataclasses import dataclass
import numpy as np
from .solver import SolverConfig, hll_flux


@dataclass
class VariableResult:
    x: np.ndarray
    bed: np.ndarray
    times: np.ndarray
    states: np.ndarray
    model: str
    diagnostics: dict

    @property
    def eta(self):
        return self.states[:, :, 0] if self.model == 'linear' else self.states[:, :, 0]+self.bed

    @property
    def velocity(self):
        return self.states[:, :, 1] if self.model == 'linear' else self.states[:, :, 1]/self.states[:, :, 0]


def _reconstruct(values, boundary, bed_reconstruction, right_wall=False):
    # Columns: surface, velocity, bed. Return each cell's left/right traces,
    # including one ghost cell at each end.
    if boundary == 'periodic':
        ext=np.concatenate((values[-2:],values,values[:2]))
    elif boundary == 'outflow':
        ext=np.pad(values,((2,2),(0,0)),mode='edge')
        if right_wall:
            ext[-2:]=values[-2:][::-1]
            ext[-2:,1]*=-1
    else:
        ext=np.concatenate((values[:2][::-1],values,values[-2:][::-1]))
        ext[:2,1]*=-1;ext[-2:,1]*=-1
    center=ext[1:-1];dm=center-ext[:-2];dp=ext[2:]-center
    middle=.5*(dm+dp)
    slope=np.where(dm*dp>0,np.sign(middle)*np.minimum.reduce(
        (2*np.abs(dm),np.abs(middle),2*np.abs(dp))),0.)
    slope[:,2]=middle[:,2] if bed_reconstruction=='centered' else 0.
    return center-.5*slope,center+.5*slope


def linear_interface_flux(left, right, g):
    """Exact stationary-interface flux with eta and D*u continuous.

    Traces contain eta,u,b with D=-b > 0. Outgoing characteristics on
    each side determine the common surface and discharge at the interface.
    """
    dl,dr=-left[:,2],-right[:,2]
    cl,cr=np.sqrt(g*dl),np.sqrt(g*dr)
    ql,qr=dl*left[:,1],dr*right[:,1]
    surface=(cl*left[:,0]+cr*right[:,0]+ql-qr)/(cl+cr)
    discharge=(cr*ql+cl*qr+cl*cr*(left[:,0]-right[:,0]))/(cl+cr)
    return np.column_stack((discharge,g*surface)),np.maximum(cl,cr)


def _rhs(state, bed, cfg, model, bed_reconstruction, boundary=None, t=0., step_coupling="hydrostatic"):
    if model=='linear':
        values=np.column_stack((state[:,0],state[:,1],bed))
    else:
        values=np.column_stack((state[:,0]+bed,state[:,1]/state[:,0],bed))
    lc,rc=_reconstruct(values,cfg.boundary,bed_reconstruction, boundary is not None and boundary.right_wall)
    left,right=rc[:-1].copy(),lc[1:].copy()
    if boundary is not None:
        eta_in=float(boundary.incoming_surface(t))
        if not np.isfinite(eta_in):raise ValueError("Incoming surface must be finite")
        left[0,2]=right[0,2]
        d=-right[0,2]
        if d<=cfg.minimum_depth or d+eta_in<=cfg.minimum_depth:raise ValueError("Incoming depth must be fully wet")
        if model=="linear":
            left[0,0]=eta_in;left[0,1]=np.sqrt(cfg.g/d)*eta_in
        else:
            h_int=right[0,0]-right[0,2]
            if h_int<=cfg.minimum_depth:raise ValueError("Boundary interior is near dry")
            rp=4*np.sqrt(cfg.g*(d+eta_in))-2*np.sqrt(cfg.g*d)
            rm=right[0,1]-2*np.sqrt(cfg.g*h_int)
            c=.25*(rp-rm);u=.5*(rp+rm)
            if c<=0 or abs(u)>=c:raise ValueError("Boundary requires positive subcritical state")
            left[0,0]=c*c/cfg.g+left[0,2];left[0,1]=u
    if model=='linear':
        if np.any(-lc[:,2]<=cfg.minimum_depth) or np.any(-rc[:,2]<=cfg.minimum_depth):
            raise ValueError('Reconstructed still-water depth is not fully wet')
        if np.any(lc[:,0]-lc[:,2]<=cfg.minimum_depth) or np.any(rc[:,0]-rc[:,2]<=cfg.minimum_depth):
            raise ValueError('Linear surface violates fully-wet scope')
        flux,speeds=linear_interface_flux(left,right,cfg.g)
        rhs=-(flux[1:]-flux[:-1])
        speed=max(float(speeds.max()),float(np.sqrt(cfg.g*(-bed)).max()))
    else:
        hl,hr=left[:,0]-left[:,2],right[:,0]-right[:,2]
        if np.any(hl<=cfg.minimum_depth) or np.any(hr<=cfg.minimum_depth):
            raise ValueError('Reconstructed physical depth is not fully wet')
        bstar=np.maximum(left[:,2],right[:,2])
        hs_l,hs_r=left[:,0]-bstar,right[:,0]-bstar
        # Intentionally restrict even hydrostatic traces to fully wet states.
        if np.any(hs_l<=cfg.minimum_depth) or np.any(hs_r<=cfg.minimum_depth):
            raise ValueError('Hydrostatic interface became near dry; no clipping')
        flux,speeds=hll_flux(np.column_stack((hs_l,hs_l*left[:,1])),
                            np.column_stack((hs_r,hs_r*right[:,1])),cfg.g)
        fl,fr=flux.copy(),flux.copy()
        fl[:,1]+=.5*cfg.g*(hl**2-hs_l**2)
        fr[:,1]+=.5*cfg.g*(hr**2-hs_r**2)
        if step_coupling=='bernoulli':
            mask=np.abs(left[:,2]-right[:,2])>1e-12
            if np.any(mask):
                coupled_l,coupled_r,coupled_speed=stationary_step_flux(left[mask],right[mask],cfg.g)
                fl[mask]=coupled_l;fr[mask]=coupled_r
                flux[mask,0]=coupled_l[:,0]
                speeds[mask]=np.maximum(speeds[mask],coupled_speed)
        rhs=-(fl[1:]-fr[:-1])
        # Centered intra-cell bed source, integrated over each physical cell.
        hcell_l=lc[1:-1,0]-lc[1:-1,2]
        hcell_r=rc[1:-1,0]-rc[1:-1,2]
        rhs[:,1]+=.5*cfg.g*(hcell_l+hcell_r)*(lc[1:-1,2]-rc[1:-1,2])
        speed=max(float(speeds.max()),float(np.max(np.abs(state[:,1]/state[:,0])+np.sqrt(cfg.g*state[:,0]))),
                  float(np.max(np.abs(left[:,1])+np.sqrt(cfg.g*hl))),
                  float(np.max(np.abs(right[:,1])+np.sqrt(cfg.g*hr))))
    return rhs, float(flux[-1,0]-flux[0,0]), speed


def _integrate(grid, bed, initial, times, model, config, bed_reconstruction, forcing, boundary, step_coupling):
    cfg=config or SolverConfig(cfl=.4,boundary='periodic')
    if step_coupling not in ('hydrostatic','bernoulli'):
        raise ValueError('Unknown step coupling')
    if boundary is not None and cfg.boundary!='outflow':
        raise ValueError('IncidentBoundary requires config.boundary=outflow')
    if bed_reconstruction not in ('centered','constant'):
        raise ValueError('bed_reconstruction must be centered or constant')
    bed=np.asarray(bed,dtype=float).copy();state=np.asarray(initial,dtype=float).copy()
    if bed.shape!=(grid.cells,) or state.shape!=(grid.cells,2):
        raise ValueError('Invalid initial array shape')
    if not np.isfinite(bed).all() or not np.isfinite(state).all():
        raise ValueError('Initial arrays must be finite')
    if model=='linear' and np.any(-bed<=cfg.minimum_depth):
        raise ValueError('Linearization about eta=0 needs D=-bed > minimum_depth')
    times=np.asarray(times,dtype=float)
    if times.ndim!=1 or len(times)<1 or not np.isfinite(times).all() or times[0]!=0 or np.any(np.diff(times)<=0):
        raise ValueError('Output times must start at zero and increase strictly')
    def depth(w):return w[:,0]-bed if model=='linear' else w[:,0]
    def check(w):
        if not np.isfinite(w).all() or np.any(depth(w)<=cfg.minimum_depth):
            raise ValueError('Fully-wet state violated; no evolved-depth clipping')
    def mass(w):return float(np.sum(depth(w))*grid.dx)
    def evaluate(w,t):
        check(w)
        rhs,net,speed=_rhs(w,bed,cfg,model,bed_reconstruction,boundary,t,step_coupling)
        rhs=rhs/grid.dx
        source_mass=0.
        if forcing is not None:
            source=np.asarray(forcing(t),dtype=float)
            if source.shape!=w.shape or not np.isfinite(source).all():
                raise ValueError('Forcing must provide finite cell-average array (N,2)')
            rhs+=source;source_mass=float(np.sum(source[:,0])*grid.dx)
        return rhs,net,source_mass,speed
    check(state);mass0=mass(state);states=[state.copy()]
    t=0.;steps=0;rejected=0;boundary_total=0.;source_total=0.;max_res=0.;max_cfl=0.
    minimum=float(depth(state).min());residuals=[0.];boundaries=[0.];sources=[0.]
    for target in times[1:]:
        while t<target:
            if steps>=cfg.max_steps:raise RuntimeError('Maximum steps exceeded')
            r0,n0,m0,s0=evaluate(state,t)
            dt=min(cfg.cfl*grid.dx/s0,float(target-t))
            for attempt in range(20):
                if dt<=0 or t+dt==t:raise RuntimeError('Time-step underflow')
                stage=state+dt*r0
                r1,n1,m1,s1=evaluate(stage,t+dt)
                if dt*s1/grid.dx>cfg.cfl*(1+1e-12):
                    dt=min(.9*dt,cfg.cfl*grid.dx/s1);rejected+=1;continue
                updated=.5*state+.5*(stage+dt*r1);check(updated)
                break
            else:raise RuntimeError('RK-stage CFL adjustment failed')
            state=updated;boundary_total+=.5*dt*(n0+n1);source_total+=.5*dt*(m0+m1)
            t=float(target) if dt==target-t else t+dt;steps+=1
            max_res=max(max_res,abs(mass(state)-mass0+boundary_total-source_total))
            max_cfl=max(max_cfl,dt*max(s0,s1)/grid.dx);minimum=min(minimum,float(depth(state).min()))
        states.append(state.copy());residuals.append(mass(state)-mass0+boundary_total-source_total)
        boundaries.append(boundary_total);sources.append(source_total)
    diagnostics={'model':model,'order':2,'bed_reconstruction':bed_reconstruction,'step_coupling':step_coupling,
        'boundary':cfg.boundary,'steps':steps,'rejected_cfl_trials':rejected,
        'maximum_actual_cfl':max_cfl,'minimum_depth_m':minimum,'depth_clipping_count':0,
        'maximum_normalized_mass_balance_residual':max_res/abs(mass0),
        'mass_balance_residual_per_output_m2':residuals,
        'boundary_integral_per_output_m2':boundaries,'forcing_mass_integral_per_output_m2':sources,
        'manufactured_forcing_enabled':forcing is not None,
        'incident_boundary_enabled':boundary is not None}
    return VariableResult(grid.x,bed,times.copy(),np.asarray(states),model,diagnostics)


def simulate_nonlinear_variable(grid, bed, h0, q0, times, config=None,
                                bed_reconstruction='centered', forcing=None, boundary=None, step_coupling='hydrostatic'):
    """Nonlinear cell averages (h,q), static bed cell averages.

    Optional forcing(t)->(N,2) cell averages is for manufactured verification.
    """
    return _integrate(grid,bed,np.column_stack((h0,q0)),times,'nonlinear',config,bed_reconstruction,forcing,boundary,step_coupling)


def simulate_linear_variable(grid, depth, eta0, u0, times, config=None,
                             bed_reconstruction='centered', forcing=None, boundary=None, step_coupling='hydrostatic'):
    """Linear cell averages (eta,u); depth contains positive D(x) averages."""
    return _integrate(grid,-np.asarray(depth),np.column_stack((eta0,u0)),times,'linear',config,bed_reconstruction,forcing,boundary,step_coupling)


@dataclass(frozen=True)
class IncidentBoundary:
    """Incoming right-going surface signal at flat left edge; optional right wall.

    Measured total surface can serve as incident input only before reflected
    waves reach that input gauge. Nonlinear input uses a simple-wave invariant.
    """
    incoming_surface: object
    right_wall: bool = True


def stationary_step_flux(left,right,g):
    """Subcritical stationary step with discharge and Bernoulli-head matching.

    Exact homogeneous wave curves connect each trace to its interface state.
    This is an inviscid energy-conserving junction model, not a guarantee for
    breaking flow, hydraulic jumps or real vertical-step losses.
    """
    hl0=left[:,0]-left[:,2];hr0=right[:,0]-right[:,2]
    ul0,ur0=left[:,1],right[:,1]
    if np.any(np.abs(ul0)>=np.sqrt(g*hl0)) or np.any(np.abs(ur0)>=np.sqrt(g*hr0)):
        raise ValueError('Stationary coupling supports only subcritical traces')
    cl,cr=np.sqrt(g*hl0),np.sqrt(g*hr0)
    eta=(cl*left[:,0]+cr*right[:,0]+hl0*ul0-hr0*ur0)/(cl+cr)
    hl=eta-left[:,2];hr=eta-right[:,2]
    def wave(h,h0):
        s=np.sqrt(.5*g*(1/h+1/h0))
        shock=(h-h0)*s
        ds=s-(h-h0)*.25*g/(h*h*s)
        return np.where(h>h0,shock,2*(np.sqrt(g*h)-np.sqrt(g*h0))),np.where(h>h0,ds,np.sqrt(g/h))
    for iteration in range(30):
        if np.any(hl<=0) or np.any(hr<=0):raise ValueError('Nonpositive junction state')
        f1,d1=wave(hl,hl0);f2,d2=wave(hr,hr0)
        ul=ul0-f1;ur=ur0+f2
        rq=hl*ul-hr*ur
        rh=hl+left[:,2]+ul*ul/(2*g)-hr-right[:,2]-ur*ur/(2*g)
        scale=np.maximum(hl0,hr0)
        if np.max(np.abs(rq)/(scale*np.sqrt(g*scale))+np.abs(rh)/scale)<1e-13:break
        a=ul-hl*d1;b=-ur-hr*d2;c=1-ul*d1/g;d=-1-ur*d2/g
        determinant=a*d-b*c
        if np.any(np.abs(determinant)<1e-14):raise RuntimeError('Singular stationary junction')
        dh_l=(d*rq-b*rh)/determinant;dh_r=(-c*rq+a*rh)/determinant
        damping=np.minimum(1.,np.minimum(np.where(dh_l>0,.8*hl/np.maximum(dh_l,1e-300),1.),np.where(dh_r>0,.8*hr/np.maximum(dh_r,1e-300),1.)))
        hl-=damping*dh_l;hr-=damping*dh_r
    else:raise RuntimeError('Stationary junction Newton solve did not converge')
    if np.any(np.abs(ul)>=np.sqrt(g*hl)) or np.any(np.abs(ur)>=np.sqrt(g*hr)):
        raise ValueError('Junction solution became supercritical')
    q=.5*(hl*ul+hr*ur)
    fl=np.column_stack((q,q*q/hl+.5*g*hl*hl))
    fr=np.column_stack((q,q*q/hr+.5*g*hr*hr))
    speed=np.maximum(np.abs(ul)+np.sqrt(g*hl),np.abs(ur)+np.sqrt(g*hr))
    return fl,fr,speed
