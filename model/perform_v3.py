"""Performance with one drag/lapse model shared with the constraint module.
All arrays use SI. Parabolic drag is a teaching approximation; no trim/wave-drag build-up.
"""
import numpy as np
from scipy.optimize import brentq
from sizing_solver import Atmosphere


def flight_curve(p, area, altitude, engine_ratio, speeds, mass_fraction=1):
    _,sigma,rho,a=Atmosphere.get_properties(altitude)
    weight=p.mass*p.g*mass_fraction
    q=.5*rho*speeds**2
    cl=weight/(q*area)
    sup=speeds/a>=1
    cd0=np.where(sup,p.cd0_sup,p.cd0_sub)
    k=np.where(sup,p.k_sup,p.k_sub)
    drag=q*area*(cd0+k*cl**2)
    lapse=sigma**.7
    if p.propulsion_type=='jet':
        thrust=np.full_like(speeds,engine_ratio*p.mass*p.g*lapse)
    else:
        thrust=engine_ratio*p.mass*1000*p.eta_prop*lapse/speeds
    roc=(thrust-drag)*speeds/weight
    clmax=getattr(p,'cl_max_clean',1.5)
    return {'speed_m_s':speeds,'drag_N':drag,'thrust_N':thrust,
            'power_required_W':drag*speeds,'power_available_W':thrust*speeds,
            'roc_m_s':roc,'cl':cl,'valid':cl<=clmax,'mach':speeds/a}


def evaluate_altitude(p,area,altitude,engine_ratio):
    _,_,rho,a=Atmosphere.get_properties(altitude)
    stall=np.sqrt(2*p.mass*p.g/(rho*area*p.cl_max_clean))
    max_speed=p.M_max*a
    if stall>=max_speed:
        return {'altitude_m':altitude,'stall_m_s':float(stall),'max_speed_m_s':None,
                'max_roc_m_s':None,'speed_at_max_roc_m_s':None,'speed_cap_limited':False}
    speeds=np.linspace(stall,max_speed,1200)
    curve=flight_curve(p,area,altitude,engine_ratio,speeds)
    roc=curve['roc_m_s']; best=int(np.argmax(roc)); good=np.flatnonzero(roc>=0)
    vmax=None; capped=False
    if len(good):
        j=int(good[-1]); vmax=float(speeds[j]); capped=j==len(speeds)-1
        if not capped:
            vmax=float(brentq(lambda v:float(flight_curve(p,area,altitude,engine_ratio,np.array([v]))['roc_m_s'][0]),speeds[j],speeds[j+1]))
    return {'altitude_m':float(altitude),'stall_m_s':float(stall),'max_speed_m_s':vmax,
            'max_roc_m_s':float(roc[best]),'speed_at_max_roc_m_s':float(speeds[best]),
            'speed_cap_limited':bool(capped)}


def ground_roll(p,area,engine_ratio):
    """Representative-speed ground-roll estimate; excludes obstacles, flare and air distance."""
    rho=1.225; w=p.mass*p.g; cl_to=p.cl_max_to; cl_l=p.cl_max_landing
    vt=.7*1.2*np.sqrt(2*w/(rho*area*cl_to))
    thrust=engine_ratio*w if p.propulsion_type=='jet' else engine_ratio*p.mass*1000*p.eta_prop/vt
    q=.5*rho*vt**2
    lift=q*area*cl_to
    drag=q*area*(p.cd0_sub+p.k_sub*cl_to**2)
    force=thrust-drag-.02*max(w-lift,0)
    to=1.44*w*w/(p.g*rho*cl_to*area*force) if force>0 else None
    wl=w*p.landing_weight_fraction
    vl=.7*1.3*np.sqrt(2*wl/(rho*area*cl_l))
    ql=.5*rho*vl**2
    retarding=ql*area*(p.cd0_sub+p.k_sub*cl_l**2)+.4*max(wl-ql*area*cl_l,0)
    land=1.69*wl*wl/(p.g*rho*cl_l*area*retarding)
    return {'takeoff_ground_roll_m':None if to is None else float(to),'landing_ground_roll_m':float(land)}


def cruise_capability(p,mission,fuel_fraction):
    """Breguet at prescribed constant cruise speed/L/D; excludes reserved fuel.
    Range/endurance alternatives use the same usable cruise fuel, not both together.
    """
    fuel_before_reserve=fuel_fraction/mission.reserve_multiplier
    cruise_remaining=(1-fuel_before_reserve)/(mission.takeoff_remaining*mission.climb_remaining
                      *mission.landing_remaining)
    # Recover cruise allocation excluding the planned loiter consumption.
    if p.propulsion_type=='jet':
        loiter=np.exp(-mission.Endurance*mission.sfc_loiter/mission.LBYD)
    else:
        loiter=np.exp(-mission.Endurance*p.g*mission.bsfc_loiter*(mission.loiter_velocity_kmh/3.6)
                      /(1000*p.eta_prop*mission.LBYD))
    cruise_remaining=min(1,cruise_remaining/loiter)
    logratio=-np.log(cruise_remaining)
    ld=mission.LBYD*mission.cruise_ld_factor
    if p.propulsion_type=='jet':
        range_km=mission.velocity*ld/mission.sfc_cruise*logratio
        endurance=ld/mission.sfc_cruise*logratio
    else:
        range_km=3.6e6*p.eta_prop*ld/(p.g*mission.bsfc_cruise)*logratio/1000
        endurance=range_km/mission.velocity
    return {'allocated_cruise_range_km':float(range_km),
            'allocated_cruise_endurance_hr':float(endurance),
            'note':'Fuel allocation check at prescribed L/D; aerodynamic/engine feasibility reported separately.'}


def cruise_range_sweep(p,mission,area,engine_ratio,fuel_fraction):
    """Integrate mass burn at constant altitude/speed using the configured drag polar.
    Cruise fuel excludes planned loiter, landing and reserve. Installed engine is
    checked at the heaviest cruise mass; no engine-efficiency/TSFC variation model.
    """
    _,sigma,rho,a=Atmosphere.get_properties(p.cruise_altitude_m)
    mi=p.mass*mission.takeoff_remaining*mission.climb_remaining
    allocation=cruise_capability(p,mission,fuel_fraction)
    if p.propulsion_type=='jet':
        logratio=allocation['allocated_cruise_endurance_hr']*mission.sfc_cruise/(mission.LBYD*mission.cruise_ld_factor)
    else:
        logratio=allocation['allocated_cruise_range_km']*1000*p.g*mission.bsfc_cruise/(3.6e6*p.eta_prop*mission.LBYD*mission.cruise_ld_factor)
    mf=mi*np.exp(-logratio)
    masses=np.linspace(mf,mi,301)
    stall=np.sqrt(2*mi*p.g/(rho*area*p.cl_max_clean))
    speeds=np.unique(np.r_[np.linspace(stall,p.M_max*a,500),mission.velocity/3.6])
    q=.5*rho*speeds[:,None]**2
    cl=masses[None,:]*p.g/(q*area)
    sup=speeds[:,None]/a>=1
    cd0=np.where(sup,p.cd0_sup,p.cd0_sub)
    k=np.where(sup,p.k_sup,p.k_sub)
    drag=q*area*(cd0+k*cl**2)
    if p.propulsion_type=='jet':
        # TSFC 1/hr: fuel mass rate = TSFC * D/g, kg/hr.
        dt_dm=p.g/(mission.sfc_cruise*drag)
        available=engine_ratio*p.mass*p.g*sigma**.7
        can_power=np.max(drag,axis=1)<=available
    else:
        shaft_kw=drag*speeds[:,None]/(1000*p.eta_prop)
        dt_dm=1/(mission.bsfc_cruise*shaft_kw)
        can_power=np.max(shaft_kw,axis=1)<=engine_ratio*p.mass*sigma**.7
    trap=getattr(np,'trapezoid',np.trapz if hasattr(np,'trapz') else None)
    endurance=trap(dt_dm,masses,axis=1)
    ranges=endurance*speeds*3.6
    valid=can_power&(np.max(cl,axis=1)<=p.cl_max_clean)&(speeds<=p.M_max*a)
    j=int(np.argmin(abs(speeds-mission.velocity/3.6)))
    eligible=np.flatnonzero(valid)
    best=int(eligible[np.argmax(ranges[eligible])]) if len(eligible) else None
    longest=int(eligible[np.argmax(endurance[eligible])]) if len(eligible) else None
    summary={'cruise_start_mass_kg':float(mi),'cruise_end_mass_kg':float(mf),
             'range_at_prescribed_speed_km':float(ranges[j]) if valid[j] else None,
             'prescribed_speed_feasible':bool(valid[j]),
             'maximum_cruise_range_km':float(ranges[best]) if best is not None else None,
             'range_optimum_speed_m_s':float(speeds[best]) if best is not None else None,
             'maximum_cruise_endurance_hr':float(endurance[longest]) if longest is not None else None,
             'endurance_optimum_speed_m_s':float(speeds[longest]) if longest is not None else None}
    rows=[{'speed_m_s':float(v),'range_km':float(r),'endurance_hr':float(e),'feasible':bool(ok)}
          for v,r,e,ok in zip(speeds,ranges,endurance,valid)]
    return summary,rows
