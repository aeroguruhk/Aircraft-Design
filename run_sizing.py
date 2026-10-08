"""Run a complete conceptual sizing example. Works from CLI or Spyder.
    python run_sizing.py --example transport
    python run_sizing.py --config examples/general_aviation.json
"""
import argparse
import csv
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.optimize import brentq
import sizing_solver as constraints
from aircraft_weights import calculate_aircraft_weights
from Mission_module_SI import Requirements, weight_calc, mission_fuel_fraction
from perform_v3 import flight_curve, evaluate_altitude, ground_roll, cruise_capability, cruise_range_sweep
from SFD_BMD_v2 import sfd_bmd, plot_distributions

BASE=Path(__file__).resolve().parent
# Change this for a plain Run in Spyder.
DEFAULT_EXAMPLE='fighter'
CASE_TYPES={name:getattr(constraints,name) for name in (
    'LevelFlightCase','SustainedTurnCase','AccelerationCase','TakeoffCase',
    'LandingCase','ClimbCase','ServiceCeilingCase','TurnRateCase')}


def load_config(path):
    with open(path,encoding='utf-8') as f:
        config=json.load(f)
    required={'aircraft','mission','cases','ws_range'}
    if not required<=config.keys():
        raise ValueError(f'Unified JSON requires sections {sorted(required)}')
    p=constraints.AircraftParameters(**config['aircraft'])
    m=Requirements.from_dict(config['mission'])
    if p.aircraft_type not in ('fighter','transport','general_aviation','ga'):
        raise ValueError('aircraft_type must be fighter, transport or general_aviation')
    if p.propulsion_type not in ('jet','turboprop'):
        raise ValueError('propulsion_type must be jet or turboprop')
    for key in ('mass','AR','e','tc_root','L_t_m','L_fuse_m','D_fuse_m',
                'W_fuse_m','AR_ht','AR_vt','N_z','N_l','L_m_in','L_n_in','N_nw',
                'engine_tw','engine_pw','cl_max_clean','cl_max_landing','cl_max_to','M_max'):
        value=getattr(p,key,None)
        if value is None or not np.isfinite(value) or value<=0:
            raise ValueError(f'aircraft.{key} must be positive')
    if not 0<p.e<=1 or not 0<p.eta_prop<=1 or not 0<p.landing_weight_fraction<=1:
        raise ValueError('Efficiency and landing fraction must be in (0,1]')
    if not 0<=p.systems_fraction<1 or not 0<p.taper<=1:
        raise ValueError('Invalid systems fraction or taper')
    if p.propulsion_type=='turboprop' and (p.maximum_power_kw is None or p.maximum_power_kw<=0):
        raise ValueError('Turboprop requires positive maximum_power_kw')
    if p.propulsion_type=='jet' and p.maximum_thrust_N<=0:
        raise ValueError('Jet requires positive maximum_thrust_N')
    for key in ('sweep','sweep_ht','sweep_vt'):
        if abs(getattr(p,key,0))>=89:
            raise ValueError(f'{key} must have magnitude below 89 degrees')
    for key in ('tc_ht','tc_vt','taper_ht','taper_vt','K_y_m','K_z_m','S_fuse_wet_m2'):
        if hasattr(p,key) and getattr(p,key)<=0:
            raise ValueError(f'{key} must be positive')
    cases=[]
    for spec in config['cases']:
        args=dict(spec); typ=args.pop('type'); name=args.pop('name')
        if typ not in CASE_TYPES:
            raise ValueError(f'Unknown case type {typ}')
        case=CASE_TYPES[typ](name,**args)
        if case.beta<=0 or case.n<=0 or (case.mach is not None and case.mach<=0):
            raise ValueError(f'Invalid beta, load factor or Mach: {name}')
        constraints.Atmosphere.get_properties(case.altitude)
        if isinstance(case,(constraints.LevelFlightCase,constraints.SustainedTurnCase)) and case.mach is None:
            raise ValueError(f'{name} requires Mach')
        cases.append(case)
    if not cases or len({c.name for c in cases})!=len(cases):
        raise ValueError('Supply nonempty constraints with unique names')
    ws_cfg=config['ws_range']; start,end,step=(ws_cfg[k] for k in ('start_pa','end_pa','step_pa'))
    if not 0<start<end or step<=0:
        raise ValueError('Invalid wing loading range')
    ws=np.arange(start,end+step*.01,step)
    return config,p,m,cases,ws


def geometry(p,mass,ws):
    area=mass*p.g/ws; span=np.sqrt(area*p.AR)
    cr=2*area/(span*(1+p.taper))
    mac=2/3*cr*(1+p.taper+p.taper**2)/(1+p.taper)
    return {'area_m2':float(area),'span_m':float(span),'mac_m':float(mac),
            'ht_area_m2':float(p.V_h*area*mac/p.L_t_m),
            'vt_area_m2':float(p.V_v*area*span/p.L_t_m)}


def available_ratio(p,mass):
    return p.maximum_thrust_N/(mass*p.g) if p.propulsion_type=='jet' else p.maximum_power_kw/mass


def close_detailed_mass(p,m,ws,required,initial,fuel):
    """Bracket physical mass closure; do not label statistical fallback as convergence."""
    history=[]
    def components(mass):
        geo=geometry(p,mass,ws)
        ratio=available_ratio(p,mass) if getattr(p,'engine_weight_basis','required')=='installed' else required
        weights=calculate_aircraft_weights(p,m,mass,geo['area_m2'],geo['ht_area_m2'],geo['vt_area_m2'],ratio)
        if any(not np.isfinite(v) or v<0 for v in weights.values()):
            raise ValueError('Nonfinite or negative component weight; check correlation inputs')
        return weights
    def residual(mass):
        weights=components(mass)
        res=mass*(1-fuel)-sum(weights.values())-m.useful_mass
        history.append({'mass_kg':float(mass),'closure_residual_kg':float(res)})
        return res
    if getattr(p,'engine_weight_basis','required') not in ('required','installed'):
        raise ValueError('engine_weight_basis must be required or installed')
    mass_grid=np.geomspace(max(m.useful_mass,1),max(initial*30,m.useful_mass*40),100)
    lo=float(mass_grid[0]); flo=residual(lo)
    for hi in mass_grid[1:]:
        fhi=residual(float(hi))
        if flo<=0<=fhi:
            mass=float(brentq(residual,lo,float(hi),xtol=1e-7))
            return mass,components(mass),history
        lo,flo=float(hi),fhi
    raise ValueError('Detailed mass closure not bracketed. Revise mission/geometry/engine assumptions.')


def write_csv(path,rows):
    with open(path,'w',newline='',encoding='utf-8') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)


def performance_outputs(p,m,geo,required,fuel,out):
    alts=np.arange(0,15001,500); rows=[]
    fig,(ar,av)=plt.subplots(1,2,figsize=(12,5))
    reports={}
    for label,ratio in [('required',required),('available',available_ratio(p,p.mass))]:
        sweep=[evaluate_altitude(p,geo['area_m2'],float(h),ratio) for h in alts]
        for row in sweep:
            rows.append({'engine':label,**row})
        roc=[r['max_roc_m_s'] if r['max_roc_m_s'] is not None else np.nan for r in sweep]
        vmax=[r['max_speed_m_s'] if r['max_speed_m_s'] is not None else np.nan for r in sweep]
        ar.plot(alts/1000,roc,label=label);av.plot(alts/1000,vmax,label=label)
        # A ceiling requires an actual bracket; no extrapolation or endpoint substitution.
        ceiling=None
        for i in range(len(sweep)-1):
            r0,r1=roc[i:i+2]
            if np.isfinite(r0+r1) and r0>=.5>=r1:
                ceiling=float(brentq(lambda h:evaluate_altitude(p,geo['area_m2'],h,ratio)['max_roc_m_s']-.5,
                                     float(alts[i]),float(alts[i+1])))
                break
        h=p.cruise_altitude_m
        v=m.velocity/3.6
        c=flight_curve(p,geo['area_m2'],h,ratio,np.array([v]))
        reports[label]={'sea_level':sweep[0],'service_ceiling_m':ceiling,
                        'ground_roll':ground_roll(p,geo['area_m2'],ratio),
                        'cruise_roc_m_s':float(c['roc_m_s'][0]),
                        'cruise_cl':float(c['cl'][0]),
                        'cruise_ld_from_polar':float(p.mass*p.g/c['drag_N'][0]),
                        'cruise_level_feasible':bool(c['valid'][0] and c['roc_m_s'][0]>=0 and v<=p.M_max*constraints.Atmosphere.get_properties(h)[3])}
    ar.axhline(.5,color='gray',ls=':');ar.set_ylabel('Maximum climb rate (m/s)')
    av.plot(alts/1000,[r['stall_m_s'] for r in sweep],label='clean stall',ls='--');av.set_ylabel('Maximum level speed (m/s)')
    for ax in (ar,av):
        ax.set_xlabel('Altitude (km)');ax.legend();ax.grid(alpha=.3)
    fig.tight_layout();fig.savefig(out/'performance_altitude.png',dpi=150);plt.close(fig)
    write_csv(out/'performance.csv',rows)
    reports['fuel_allocation']=cruise_capability(p,m,fuel)
    fig,ax=plt.subplots(figsize=(8,5))
    for label,ratio in [('required',required),('available',available_ratio(p,p.mass))]:
        summary,range_rows=cruise_range_sweep(p,m,geo['area_m2'],ratio,fuel)
        reports[label]['cruise_range']=summary
        write_csv(out/f'cruise_range_{label}.csv',range_rows)
        ax.plot([r['speed_m_s'] for r in range_rows],
                [r['range_km'] if r['feasible'] else np.nan for r in range_rows],label=label)
    ax.axhline(m.Range,color='gray',ls=':',label='mission range')
    ax.set_xlabel('True airspeed (m/s)');ax.set_ylabel('Cruise range (km)')
    ax.set_title('Constant-altitude cruise; allocated fuel; drag polar')
    ax.legend();ax.grid(alpha=.3);fig.tight_layout()
    fig.savefig(out/'cruise_range.png',dpi=150);plt.close(fig)
    reports['glide_ld_max']=float(1/(2*np.sqrt(p.cd0_sub*p.k_sub)))
    return reports


def run(config_path,output_dir=None):
    config,p,m,cases,ws=load_config(config_path)
    out=Path(output_dir) if output_dir else BASE/'outputs'/Path(config_path).stem
    out.mkdir(parents=True,exist_ok=True)
    solver=constraints.SizingSolver(p,cases)
    curves=solver.solve(ws); boundary=solver.get_design_point(ws,curves)
    # Explicit clean-stall design limit, as distinct from an approach/go-around case.
    stall_limit=config.get('design',{}).get('stall_speed_limit_m_s')
    if stall_limit is not None:
        limit=.5*1.225*stall_limit**2*p.cl_max_clean
        boundary=np.where(ws<=limit,boundary,np.inf)
    feasible=np.flatnonzero(np.isfinite(boundary))
    if not len(feasible):
        raise ValueError('No wing loading satisfies all constraints and stall limit')
    specified=config.get('design',{}).get('wing_loading_pa')
    if specified is None:
        idx=int(feasible[np.argmin(boundary[feasible])])
    else:
        if not ws[0]<=specified<=ws[-1]:
            raise ValueError('Specified wing loading is outside grid')
        idx=int(np.argmin(abs(ws-specified)))
        if not np.isfinite(boundary[idx]):
            raise ValueError('Specified design wing loading is infeasible')
    design_ws,required=float(ws[idx]),float(boundary[idx])
    initial=weight_calc(m,2 if p.propulsion_type=='turboprop' else 1,p.eta_prop)
    fuel=mission_fuel_fraction(m,p.propulsion_type,p.eta_prop,p.g)
    mass,weights,history=close_detailed_mass(p,m,design_ws,required,initial,fuel)
    p.mass=mass;geo=geometry(p,mass,design_ws)
    available=available_ratio(p,mass)
    fig,_=solver.plot_diagram(ws,curves,np.where(np.isfinite(boundary),boundary,np.nan),design_ws,required,
                              available,save_path=out/'constraint_diagram.png');plt.close(fig)
    loads=sfd_bmd(weight=mass*p.g,span=geo['span_m'],wing_area=geo['area_m2'],taper=p.taper,
                  n=p.N_z,safety_factor=1,num_points=501)
    plot_distributions(loads,save_path=out/'wing_loads.png');plt.close('all')
    write_csv(out/'wing_loads.csv',[{'y_m':float(y),'shear_N':float(v),'moment_N_m':float(b)}
              for y,v,b in zip(loads['y'],loads['shear_force'],loads['bending_moment'])])
    write_csv(out/'mass_closure.csv',history)
    write_csv(out/'component_weights.csv',[{'component':k,'mass_kg':float(v)} for k,v in weights.items()])
    write_csv(out/'constraints.csv',[{'wing_loading_pa':float(x),**{k:float(v[i]) if np.isfinite(v[i]) else '' for k,v in curves.items()},
               'critical_ratio':float(boundary[i]) if np.isfinite(boundary[i]) else ''} for i,x in enumerate(ws)])
    perf=performance_outputs(p,m,geo,required,fuel,out)
    warnings=['Systems/nacelles/propellers/accessories use an aggregate systems allowance, not detailed Raymer equations.',
              'Schrenk loads are aerodynamic loads only; fuel, wing and engine inertia relief are omitted.',
              'Fixed fuselage dimensions do not resize with mass; volume and packaging closure are not checked.',
              'Takeoff is an approximate ground-run constraint; landing case is approach/go-around, not a landing-distance constraint.']
    if required>available:
        warnings.append('Available engine cannot satisfy the constraint design. Required engine is a conceptual resized engine.')
    if p.aircraft_type=='fighter':
        warnings.append('Fighter landing-gear equations retained from original code; verify against your source page.')
    if not perf['available']['cruise_level_feasible']:
        warnings.append('Prescribed cruise is infeasible at gross weight with available engine under the adopted drag model.')
    actual_range=perf['available']['cruise_range']['range_at_prescribed_speed_km']
    if actual_range is None or actual_range < m.Range*(1-1e-6):
        warnings.append('Available-engine cruise range at prescribed speed does not meet the mission under the drag polar; prescribed mission L/D is independent of that polar.')
    report={'title':config.get('title',Path(config_path).stem),'aircraft_type':p.aircraft_type,
            'propulsion_type':p.propulsion_type,'initial_statistical_mass_kg':initial,
            'gross_mass_kg':mass,'fuel_fraction':fuel,'fuel_mass_kg':mass*fuel,
            'useful_mass_kg':m.useful_mass,'component_masses_kg':weights,
            'mass_closure_residual_kg':mass-sum(weights.values())-mass*fuel-m.useful_mass,
            'wing_loading_pa':design_ws,'ratio_units':'T/W' if p.propulsion_type=='jet' else 'kW/kg',
            'required_engine_ratio':required,'available_engine_ratio':available,
            'available_engine_constraint_feasible':bool(available>=required),
            'engine_weight_basis':getattr(p,'engine_weight_basis','required'),
            'geometry':geo,'ultimate_load_factor':p.N_z,
            'root_shear_N':float(loads['shear_force'][-1]),
            'root_bending_N_m':float(loads['bending_moment'][-1]),'performance':perf,'warnings':warnings}
    with open(out/'summary.json','w') as f:
        json.dump(report,f,indent=2,allow_nan=False)
    print(f"{report['title']}: mass={mass:.2f} kg, S={geo['area_m2']:.2f} m², W/S={design_ws:.0f} Pa")
    print(f"Required {report['ratio_units']}={required:.4f}; available={available:.4f}; mass residual={report['mass_closure_residual_kg']:.3e} kg")
    for warning in warnings:
        print('NOTE:',warning)
    print('Results:',out.resolve())
    return report


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    group=parser.add_mutually_exclusive_group()
    group.add_argument('--example',choices=['fighter','transport','general_aviation','transport_turboprop'],default=None)
    group.add_argument('--config',type=Path)
    parser.add_argument('--output-dir',type=Path)
    args=parser.parse_args(argv)
    path=args.config or BASE/'examples'/f'{args.example or DEFAULT_EXAMPLE}.json'
    return run(path,args.output_dir)

if __name__=='__main__':
    main()
