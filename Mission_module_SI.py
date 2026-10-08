"""SI mission sizing. Jet TSFC: 1/hour; turboprop BSFC: kg/(kW hour)."""
from dataclasses import dataclass, fields
import json
import numpy as np
from scipy.optimize import brentq

@dataclass
class Requirements:
    Pax: float = 1
    Pax_mass: float = 100
    payload: float = 4000
    Range: float = 500                 # km
    Endurance: float = .5              # loiter hours, additional to cruise
    Engine: int = 1
    sfc_cruise: float = .6             # jet only, 1/hour
    sfc_loiter: float = .5             # jet only, 1/hour
    bsfc_cruise: float = .30           # turboprop only, kg/(kW hour)
    bsfc_loiter: float = .30
    loiter_velocity_kmh: float = 200
    velocity: float = 600              # true airspeed, km/hour
    LBYD: float = 14
    cruise_ld_factor: float = .866
    reserve_multiplier: float = 1.06
    takeoff_remaining: float = .97
    climb_remaining: float = .985
    landing_remaining: float = .995
    empty_fraction_reference: float = .55  # teaching initial-guess fit
    reference_mass_kg: float = 10000
    empty_fraction_exponent: float = -.06

    @classmethod
    def from_dict(cls, data):
        unknown = set(data)-{f.name for f in fields(cls)}
        if unknown:
            raise ValueError(f'Unknown mission keys: {sorted(unknown)}')
        return cls(**data)

    @classmethod
    def from_json(cls, filepath):
        with open(filepath) as f:
            config=json.load(f)
        return cls.from_dict(config.get('mission',config.get('data',config)))

    @property
    def useful_mass(self):
        return self.payload + self.Pax*self.Pax_mass


def mission_fuel_fraction(ac, propulsion_type='jet', eta_prop=.8, g=9.81):
    if ac.Range < 0 or ac.Endurance < 0 or ac.velocity <= 0 or ac.LBYD <= 0:
        raise ValueError('Mission requires nonnegative range/loiter and positive speed/L/D')
    if ac.payload < 0 or ac.Pax < 0 or ac.Pax_mass < 0 or ac.useful_mass <= 0:
        raise ValueError('Mission useful mass must be positive')
    if not 0 < ac.cruise_ld_factor <= 1 or ac.reserve_multiplier < 1:
        raise ValueError('Invalid cruise_ld_factor or reserve_multiplier')
    if any(not 0 < x <= 1 for x in (ac.takeoff_remaining,ac.climb_remaining,ac.landing_remaining)):
        raise ValueError('Segment remaining fractions must be in (0,1]')
    ld=ac.LBYD*ac.cruise_ld_factor
    if propulsion_type=='jet':
        if min(ac.sfc_cruise,ac.sfc_loiter)<=0:
            raise ValueError('Jet TSFC must be positive, in 1/hour')
        cruise=np.exp(-ac.Range*ac.sfc_cruise/(ac.velocity*ld))
        loiter=np.exp(-ac.Endurance*ac.sfc_loiter/ac.LBYD)
    elif propulsion_type=='turboprop':
        if not 0 < eta_prop <= 1 or min(ac.bsfc_cruise,ac.bsfc_loiter,ac.loiter_velocity_kmh)<=0:
            raise ValueError('Turboprop requires positive BSFC, loiter speed and efficiency')
        # Fuel mass rate = BSFC * shaft power. BSFC conversion to kg/J.
        cruise=np.exp(-ac.Range*1000*g*ac.bsfc_cruise/(3.6e6*eta_prop*ld))
        loiter=np.exp(-ac.Endurance*g*ac.bsfc_loiter*(ac.loiter_velocity_kmh/3.6)
                      /(1000*eta_prop*ac.LBYD))
    else:
        raise ValueError('Unknown propulsion type')
    fraction=ac.reserve_multiplier*(1-ac.takeoff_remaining*ac.climb_remaining*cruise*loiter*ac.landing_remaining)
    if not 0 <= fraction < 1:
        raise ValueError('Mission fuel fraction is outside [0,1)')
    ac.fuel_fraction=float(fraction)
    return float(fraction)


def empty_weight_fraction(W0, ac=None):
    ac=ac or Requirements()
    return ac.empty_fraction_reference*(W0/ac.reference_mass_kg)**ac.empty_fraction_exponent


def weight_calc(ac, etype=1, eta_prop=.8):
    fuel=mission_fuel_fraction(ac,'jet' if etype==1 else 'turboprop',eta_prop)
    residual=lambda mass: mass*(1-empty_weight_fraction(mass,ac)-fuel)-ac.useful_mass
    masses=np.geomspace(max(ac.useful_mass,1),1e7,400)
    for lo,hi in zip(masses[:-1],masses[1:]):
        if residual(lo)<=0 and residual(hi)>=0:
            return float(brentq(residual,lo,hi))
    raise ValueError('No statistical mass closure below 10 million kg; revise mission/empty-fraction fit')


def range_weight(R,c,v,LD,etype=1):
    """Legacy jet helper: R metres, c 1/s, v m/s. Propeller needs BSFC instead."""
    if etype != 1:
        raise ValueError('Use mission_fuel_fraction for a turboprop mission')
    return np.exp(-R*c/(v*LD))


def loiter_weight(E,c,LD,etype=1):
    if etype != 1:
        raise ValueError('Use mission_fuel_fraction for a turboprop mission')
    return np.exp(-E*c/LD)
