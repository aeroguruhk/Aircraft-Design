# -*- coding: utf-8 -*-
"""
Created on Thu Jul  9 11:31:51 2026

@author: hknar
"""
import numpy as np
from dataclasses import dataclass

@dataclass 
 
class Aircraft:
     W_dg: float = 0.0
     S_w: float = 0.0
     AR: float = 0.0
     taper: float = 0.0
     sweep: float = 0.0
     tc_root: float = 0.0
     F_w : float = 0.0   
     B_h : float = 0.0 # Horizontail span
     S_ht : float = 0.0 # Horizontail area
     H_t : float = 0.0  # Horizontail height
     H_v : float = 0.0  # Vertical Tail height
     S_vt : float = 0.0 # Vertical Tail area
     M : float = 0.0    #  Mach number
     L_t : float = 0.0 # Vertical Tail arm
     S_r :float = 0.0 # Vertical Tail rudder area
     A_vt :float = 0.0 # Vertical Tail aspect ratio
     sweep_vt : float = 0.0 # Vertical Tail sweep
     L : float = 0.0  # fuselage length
     D : float = 0.0  # fuselage depth
     W : float = 0.0  # fuselage widt

    
''' 
Note  : all equations are in US units, so we input in US units get the weights in lbs 
 and then convert back to SI
'''

def wing_weight_fighter(ac, K_dw, K_vs, N_z, S_csw):
    """
    Raymer Eq. xx.x

    Inputs
    ------
    W_dg : lb
    S_w : ft²
    AR : -
    tc_root : -
    Sweep : radians
    N_z : ultimate load factor

    Returns
    -------
    Wing structural weight for a combat aircraft(lb)
    """
    W_wing_F = (
        0.0103
        * K_dw
        * K_vs
        * (ac.W_dg * N_z)**0.5
        * ac.S_w**0.622
        * ac.AR**0.785
        * ac.tc_root**(-0.4)
        * (1 + ac.taper)**0.05
        * np.cos(ac.sweep)**(-1.0)
        * S_csw**0.04
    )
    return W_wing_F
def wing_weight_transport(ac, N_z, S_csw):
    """
    Transport aircraft wing weight.

    Source
    ------
    Raymer Eq. xx.x

    Applicable to:
        Civil transports

    Inputs
    ------
    W_dg : lb
    S_w  : ft²
    ...

    Returns
    -------
    Wing structural weight (lb)
    """
    Wdg = ac.W_dg
    Sw = ac.S_w
    AR = ac.AR
    tc = ac.tc_root
    lam = ac.taper
    sweep = ac.sweep
    
    W_wing_T =(
        0.0051
        * (Wdg * N_z) ** 0.557
        * Sw ** 0.649
        * AR ** 0.50
        * tc ** (-0.40)
       * (1 + lam) ** 0.10
       * np.cos(sweep) ** (-1.0)
       * S_csw ** 0.10
    )
    return W_wing_T

def horizontal_tail_weight(ac,N_z):
    '''
    F_w = Fuselage width, B_H = HT span , S_ht = H tail area
    NZ = Limit load factor*1.5
    '''
    W_ht = (
        3.316
        * (1 + ac.F_w/ac.B_h)**(-2.0)
        * ((ac.W_dg * N_z)/1000.0)**0.260
        * ac.S_ht**0.806
    )

    return W_ht
def vertical_tail_weight(ac,K_rht,N_z ):
                        
    ''':
    K_rht = = 1.0, L_t = tail arm , A_vt = Vtail aspect ratio
    H_t = H tail height, H_v = Vt tail height, both above fuselage
    M = Mach number,NZ = Limit load factor*1.5
    '''
    W_vt = (
        0.452
        * K_rht
        * (1 + ac.H_t/ac.H_v)**0.5
        * (ac.W_dg * N_z)**0.488
        * ac.S_vt**0.718
        * ac.M**0.341
        * ac.L_t**(-1.0)
        * (1 + ac.S_r/ac.S_vt)**0.348
        * ac.A_vt**0.223
        * (1 + ac.taper)**0.25
        * (np.cos(ac.sweep_vt))**(-0.323)
    )

    return W_vt
def fuselage_weight(ac,K_dwf,N_z):
                    
    '''L = fuselage length
    D = fuselage depth
    W = fuselage width
    '''
    W_fus = (
        0.499
        * K_dwf
        * ac.W_dg**0.35
        * N_z**0.25
        * ac.L**0.5
        * ac.D**0.849
        * ac.W**0.685
    )

    return W_fus
def main_landing_gear_weight(K_cb,
                             K_tpg,
                             W_l,
                             N_l,
                             L_m):

    '''
    K_cb = 2.25 for cross beam, else 1. K_tpg = 0.826 for tri-pod else 1.0
    l_m - length of landing gear
    W_l = Design landing weight, N_l = landing load factor
    '''
    W_mlg = (
        K_cb
        * K_tpg
        * (W_l/N_l)**0.25
        * L_m**0.973
    )

    return W_mlg
def nose_landing_gear_weight(W_l,
                             N_l,
                             L_n,
                             N_nw):

    W_nlg = (
        (W_l/N_l)**0.290
        * L_n**0.5
        * N_nw**0.525
    )

    return W_nlg

def calculate_aircraft_weights(params, mission_reqs, w0_kg, S_m2, S_ht_m2, S_vt_m2, required_tw):
    import math
    # Constants
    LB_TO_KG = 0.45359237
    FT_TO_M = 0.3048
    FT2_TO_M2 = 0.09290304
    KG_TO_LB = 1.0 / LB_TO_KG
    M_TO_FT = 1.0 / FT_TO_M
    M2_TO_FT2 = 1.0 / FT2_TO_M2

    # Aspect ratio & wing parameters
    AR = getattr(params, 'AR', 2.0)
    taper = getattr(params, 'taper', 0.35)
    sweep_rad = math.radians(getattr(params, 'sweep', 29.0))
    tc_root = getattr(params, 'tc_root', 0.04)

    # Geometry estimate
    S_w_ft2 = S_m2 * M2_TO_FT2
    b_m = math.sqrt(S_m2 * AR)
    b_ft = b_m * M_TO_FT

    # HT span B_h
    B_h_m = getattr(params, 'B_h_m', 0.0)
    if B_h_m <= 0.0:
        AR_ht = getattr(params, 'AR_ht', 3.0)
        B_h_m = math.sqrt(S_ht_m2 * AR_ht)
    B_h_ft = B_h_m * M_TO_FT

    # HT height above fuselage
    H_t_ft = getattr(params, 'H_t_m', 0.0) * M_TO_FT

    # VT height above fuselage
    H_v_m = getattr(params, 'H_v_m', 0.0)
    if H_v_m <= 0.0:
        AR_vt = getattr(params, 'AR_vt', 1.5)
        H_v_m = math.sqrt(S_vt_m2 * AR_vt)
    H_v_ft = H_v_m * M_TO_FT

    # Rudder area
    S_r_m2 = getattr(params, 'S_r_m2', 0.0)
    if S_r_m2 <= 0.0:
        S_r_m2 = 0.3 * S_vt_m2
    S_r_ft2 = S_r_m2 * M2_TO_FT2

    # Mach and Arm
    M = getattr(params, 'M_max', 1.5)
    L_t_m = getattr(params, 'L_t_m', 9.0)  # Use L_t_m or default to 9.0
    L_t_ft = L_t_m * M_TO_FT

    # Fuselage
    L_fuse_ft = getattr(params, 'L_fuse_m', 14.0) * M_TO_FT
    D_fuse_ft = getattr(params, 'D_fuse_m', 1.5) * M_TO_FT
    W_fuse_ft = getattr(params, 'W_fuse_m', 1.5) * M_TO_FT
    F_w_ft = getattr(params, 'F_w_m', 1.2) * M_TO_FT

    # Factors
    K_dw = getattr(params, 'K_dw', 1.0)
    K_vs = getattr(params, 'K_vs', 1.0)
    S_csw_m2 = getattr(params, 'S_csw_m2', 0.0)
    if S_csw_m2 <= 0.0:
        S_csw_ft2 = 0.1 * S_w_ft2  # 10% of wing area
    else:
        S_csw_ft2 = S_csw_m2 * M2_TO_FT2
    N_z = getattr(params, 'N_z', 7.5)
    K_rht = getattr(params, 'K_rht', 1.0)
    K_dwf = getattr(params, 'K_dwf', 1.0)

    # Instantiate Aircraft object
    ac = Aircraft(
        W_dg = w0_kg * KG_TO_LB,
        S_w = S_w_ft2,
        AR = AR,
        taper = taper,
        sweep = sweep_rad,
        tc_root = tc_root,
        F_w = F_w_ft,
        B_h = B_h_ft,
        S_ht = S_ht_m2 * M2_TO_FT2,
        H_t = H_t_ft,
        H_v = H_v_ft,
        S_vt = S_vt_m2 * M2_TO_FT2,
        M = M,
        L_t = L_t_ft,
        S_r = S_r_ft2,
        A_vt = getattr(params, 'AR_vt', 1.5),
        sweep_vt = math.radians(getattr(params, 'sweep_vt', 35.0)),
        L = L_fuse_ft,
        D = D_fuse_ft,
        W = W_fuse_ft
    )

    # Each category selects all six equations, never just the wing equation.
    aircraft_type = params.aircraft_type.lower()
    W_l_lb = ac.W_dg * params.landing_weight_fraction
    N_l, L_m_in, L_n_in = params.N_l, params.L_m_in, params.L_n_in
    if aircraft_type == 'fighter':
        W_wing_lb = wing_weight_fighter(ac, K_dw, K_vs, N_z, S_csw_ft2)
        W_ht_lb = horizontal_tail_weight(ac, N_z)
        W_vt_lb = vertical_tail_weight(ac, K_rht, N_z)
        W_fus_lb = fuselage_weight(ac, K_dwf, N_z)
        # Retained user fighter gear correlations; source page not supplied.
        W_mlg_lb = main_landing_gear_weight(params.K_cb, params.K_tpg, W_l_lb, N_l, L_m_in)
        W_nlg_lb = nose_landing_gear_weight(W_l_lb, N_l, L_n_in, params.N_nw)
    elif aircraft_type in ('transport', 'general_aviation', 'ga'):
        W_wing_lb, W_ht_lb, W_vt_lb, W_fus_lb, W_mlg_lb, W_nlg_lb = category_structure(
            ac, params, mission_reqs, N_z, S_csw_ft2, W_l_lb)
    else:
        raise ValueError(f'Unknown aircraft_type: {aircraft_type}')

    # Engine Weight
    if getattr(params, 'propulsion_type', 'jet') == 'jet':
        # Estimate engine weight from required sea-level thrust
        thrust_sl_N = required_tw * w0_kg * getattr(params, 'g', 9.81)
        thrust_sl_lb = thrust_sl_N * 0.22480894
        engine_tw = getattr(params, 'engine_tw', 6.0)
        W_engine_lb = thrust_sl_lb / engine_tw
    else:
        # Turboprop engine weight from required power
        power_sl_kW = required_tw * w0_kg
        engine_pw = getattr(params, 'engine_pw', 3.0) # kW/kg
        W_engine_lb = (power_sl_kW / engine_pw) * KG_TO_LB

    # Systems & Avionics
    systems_fraction = getattr(params, 'systems_fraction', 0.15)
    W_sys_lb = ac.W_dg * systems_fraction

    # Convert to kg
    weights = {
        'W_wing': W_wing_lb * LB_TO_KG,
        'W_ht': W_ht_lb * LB_TO_KG,
        'W_vt': W_vt_lb * LB_TO_KG,
        'W_fus': W_fus_lb * LB_TO_KG,
        'W_mlg': W_mlg_lb * LB_TO_KG,
        'W_nlg': W_nlg_lb * LB_TO_KG,
        'W_engine': W_engine_lb * LB_TO_KG,
        'W_systems': W_sys_lb * LB_TO_KG
    }
    
    return weights



def category_structure(ac, p, mission, nz, scsw, wl):
    """Screenshot Eqs. 15.25–30 and 15.46–51. Internal US units:
    lb, ft, ft², q in lbf/ft², gear length inches, stall speed knots.
    W_press_lb is an explicit input because its defining equation was not supplied.
    """
    from sizing_solver import Atmosphere
    import math
    f = lambda key, default: getattr(p, key, default)
    ht_sweep = math.radians(f('sweep_ht', p.sweep))
    cv, ch, cw = math.cos(ac.sweep_vt), math.cos(ht_sweep), math.cos(ac.sweep)
    tc_h, tc_v = f('tc_ht', p.tc_root), f('tc_vt', p.tc_root)
    sf = f('S_fuse_wet_m2', math.pi * p.L_fuse_m * (p.D_fuse_m+p.W_fuse_m)/2) / .09290304
    nl, lm, ln = p.N_l, p.L_m_in, p.L_n_in
    if p.aircraft_type.lower() == 'transport':
        # Defaults Ky=.3 Lt, Kz=Lt are radii of gyration, not dimensionless factors.
        ky, kz = f('K_y_m', .3*p.L_t_m)/.3048, f('K_z_m', p.L_t_m)/.3048
        se = f('S_e_m2', .25*ac.S_ht*.09290304)/.09290304
        kws = f('K_ws', .75*(1+2*ac.taper)/(1+ac.taper)
                * math.sqrt(ac.S_w*ac.AR)*math.tan(ac.sweep)/ac.L)
        wing = wing_weight_transport(ac, nz, scsw)                  # 15.25
        ht = (.0379*f('K_uht',1)*(1+ac.F_w/ac.B_h)**(-.25)
              * ac.W_dg**.639*nz**.10*ac.S_ht**.75/ac.L_t
              * ky**.704/ch*p.AR_ht**.166*(1+se/ac.S_ht)**.1)       # 15.26
        vt = (.0026*(1+ac.H_t/ac.H_v)**.225*ac.W_dg**.556
              * nz**.536*ac.L_t**(-.5)*ac.S_vt**.5*kz**.875
              /cv*ac.A_vt**.35*tc_v**(-.5))                        # 15.27
        fus = (.3280*f('K_door',1)*f('K_Lg',1)*(ac.W_dg*nz)**.5
               *ac.L**.25*sf**.302*(1+kws)**.04*(ac.L/ac.D)**.10)  # 15.28
        vs = f('V_stall_landing_m_s', math.sqrt(2*p.g*(wl*.45359237)
                 /(1.225*(ac.S_w*.09290304)*f('cl_max_landing',2.5))))/.514444444
        mlg = (.0106*f('K_mp',1)*wl**.888*nl**.25*lm**.4
               *f('N_mw',4)**.321*f('N_mss',2)**(-.5)*vs**.1)      # 15.29
        nlg = .032*f('K_np',1)*wl**.646*nl**.2*ln**.5*p.N_nw**.45  # 15.30
    else:
        rho = Atmosphere.get_properties(f('cruise_altitude_m',0))[2]
        q = f('q_cruise_pa', .5*rho*(mission.velocity/3.6)**2)/47.88025898
        if q <= 0:
            raise ValueError('GA weights require positive cruise dynamic pressure')
        fuel = f('wing_fuel_mass_kg', getattr(mission,'fuel_fraction',0)*ac.W_dg*.45359237)/.45359237
        if fuel < 0:
            raise ValueError('wing_fuel_mass_kg cannot be negative')
        # At zero wing fuel, the tiny fuel correction is unity (dry-wing extension).
        fuel_factor = fuel**.0035 if fuel > 0 else 1.0
        wing = (.036*ac.S_w**.758*fuel_factor*(ac.AR/cw**2)**.6
                *q**.006*ac.taper**.04*(100*ac.tc_root/cw)**(-.3)
                *(nz*ac.W_dg)**.49)                               # 15.46
        ht = (.016*(nz*ac.W_dg)**.414*q**.168*ac.S_ht**.896
              *(100*tc_h/ch)**(-.12)*(p.AR_ht/ch**2)**.043
              *f('taper_ht',ac.taper)**(-.02))                     # 15.47
        vt = (.073*(1+.2*ac.H_t/ac.H_v)*(nz*ac.W_dg)**.376
              *q**.122*ac.S_vt**.873*(100*tc_v/cv)**(-.49)
              *(ac.A_vt/cv**2)**.357*f('taper_vt',ac.taper)**.039)  # 15.48
        fus = (.052*sf**1.086*(nz*ac.W_dg)**.177*ac.L_t**(-.051)
               *(ac.L/ac.D)**(-.072)*q**.241+f('W_press_lb',0))     # 15.49
        mlg = .095*(nl*wl)**.768*(lm/12)**.409                     # 15.50
        nlg = .125*(nl*wl)**.566*(ln/12)**.845                     # 15.51
    return wing, ht, vt, fus, mlg, nlg
