import numpy as np
import matplotlib.pyplot as plt

class Atmosphere:
   
    @staticmethod
    def get_properties(h_m):
        """
        Returns:
            temp_K: Temperature in Kelvin
            density_ratio_sigma: Density ratio (rho / rho_sl)
            density_kg_m3: Density in kg/m3
            sound_speed_m_s: Speed of sound in m/s
        """
        if not 0 <= h_m <= 20000:
            raise ValueError('Atmosphere supports altitudes from 0 to 20,000 m')
        R, gravity = 287.05, 9.80665
        if h_m <= 11000:
            temp_K = 288.15 - .0065*h_m
            pressure = 101325*(temp_K/288.15)**(gravity/(R*.0065))
        else:
            temp_K = 216.65
            pressure = 22632.04*np.exp(-gravity*(h_m-11000)/(R*temp_K))
        rho = pressure/(R*temp_K)
        return temp_K, rho/1.225, rho, np.sqrt(1.4*R*temp_K)



class AircraftParameters:
    """
    Holds the aircraft design constants and coefficients.
    """
    def __init__(self, mass_kg=10000.0, aspect_ratio=2.0, efficiency=0.8,
                 cd0_sub=0.016, cd0_sup=0.035, k_sup=0.4, k_sub=0.16,
                 cl_max_to=2.2, s_g_m=900.0, maximum_thrust_N=150000.0,
                 maximum_power_kw=None, g=9.81, rho_sl=1.225, beta_acc=0.8,
                 propulsion_type="jet", eta_prop=0.85, sweep=0.0,
                 taper=0.35, tc_root=0.04, V_h=0.4, V_v=0.07, L_t_m=9.0,
                 AR_ht=3.0, AR_vt=1.5, sweep_vt=35.0, H_t_m=0.0, H_v_m=0.0,
                 S_r_m2=0.0, L_fuse_m=14.0, D_fuse_m=1.5, W_fuse_m=1.5,
                 F_w_m=1.2, N_z=7.5, K_dw=1.0, K_vs=1.0, K_rht=1.0, K_dwf=1.0,
                 landing_weight_fraction=0.78, N_l=3.0, L_m_in=30.0, L_n_in=30.0,
                 K_cb=1.0, K_tpg=1.0, N_nw=2.0, aircraft_type="fighter",
                 systems_fraction=0.15, engine_tw=6.0, engine_pw=3.0, M_max=1.5,
                 **kwargs):
        self.mass = mass_kg
        self.AR = aspect_ratio
        self.e = efficiency
        self.cd0_sub = cd0_sub
        self.cd0_sup = cd0_sup
        self.k_sup = k_sup
        self.k_sub = k_sub
        self.cl_max_to = cl_max_to
        self.s_g = s_g_m
        self.g = g
        self.rho_sl = rho_sl
        self.beta_acc = beta_acc
        self.maximum_thrust_N = maximum_thrust_N
        self.maximum_power_kw = maximum_power_kw
        self.maximum_power_W = maximum_power_kw * 1000.0 if maximum_power_kw is not None else None
        self.propulsion_type = propulsion_type  # "jet" or "turboprop"
        self.eta_prop = eta_prop
        self.sweep = sweep
        
        # Structural weight estimation parameters
        self.taper = taper
        self.tc_root = tc_root
        self.V_h = V_h
        self.V_v = V_v
        self.L_t_m = L_t_m
        self.AR_ht = AR_ht
        self.AR_vt = AR_vt
        self.sweep_vt = sweep_vt
        self.H_t_m = H_t_m
        self.H_v_m = H_v_m
        self.S_r_m2 = S_r_m2
        self.L_fuse_m = L_fuse_m
        self.D_fuse_m = D_fuse_m
        self.W_fuse_m = W_fuse_m
        self.F_w_m = F_w_m
        self.N_z = N_z
        self.K_dw = K_dw
        self.K_vs = K_vs
        self.K_rht = K_rht
        self.K_dwf = K_dwf
        self.landing_weight_fraction = landing_weight_fraction
        self.N_l = N_l
        self.L_m_in = L_m_in
        self.L_n_in = L_n_in
        self.K_cb = K_cb
        self.K_tpg = K_tpg
        self.N_nw = N_nw
        self.aircraft_type = aircraft_type
        self.systems_fraction = systems_fraction
        self.engine_tw = engine_tw
        self.engine_pw = engine_pw
        self.M_max = M_max

        # Set any additional kwargs passed in
        for key, value in kwargs.items():
            setattr(self, key, value)

        # Subsonic K calculation
        self.k_sub = 1.0 / (np.pi * self.AR * self.e)

    @classmethod
    def from_json(cls, filepath):
        import json
        with open(filepath, 'r') as f:
            config = json.load(f)
        data = config.get('data', config)  # Support both old and new format
        return cls(**data)


class SizingCase:
    """
    Base class representing a constraint case.
    """
    def __init__(self, name, altitude_m, mach=None, load_factor=1.0, beta=0.78):
        self.name = name
        self.altitude = altitude_m
        self.mach = mach
        self.n = load_factor
        self.beta = beta

    def calculate_tw(self, ws_range, params):
        """
        Calculate required T/W at sea level for a range of wing loadings.
        """
        raise NotImplementedError


class LevelFlightCase(SizingCase):
    """
    Case 1 & 2: Constant altitude, constant speed level flight (or dash).
    """
    def calculate_tw(self, ws_range, params):
        temp, sigma, rho, a = Atmosphere.get_properties(self.altitude)
        V = self.mach * a
        q = 0.5 * rho * V ** 2
        alpha = sigma ** 0.7
        
        # Determine aerodynamic coefficients
        is_supersonic = self.mach >= 1.0
        
       
        cd0 = params.cd0_sup if is_supersonic else params.cd0_sub
        K = params.k_sup if is_supersonic else params.k_sub
            
        const1 = q * cd0 / self.beta
        const2 = K * self.beta * (self.n ** 2) / q
        
        # T/W calculation
        term1 = const1 / ws_range
        term2 = const2 * ws_range
        
        tw_sl = (term1 + term2) * self.beta / alpha
        if params.propulsion_type == "jet":
            return tw_sl
        elif params.propulsion_type == "turboprop":
            return tw_sl * params.g * V / (params.eta_prop * 1000.0)
        

class SustainedTurnCase(SizingCase):
    """
    Case 3: Sustained g-turn at constant altitude and speed.
    """
    def calculate_tw(self, ws_range, params):
        temp, sigma, rho, a = Atmosphere.get_properties(self.altitude)
        V = self.mach * a
        q = 0.5 * rho * V ** 2
        alpha = sigma ** 0.7
        
        is_supersonic = self.mach >= 1.0
        cd0 = params.cd0_sup if is_supersonic else params.cd0_sub
        K = params.k_sup if is_supersonic else params.k_sub
        
        const1 = q * cd0 / self.beta
        const2 = K * self.beta * (self.n ** 2) / q
        
        term1 = const1 / ws_range
        term2 = const2 * ws_range
        
        tw_sl = (term1 + term2) * self.beta / alpha
        if params.propulsion_type == "jet":
            return tw_sl
        elif params.propulsion_type == "turboprop":
            return tw_sl * params.g * V / (params.eta_prop * 1000.0)


class AccelerationCase(SizingCase):
    """
    Case 4: Horizontal acceleration from Mach M1 to M2 in time dt.
    """
    def __init__(self, name, altitude_m, m1=0.8, m2=1.2, dt_s=40.0, load_factor=1.0, beta=0.78):
        # Average Mach is evaluated
        m_avg = 0.5 * (m1 + m2)
        super().__init__(name, altitude_m, mach=m_avg, load_factor=load_factor, beta=beta)
        self.m1 = m1
        self.m2 = m2
        self.dt = dt_s

    def calculate_tw(self, ws_range, params):
        temp, sigma, rho, a = Atmosphere.get_properties(self.altitude)
        V_avg = self.mach * a
        q = 0.5 * rho * V_avg ** 2
        alpha = sigma ** 0.7
        
        # SPREADSHEET BUG: Case 4 (Mach average 1.0) uses subsonic coefficients
       
        is_supersonic = self.mach >= 1.0
        cd0 = params.cd0_sup if is_supersonic else params.cd0_sub
        K = params.k_sup if is_supersonic else params.k_sub
            # Correct acceleration term (no division by alpha or multiplication by beta_acc inside const3)
        const3 = (a / params.g) * (self.m2 - self.m1) / self.dt
            
        const1 = q * cd0 / self.beta
        const2 = K * self.beta * (self.n ** 2) / q
        
        term1 = const1 / ws_range
        term2 = const2 * ws_range
        
        
            # Correct scaling: only scale the total required local T/W by beta / alpha
            # local T/W = drag/weight + acc/g = (term1*beta + term2*beta + const3)
            # wait, local drag/weight = q*Cd0/(W_local/S) + K*(W_local/S)*n^2/q
            # Since W_local = beta * W_TO:
            # local drag/weight = q*Cd0/(beta * W/S) + K*beta*(W/S)*n^2/q
            # required local T/W = q*Cd0/(beta * W/S) + K*beta*(W/S)*n^2/q + (1/g)*dV/dt
            # required SL T/W = local T/W * beta / alpha
        tw_sl = (term1 + term2 + const3) * self.beta / alpha
        if params.propulsion_type == "jet":
            return tw_sl
        elif params.propulsion_type == "turboprop":
            return tw_sl * params.g * V_avg / (params.eta_prop * 1000.0)
        else:
            raise ValueError(f"Unsupported propulsion_type: {params.propulsion_type}")
            
        


class TakeoffCase(SizingCase):
    """
    Takeoff ground run S_G.
    Allows specifying custom takeoff roll distance (s_g_m) and max lift coefficient (cl_max_to).
    If omitted, falls back to global aircraft parameters.
    """
    def __init__(self, name, altitude_m=0.0, s_g_m=None, cl_max_to=None, load_factor=1.0, beta=0.78):
        super().__init__(name, altitude_m, mach=None, load_factor=load_factor, beta=beta)
        self.s_g = s_g_m
        self.cl_max_to = cl_max_to

    def calculate_tw(self, ws_range, params):
        temp, sigma, rho, a = Atmosphere.get_properties(self.altitude)
        s_g = self.s_g if self.s_g is not None else params.s_g
        cl_max_to = self.cl_max_to if self.cl_max_to is not None else params.cl_max_to
        
        tw_sl = 1.4 * ws_range / (params.g * cl_max_to * s_g)
        if params.propulsion_type == "jet":
            return tw_sl
        elif params.propulsion_type == "turboprop":
            V = 1.2 * np.sqrt(2.0 * ws_range / (rho * cl_max_to))
            return tw_sl * params.g * V / (params.eta_prop * 1000.0)
        else:
            raise ValueError(f"Unsupported propulsion_type: {params.propulsion_type}")


class LandingCase(SizingCase):
    """
    Landing and approach case.
    Uses approach speed and lift limit to estimate low-speed thrust requirement.
    """
    def __init__(self, name, altitude_m=0.0, mach=None, cl_max_landing=2.8,
                 approach_factor=1.3, climb_grad=0.03, load_factor=1.0, beta=0.78):
        super().__init__(name, altitude_m, mach=mach, load_factor=load_factor, beta=beta)
        self.cl_max_landing = cl_max_landing
        self.approach_factor = approach_factor
        self.climb_grad = climb_grad

    def calculate_tw(self, ws_range, params):
        temp, sigma, rho, a = Atmosphere.get_properties(self.altitude)
        alpha = sigma ** 0.7

        if self.mach is None:
            # Use landing approach speed based on stall speed
            v_stall = np.sqrt(2.0 * self.beta * ws_range / (rho * self.cl_max_landing))
            V = self.approach_factor * v_stall
        else:
            V = self.mach * a

        q = 0.5 * rho * V ** 2
        cd0 = params.cd0_sub
        K = params.k_sub

        Cl = self.beta * ws_range / q
        if self.mach is not None:
            Cl = np.where(Cl <= self.cl_max_landing, Cl, np.nan)

        cd = cd0 + K * Cl ** 2
        tw_sl = (q * cd / ws_range + self.beta*self.climb_grad) / alpha
        if params.propulsion_type == "jet":
            return tw_sl
        elif params.propulsion_type == "turboprop":
            return tw_sl * params.g * V / (params.eta_prop * 1000.0)
       

class ClimbCase(SizingCase):
    """
    Rate of Climb (ROC) constraint.
    Supports climb at maximum L/D speed (when mach=None) or at a specified Mach number.
    Works for both Jet (T/W) and Turboprop (P/W).
    """
    def __init__(self, name, altitude_m, roc_m_s=200.0, mach=None, load_factor=1.0, beta=1.0):
        super().__init__(name, altitude_m, mach=mach, load_factor=load_factor, beta=beta)
        self.roc = roc_m_s

    def calculate_tw(self, ws_range, params):
        temp, sigma, rho, a = Atmosphere.get_properties(self.altitude)
        alpha = sigma ** 0.7

        is_supersonic = (self.mach >= 1.0) if (self.mach is not None) else False
        cd0 = params.cd0_sup if is_supersonic else params.cd0_sub
        K = params.k_sup if is_supersonic else params.k_sub

        if self.mach is None:
            # Climb at max L/D speed
            q = self.beta * self.n * ws_range / min(np.sqrt(cd0 / K), getattr(params, "cl_max_clean", 1.5))
            V = np.sqrt(2*q/rho)
            drag_ratio = q*cd0/ws_range + K*(self.beta*self.n)**2*ws_range/q
            tw_sl = (drag_ratio + self.beta*self.roc/V)/alpha
        else:
            # Climb at a constant specified Mach number
            V = self.mach * a
            q = 0.5 * rho * V ** 2

            term1 = q * cd0 / (self.beta * ws_range)
            term2 = K * self.beta * (self.n ** 2) * ws_range / q
            climb_rate_term = self.roc / V

            tw_sl = (term1 + term2 + climb_rate_term) * self.beta / alpha

        # 🔑 Switch between Jet and Turboprop
        if params.propulsion_type == "jet":
            return tw_sl
        elif params.propulsion_type == "turboprop":
            return tw_sl * params.g * V / (params.eta_prop * 1000.0)

class ServiceCeilingCase(ClimbCase):
    """
    Service ceiling sizing case.
    Uses a small required rate of climb at high altitude to build a ceiling constraint.
    """
    def __init__(self, name, altitude_m, roc_m_s=0.5, mach=None, load_factor=1.0, beta=0.78):
        super().__init__(name, altitude_m, roc_m_s=roc_m_s, mach=mach, load_factor=load_factor, beta=beta)


class TurnRateCase(SizingCase):
    """
    Sizing constraint for a specified Turn Rate (deg/s).
    Supports Sustained Turn Rate (STR, when ps_m_s = 0) and
    Attained Turn Rate (ATR, trading energy for turn rate when ps_m_s < 0).
    """
    def __init__(self, name, altitude_m, mach, turn_rate_deg_s, ps_m_s=0.0, cl_max=None, beta=0.78):
        super().__init__(name, altitude_m, mach=mach, load_factor=1.0, beta=beta)
        self.turn_rate = turn_rate_deg_s
        self.ps = ps_m_s
        self.cl_max = cl_max

    def calculate_tw(self, ws_range, params):
        temp, sigma, rho, a = Atmosphere.get_properties(self.altitude)
        V = self.mach * a
        q = 0.5 * rho * V ** 2
        alpha = sigma ** 0.7
        
        # Convert turn rate to rad/s
        omega_rad = np.radians(self.turn_rate)
        
        # Calculate required load factor n
        n = np.sqrt(1.0 + (omega_rad * V / params.g)**2)
        self.n = n # Save calculated load factor for reference
        
        is_supersonic = self.mach >= 1.0
        cd0 = params.cd0_sup if is_supersonic else params.cd0_sub
        K = params.k_sup if is_supersonic else params.k_sub
        
        const1 = q * cd0 / self.beta
        const2 = K * self.beta * (n ** 2) / q
        
        term1 = const1 / ws_range
        term2 = const2 * ws_range
        ps_term = self.ps / V
        
        tw_sl = (term1 + term2 + ps_term) * self.beta / alpha
        
        # Aerodynamic lift limit check
        cl = n * self.beta * ws_range / q
        if self.cl_max is not None:
            # Mask out values exceeding CL_max by setting them to NaN
            tw_sl = np.where(cl <= self.cl_max, tw_sl, np.nan)
        if params.propulsion_type == "jet":
            return tw_sl
        elif params.propulsion_type == "turboprop":
            return tw_sl * params.g * V / (params.eta_prop * 1000.0)
        


class SizingSolver:
    """
    Main solver class that aggregates the design parameters, runs the case constraints,
    solves for the design point, and generates plots.
    """
    def __init__(self, params, cases=None):
        self.params = params
        self.cases = cases or []

    def solve(self, ws_range):
        """
        Solves all constraint cases for the given range of wing loadings.
        Returns:
            dict of {case_name: tw_array}
        """
        results = {}
        for case in self.cases:
            tw_array = case.calculate_tw(ws_range, self.params)
            # Every airborne constraint must meet its lift and structural limits.
            if not isinstance(case, (TakeoffCase, LandingCase)):
                if case.n > self.params.N_z / 1.5:
                    tw_array = np.full_like(ws_range, np.nan, dtype=float)
                elif case.mach is not None:
                    _, _, rho, a = Atmosphere.get_properties(case.altitude)
                    q = .5*rho*(case.mach*a)**2
                    cl_max = getattr(case, 'cl_max', None) or getattr(self.params, 'cl_max_clean', 1.5)
                    tw_array = np.where(case.n*case.beta*ws_range/q <= cl_max, tw_array, np.nan)
            results[case.name] = tw_array
        return results

    def get_design_point(self, ws_range, tw_results):
        """
        Finds the critical (maximum required T/W) curve at each W/S.
        Then, identifies the design point (minimum T/W that satisfies all constraints).
        Since higher wing loading is generally preferred for smaller wing area (lower weight/drag),
        the design point is usually determined by the intersection of the critical constraint curve
        with the wing loading limit.
        If we select a target W/S (e.g. 3250 Pa as in the spreadsheet), we find the required T/W.
        """
        # Build array of maximum required T/W across all cases at each W/S
        all_tw = np.array([tw_results[name] for name in tw_results])
        critical_tw = np.max(np.where(np.isfinite(all_tw), all_tw, np.inf), axis=0)
        return critical_tw

    def plot_diagram(self, ws_range, tw_results, critical_tw, design_ws=3250.0, design_tw=1.2, available_tw=None, save_path=None):
        """
        Generates a premium matplotlib chart with case identifiers and shaded design space.
        """
        plt.style.use('seaborn-v0_8-whitegrid' if 'seaborn-v0_8-whitegrid' in plt.style.available else 'default')
        fig, ax = plt.subplots(figsize=(11, 7.5), dpi=150)

        # Define premium color palette
        colors = {
            "Case 1": "#1f77b4",          # Blue
            "Case 2": "#aec7e8",          # Light blue
            "Case 3": "#ff7f0e",          # Orange
            "Case 4": "#ffbb78",          # Light orange
            "Case 5": "#2ca02c",          # Green
            "Case 6": "#d62728",          # Red
            "Case 7": "#9467bd",          # Purple
            "Case 8": "#bcbd22",          # Olive-green
            "Case 9": "#17becf",          # Teal
            "Critical": "#222222"         # Dark Grey/Black for critical boundary
        }

        # Plot each case, mapping colors using substring matching to accommodate descriptive names
        for name, tw in tw_results.items():
            color = None
            for key, val in colors.items():
                if key.lower() in name.lower():
                    color = val
                    break
            ax.plot(ws_range, tw, label=name, linewidth=2, color=color)

        tw_label = "T/W" if self.params.propulsion_type == "jet" else "P/W"

        # Plot critical boundary
        ax.plot(ws_range, critical_tw, label=f"Critical Boundary ({tw_label} Min)", color=colors["Critical"], linestyle="--", linewidth=2.5)

        # Use a scale that matches the actual data range instead of forcing a large fixed y-axis floor.
        max_crit = float(np.nanmax(critical_tw)) if np.size(critical_tw) else 1.0
        max_case = float(np.nanmax([np.nanmax(tw) for tw in tw_results.values()])) if tw_results else max_crit
        ylim_max = max(max_crit, max_case)

        if available_tw is not None:
            ylim_max = max(ylim_max, float(available_tw))

        if not np.isfinite(ylim_max) or ylim_max <= 0:
            ylim_max = 1.0

        # Add 15% headroom above the largest value, but avoid artificial floors like 3.0.
        ylim_max = float(ylim_max * 1.15)

        ax.set_ylim(0.0, ylim_max)

        if available_tw is not None:
            # Plot upper bound and shade between critical and available T/W or P/W
            ax.axhline(y=available_tw, color='darkred', linestyle='-.', linewidth=2, label=f"Max Available {tw_label} ({available_tw:.2f})")
            ax.fill_between(ws_range, critical_tw, available_tw, where=(critical_tw <= available_tw), color='#2ca02c', alpha=0.12, label='Feasible Design Space')
        else:
            # Shade the feasible design space (region above the critical constraint boundary)
            ax.fill_between(ws_range, critical_tw, ylim_max, color='#2ca02c', alpha=0.12, label='Feasible Design Space')

        # Plot design point
        ax.plot(design_ws, design_tw, marker='o', color='magenta', markersize=10, label=f"Design Point ({design_ws:.0f} Pa, {design_tw:.3f})")
        ax.annotate(f"Design Point\nW/S = {design_ws:.0f} Pa\n{tw_label} = {design_tw:.3f}",
                    xy=(design_ws, design_tw), xytext=(design_ws - 600, design_tw + (ylim_max * 0.05)),
                    arrowprops=dict(facecolor='black', shrink=0.05, width=1, headwidth=6),
                    fontsize=10, fontweight='bold', bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="magenta", lw=1.5))

        # Style chart
        ax.set_title(f"Aircraft Sizing Constraint Diagram ({tw_label} vs W/S)", fontsize=14, fontweight='bold', pad=15)
        ax.set_xlabel("Wing Loading, W/S (N/m²)", fontsize=11, labelpad=8)
        if self.params.propulsion_type == "jet":
            ax.set_ylabel("Thrust-to-Weight Ratio, T/W", fontsize=11, labelpad=8)
        else:
            ax.set_ylabel("Power-to-Weight Ratio, P/W (kW/kg)", fontsize=11, labelpad=8)
        ax.set_xlim(np.min(ws_range), np.max(ws_range))

        # Add legend and grid
        ax.legend(loc="upper right", frameon=True, facecolor="white", edgecolor="none", shadow=True)
        ax.grid(True, linestyle=":", alpha=0.6)

        # Add a water mark / design info
        ax.text(0.02, 0.05, "Conceptual Aircraft Sizing Tool\nMattingly Equations",
                transform=ax.transAxes, fontsize=9, alpha=0.4, verticalalignment='bottom')

        plt.tight_layout()
        if save_path:
            plt.savefig(save_path)
            print(f"Plot saved successfully to: {save_path}")

        return fig, ax
