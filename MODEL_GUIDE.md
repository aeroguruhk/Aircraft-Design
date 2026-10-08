# Aircraft Sizing V2 — teaching edition

A complete, runnable revision of the supplied fragments. Python 3.10 or later.
The original uploaded files are unchanged. This package is a separate revised version.

## Quick start

Extract the ZIP, and keep all Python modules in the same `aircraft_sizing_v2` folder.
NumPy, SciPy and Matplotlib are the only external dependencies; these are normally
available in Anaconda. If needed: `python -m pip install -r requirements.txt`.

From a terminal in that folder:

```bash
python run_sizing.py --example fighter
python run_sizing.py --example transport
python run_sizing.py --example general_aviation
python run_sizing.py --example transport_turboprop
python run_sizing.py --config examples/transport.json --output-dir my_results
```

In Spyder, open `run_sizing.py`, change `DEFAULT_EXAMPLE` near the top, then Run.
Alternatively, type `main(['--example', 'transport'])` after loading the module.
The default example is fighter; output paths are based on the script folder and
are independent of Spyder's working directory. Plots are saved using the Agg
backend; open the PNGs to inspect them. A repeat run overwrites its own outputs.

Each example is a **single JSON** containing `aircraft`, `mission`, `cases`,
`ws_range`, and an optional `design` section. The previous separate aircraft,
mission and cases JSONs must be migrated to this structure. The two supplied
legacy aircraft JSONs are not drop-in inputs for this revised driver.

## Modules and workflow

1. `Mission_module_SI.py`: mission fuel fraction and statistical initial mass.
2. `sizing_solver.py`: constraint curves and minimum feasible engine ratio.
3. `aircraft_weights.py`: category-specific structural masses and simplified
   engine/system masses; a bracketed root solve closes gross mass.
4. `perform_v3.py`: performance using the same atmosphere, drag polar and engine
   lapse as the constraints; range integrates fuel burn over changing cruise mass.
5. `SFD_BMD_v2.py`: Schrenk aerodynamic wing loads, shear force and bending moment.
6. `run_sizing.py`: configuration, geometry, orchestration and exports.

`Fighter_aicraft_weights_v3.py` remains as a compatibility import for the original
weight-function names. Other module APIs, particularly performance `quick()`,
have deliberately changed; use the revised driver rather than mixing versions.

## What the review changed

- Removed the 1,198-line monolithic driver, hardcoded turboprop filename,
  missing auxiliary-file dependencies and silent default-file creation.
- All six structural equations now dispatch by category. Previously only the
  transport wing dispatched, and its control-surface area was a 20 ft² placeholder.
- Implemented transport screenshot equations **15.25–15.30** and GA screenshot
  equations **15.46–15.51**; the screenshots are included for comparison.
- Corrected mission `Range` from an NM conversion to its documented km unit.
- Replaced propeller use of jet TSFC with BSFC in kg/(kW hour); changing BSFC now
  changes mission fuel burn and performance range.
- Replaced an invalid “units do not matter” statistical-weight assumption with
  an explicitly normalized, configurable teaching fit. Its coefficients are
  example assumptions, not newly transcribed Raymer Table 3.1 coefficients.
- Replaced loosely monitored mass iteration and statistical fallbacks with a
  bracketed physical mass closure and a final closure residual. Failure raises
  a clear error; the software never labels a fallback mass as converged.
- Use trapezoidal mean aerodynamic chord for horizontal-tail volume sizing.
- `N_z` is ultimate throughout. The wing-load module receives it once, without
  a second 1.5 multiplier. Cases exceeding `N_z/1.5` are structurally invalid.
- Invalid/NaN constraints remain binding instead of disappearing through nanmax.
  Clean lift limits are checked on airborne fixed-Mach constraints. A clean
  stall-speed limit is distinct from the landing approach/go-around case.
- Atmosphere uses a troposphere and isothermal layer, valid here to 20 km.
- Unified thrust/shaft-power lapse as `(rho/rho_SL)^0.7`, an explicit simplified
  assumption. Subsonic induced-drag factor is derived from AR and efficiency.
- Removed propeller thrust/power inconsistency, arbitrary fuel clamping, speeds
  below stall in level/climb searches and incorrectly labelled range units.
- A speed capped by `M_max` is labelled `speed_cap_limited`; it is not presented
  as an unconstrained aerodynamic maximum. Service ceiling is reported only
  when an actual 0.5 m/s crossing is bracketed within the sampled 0–15 km range.
- Normalized Schrenk loads to exactly half the total lift on each semispan.

## Inputs and units

JSON geometry is SI unless a name explicitly gives another unit. Component
correlations run internally in their original US units and return **kg masses**.
The `W_*` dictionary names are retained for compatibility; these values are not N.

| Input | Meaning / unit |
|---|---|
| `aircraft_type` | `fighter`, `transport`, `general_aviation` (`ga` alias) |
| `propulsion_type` | `jet` or `turboprop`, independent of aircraft category |
| `mass_kg` | Seed/reference input; final mass is solved |
| `aspect_ratio`, `efficiency` | Wing AR and Oswald efficiency; K = 1/(pi AR e) |
| `N_z`, `N_l` | Ultimate flight and landing load factors |
| `sweep`, `sweep_ht`, `sweep_vt` | Degrees; quarter-chord sweep for correlations |
| `tc_root`, `tc_ht`, `tc_vt` | Thickness/chord ratios; tail defaults use wing value |
| `taper`, `taper_ht`, `taper_vt` | Taper ratios; tail defaults use wing value |
| `L_t_m`, `L_fuse_m`, `D_fuse_m`, `W_fuse_m`, `F_w_m` | Metres; tail arm, fuselage structural length, depth, width, width at tail |
| `V_h`, `V_v` | Horizontal/vertical tail volume coefficients |
| `S_csw_m2` | Wing control-surface area; omitted/zero defaults to 10% wing area |
| `S_e_m2` | Elevator area; defaults to 25% HT area |
| `S_fuse_wet_m2` | Fuselage wetted area; default pi L (depth+width)/2 |
| `K_y_m`, `K_z_m` | Transport radii of gyration, metres; defaults 0.3 Lt and Lt |
| `K_uht`, `K_door`, `K_Lg`, `K_mp`, `K_np` | Explicit transport correction factors; default 1 |
| `K_ws` | Transport wing-sweep factor; calculated from wing geometry when omitted |
| `N_mw`, `N_mss`, `N_nw` | Main wheels, main shock struts, nose wheels; transport defaults 4, 2, 2 |
| `L_m_in`, `L_n_in` | Gear lengths **inches**, kept explicit from original code |
| `V_stall_landing_m_s` | Transport gear correlation stall speed; derived at landing mass if omitted; converted to knots internally |
| `landing_weight_fraction` | Design landing mass / gross mass |
| `H_t_m`, `H_v_m` | Tail heights; VT height estimated from VT area and AR when zero |
| `q_cruise_pa` | GA cruise dynamic pressure override; default 0.5 rho V² |
| `cruise_altitude_m` | Cruise altitude for GA weights and performance |
| `wing_fuel_mass_kg` | GA wing fuel; default all mission fuel in wing, including reserve |
| `W_press_lb` | GA fuselage pressure allowance **lb**, default zero for unpressurized example |
| `maximum_thrust_N`, `maximum_power_kw` | Total installed propulsion capacity, all engines combined |
| `eta_prop` | Propeller efficiency |
| `engine_tw`, `engine_pw` | Engine thrust/weight ratio; engine specific power kW/kg |
| `engine_weight_basis` | `installed` (examples) or `required` (conceptual engine resized to constraints) |
| `systems_fraction` | Aggregate missing systems/accessories mass fraction of gross mass |
| `cl_max_clean`, `cl_max_to`, `cl_max_landing` | Clean/takeoff/landing lift limits |
| `M_max` | Performance sweep Mach cap; not a transonic drag model |
| mission `Range`, `velocity` | Kilometres, true airspeed km/hour |
| mission `Endurance` | Additional loiter hours, not total flight time |
| mission `sfc_cruise`, `sfc_loiter` | Jet TSFC, 1/hour |
| mission `bsfc_cruise`, `bsfc_loiter` | Propeller shaft BSFC, kg/(kW hour) |
| mission `loiter_velocity_kmh` | Propeller loiter true airspeed |
| mission `LBYD`, `cruise_ld_factor` | Prescribed loiter L/D and cruise multiplier |
| mission `Pax`, `Pax_mass`, `payload` | Total people, mass/person kg, additional cargo kg; avoid double counting |
| mission `empty_fraction_reference`, `reference_mass_kg`, `empty_fraction_exponent` | Initial fit: We/W0 = f_ref (mass/mass_ref)^exponent |
| design `wing_loading_pa` | Optional specified point, snapped to nearest grid station |
| design `stall_speed_limit_m_s` | Optional sea-level clean gross-weight stall cap |

GA tail A values use the corresponding tail AR; tail t/c and sweep use the
corresponding surface. Fuselage L/D is interpreted as structural length/depth,
as in the screenshot's geometry notation. The pressure allowance equation was
not supplied, so it remains an explicit input rather than an invented formula.
For a dry GA wing the very small wing-fuel correction is set to unity; this is
an extension outside the positive-fuel formula, stated in the function docstring.

## Example results

See `EXAMPLE_RESULTS.md` for the actual tested masses, geometry and range checks.
These are generic teaching configurations, not calibrated representations of
named real aircraft. Fixed fuselage geometry can produce packaging inconsistencies
if the payload or mission is changed substantially.

The included examples use installed engine mass, so their available-engine
performance belongs to the same physical mass build-up. The required-engine
performance is a comparison at that same aircraft mass; to close a design with
that smaller engine instead, select `engine_weight_basis: "required"` and rerun.

## Outputs

Each example writes to `outputs/<example>/`:

- `summary.json`: mass closure, component masses, engine feasibility, geometry,
  root loads and performance metrics, plus assumptions/warnings.
- `component_weights.csv`, `mass_closure.csv`, `constraints.csv`.
- `wing_loads.csv` and `wing_loads.png`.
- `constraint_diagram.png` and `performance_altitude.png`.
- `performance.csv`, `cruise_range_required.csv`, `cruise_range_available.csv`,
  and `cruise_range.png`.

Mass closure and constraint feasibility do **not** prove mission range. Mission
fuel uses prescribed L/D; performance integrates the actual configured drag polar.
Both the fuel-allocation identity and achievable cruise range are reported, with
warnings if the latter falls short. Cruise range/endurance are alternatives for
the allocated cruise fuel, excluding planned loiter/landing/reserve.

## Limits and remaining source checks

Only the six structural correlations per category were supplied. Engine mass is
an estimate based on thrust/weight or specific power. Systems, nacelles, propellers,
controls and accessories use the aggregate `systems_fraction`; they are not a
complete Raymer component-by-component empty-weight build-up. Review that allowance
before interpreting any absolute mass as an aircraft prediction.

The fighter structural equations, including its landing-gear W_l/N_l factors,
were retained from the uploaded code. Their source page was not supplied. The
fighter main/nose gear formulas need a separate source check; the new transport
and GA gear formulas are transcribed exactly from their supplied screenshots.

Schrenk loads include aerodynamic lift only: no fuel/wing/engine inertia relief,
point loads, torsion, aeroelastic redistribution or spar sizing. Wing-weight and
SFD/BMD calculations are linked by geometry and mass, not a structural stress loop.

Climb uses specific excess power, a small-flight-path-angle approximation; very
large climb rates are illustrative. Jet lapse, propeller efficiency and SFC are
simplified. Supersonic drag switches between supplied polar coefficients at Mach 1;
there is no continuous transonic drag rise or compressibility lift model. An
acceleration case uses mean Mach, not time-integrated acceleration.

Takeoff is a simplified ground-run model, not a full runway/obstacle analysis.
The landing constraint models approach/go-around; the separate landing performance
estimate is ground roll only. Required clean stall limits must be entered explicitly.

## Verification and supplementary references

Run `python -m unittest discover -s tests -v` from the package folder. Twelve tests
cover all twelve new screenshot equations via independent unit-value benchmarks,
GA pressure scaling, km conversion, BSFC response, atmosphere, Schrenk conservation
and analytical bending moment, binding invalid constraints, shared climb equations,
all-example mass closure and an impossible design failure.

Primary equation sources are the three supplied screenshots. Supplementary US-unit
and radius-of-gyration checks used the implementation authors' documentation:

- https://aerosandbox.readthedocs.io/en/master/_modules/aerosandbox/library/weights/raymer_cargo_transport_weights.html
- https://aerosandbox.readthedocs.io/en/master/_modules/aerosandbox/library/weights/raymer_general_aviation_weights.html

No external package from these references is required at runtime.
