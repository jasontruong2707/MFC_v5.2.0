#!/usr/bin/env python3
"""
3D Lid-Driven Cavity Flow with Herschel-Bulkley Non-Newtonian Fluid

Replicates the benchmark from Chen & Shu (Wiley) for power-law fluids.
Re_gen = rho * U^(2-n) * L^n / K = 400

Setup:
  - Unit cube [0,1]^3, lid at z=1 moves in x-direction at U=1
  - No-slip walls on all 6 faces
  - Power-law fluid: tau0=0, mu = K * gdot^(n-1)
  - n=0.5 (shear-thinning) or n=1.5 (shear-thickening)

Grid: 64^3 (comparable to 61^3 coarsest grid in reference)
"""
import json

eps = 1e-6

# === Flow parameters ===
Re_gen = 400        # Generalized Reynolds number
lid_velocity = 1.0  # Lid velocity (m/s)

# === HB model parameters ===
tau0 = 0.0          # Yield stress (0 for power-law)
nn = 1.5            # Flow behavior index (0.5=shear-thinning, 1.5=shear-thickening)
K = 1.0 / Re_gen    # Consistency index: K = rho * U^(2-n) * L^n / Re = 1/400 = 0.0025
hb_m = 1000.0       # Papanastasiou regularization parameter
mu_bulk = 0.0

# Viscosity bounds depend on n
if nn < 1.0:
    # Shear-thinning: viscosity decreases with shear rate
    # mu = K * gdot^(n-1), n-1 < 0 => high mu at low gdot, low mu at high gdot
    mu_min = K * (1e5)**(nn - 1)   # ~7.9e-6 for n=0.5
    mu_max = K * (0.01)**(nn - 1)  # ~0.025 for n=0.5
else:
    # Shear-thickening: viscosity increases with shear rate
    # mu = K * gdot^(n-1), n-1 > 0 => low mu at low gdot, high mu at high gdot
    mu_min = K * (0.01)**(nn - 1)  # ~2.5e-4 for n=1.5
    mu_max = K * (1e5)**(nn - 1)   # ~0.79 for n=1.5

# Reference Re for MFC (based on consistency index)
Re_ref = 1.0 / K   # = Re_gen = 400

# === Grid ===
N = 63  # 64^3 grid (0-indexed: 0..63 = 64 cells)

# Configuring case dictionary
print(
    json.dumps(
        {
            # Logistics
            "run_time_info": "T",
            # Computational Domain Parameters
            "x_domain%beg": 0.0,
            "x_domain%end": 1.0,
            "y_domain%beg": 0.0,
            "y_domain%end": 1.0,
            "z_domain%beg": 0.0,
            "z_domain%end": 1.0,
            "m": N,
            "n": N,
            "p": N,
            "cfl_adap_dt": "T",
            "cfl_target": 0.5,
            "n_start": 0,
            "t_stop": 40.0,
            "t_save": 20.0,
            # Simulation Algorithm Parameters
            "num_patches": 1,
            "model_eqns": 2,
            "alt_soundspeed": "F",
            "num_fluids": 2,
            "mpp_lim": "F",
            "mixture_err": "T",
            "time_stepper": 3,
            "weno_order": 5,
            "weno_eps": 1e-16,
            "mapped_weno": "T",
            "weno_Re_flux": "T",
            "mp_weno": "T",
            "weno_avg": "T",
            "riemann_solver": 2,
            "wave_speeds": 1,
            "avg_state": 2,
            # Boundary Conditions: no-slip walls on all faces
            "bc_x%beg": -16,
            "bc_x%end": -16,
            "bc_y%beg": -16,
            "bc_y%end": -16,
            "bc_z%beg": -16,
            "bc_z%end": -16,
            # Lid: top wall (z=1) moves in x-direction
            "bc_z%ve1": lid_velocity,
            # Viscous
            "viscous": "T",
            # Formatted Database Files Structure Parameters
            "format": 1,
            "precision": 2,
            "prim_vars_wrt": "T",
            "fd_order": 4,
            "parallel_io": "T",
            # Patch 1: Entire cube (initially at rest)
            "patch_icpp(1)%geometry": 9,
            "patch_icpp(1)%x_centroid": 0.5,
            "patch_icpp(1)%y_centroid": 0.5,
            "patch_icpp(1)%z_centroid": 0.5,
            "patch_icpp(1)%length_x": 1.0,
            "patch_icpp(1)%length_y": 1.0,
            "patch_icpp(1)%length_z": 1.0,
            "patch_icpp(1)%vel(1)": 0.0,
            "patch_icpp(1)%vel(2)": 0.0,
            "patch_icpp(1)%vel(3)": 0.0,
            "patch_icpp(1)%pres": 5e2,
            "patch_icpp(1)%alpha_rho(1)": 0.5,
            "patch_icpp(1)%alpha(1)": 0.5,
            "patch_icpp(1)%alpha_rho(2)": 0.5,
            "patch_icpp(1)%alpha(2)": 0.5,
            # Fluid 1: HB non-Newtonian
            "fluid_pp(1)%gamma": 1.0 / (1.4 - 1.0),
            "fluid_pp(1)%pi_inf": 0.0,
            "fluid_pp(1)%Re(1)": Re_ref,
            "fluid_pp(1)%non_newtonian": "T",
            "fluid_pp(1)%tau0": tau0,
            "fluid_pp(1)%K": K,
            "fluid_pp(1)%nn": nn,
            "fluid_pp(1)%hb_m": hb_m,
            "fluid_pp(1)%mu_min": mu_min,
            "fluid_pp(1)%mu_max": mu_max,
            "fluid_pp(1)%mu_bulk": mu_bulk,
            # Fluid 2: same properties (single-phase effectively)
            "fluid_pp(2)%gamma": 1.0 / (1.4 - 1.0),
            "fluid_pp(2)%pi_inf": 0.0,
            "fluid_pp(2)%Re(1)": Re_ref,
            "fluid_pp(2)%non_newtonian": "T",
            "fluid_pp(2)%tau0": tau0,
            "fluid_pp(2)%K": K,
            "fluid_pp(2)%nn": nn,
            "fluid_pp(2)%hb_m": hb_m,
            "fluid_pp(2)%mu_min": mu_min,
            "fluid_pp(2)%mu_max": mu_max,
            "fluid_pp(2)%mu_bulk": mu_bulk,
        }
    )
)
