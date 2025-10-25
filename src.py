"""
LAL Digital Demonstrations — ME323
Author: <your name>

This script generates four simple, defensible “digital experiments” for
Laser-Assisted Lubrication (LAL) during high-speed machining of stainless steels.

Outputs (saved in current directory):
  A1_LAL_T_vs_t.png              # E1: T(z=50 µm) vs time for two laser leads
  A2_LAL_Tmax_vs_lead.png        # E1: Peak T vs lead distance sweep
  B1_LAL_force_vs_T.png          # E2: Relative cutting force vs temperature (JC softening)
  C1_MQL_film_vs_flow.png        # E3: Film thickness vs MQL oil flow
  C1_film_thickness_table.csv    # E3 table
  D1_tof_vs_diameter.png         # E4: Droplet time-of-flight vs diameter
  D1_droplet_feasibility.csv     # E4 table

Libraries: numpy, matplotlib, pandas
"""

from math import pi, log
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# =========================
# E1 — 1D transient preheat
# =========================
def e1_preheat_simulation():
    """
    Simple 1D in-depth transient conduction with a moving laser preheat modeled
    as a Gaussian-in-time surface heat flux at a fixed spatial point.

    PDE:  ∂T/∂t = α ∂²T/∂z²
    BC:   -k (∂T/∂z)|_{z=0} = q(t)    (Neumann: imposed surface heat flux)
          ∂T/∂z|_{z=z_max} = 0        (insulated deep boundary)
    IC:   T(z,0) = T0

    Discretization: Explicit FTCS with a ghost node for the Neumann BC.
    Stability (1D heat eq.): dt ≤ dz²/(2α). We use a safety factor.

    Returns: saves A1/A2 figures.
    """
    # --- Material (304/316-ish; constant for simplicity) ---
    k = 15.0                          # W/m-K  (thermal conductivity)
    rho_c = 4.0e6                     # J/m^3-K (ρ c_p)
    alpha = k / rho_c                 # m^2/s   (thermal diffusivity)

    # --- Laser/process parameters (illustrative) ---
    P = 300.0                         # W       (laser power)
    eta = 0.35                        # (–)     (absorptivity)
    r = 0.5e-3                        # m       (Gaussian 1/e^2 radius)
    v = 0.15                          # m/s     (relative speed along lead direction)
    z_shear = 50e-6                   # m       (depth of interest)
    T0 = 298.15                       # K       (ambient)

    # --- Numerical grid ---
    zmax = 1.5e-3                     # m
    dz = 5e-6                         # m
    nz = int(zmax/dz) + 1
    z = np.linspace(0.0, zmax, nz)

    # Explicit stability
    dt = 0.2 * dz**2 / alpha          # s (safety factor 0.2 with the 1/(2α) limit)
    tmax = 0.12                       # s
    nt = int(tmax/dt) + 1
    t = np.linspace(0.0, tmax, nt)

    def pulse_flux_time(t_now, lead):
        """
        Time history of surface heat flux at the observation point as the laser
        spot sweeps past. Gaussian in time with FWHM ~ lead / v.
        Peak flux approximated by 2 η P / (π r^2) for a Gaussian spot.
        """
        t0 = lead / v
        sigma = (lead / v) / 2.355                 # FWHM = 2.355 σ
        qpeak = 2.0 * eta * P / (pi * r**2)        # W/m^2
        return qpeak * np.exp(-0.5 * ((t_now - t0) / sigma)**2)

    def run_case(lead):
        T = np.full(nz, T0)
        Tshear = []
        lam = alpha * dt / dz**2
        for ti in range(nt):
            q = pulse_flux_time(t[ti], lead)
            Tnew = T.copy()
            # interior update
            Tnew[1:-1] = T[1:-1] + lam * (T[2:] - 2*T[1:-1] + T[:-2])
            # Neumann at surface: -k dT/dz = q(t) -> ghost node T(-dz) = T(1) + dz*q/k
            Tghost = T[1] + dz * q / k
            Tnew[0] = T[0] + lam * (T[1] - 2*T[0] + Tghost)
            # insulated deep boundary
            Tnew[-1] = T[-1] + lam * (T[-2] - T[-1])
            T = Tnew
            Tshear.append(np.interp(z_shear, z, T))
        return np.array(Tshear)

    # --- A1: time histories for two leads ---
    leads = [0.3e-3, 0.8e-3]  # 0.3 mm, 0.8 mm
    plt.figure()
    for L in leads:
        TT = run_case(L)
        plt.plot(t, TT - 273.15, label=f'L = {L*1e3:.1f} mm')
    plt.xlabel('Time (s)')
    plt.ylabel('T at z = 50 μm (°C)')
    plt.legend()
    plt.tight_layout()
    plt.savefig('A1_LAL_T_vs_t.png', dpi=200)

    # --- A2: peak temperature vs lead sweep ---
    Ls = np.linspace(0.2e-3, 1.2e-3, 9)
    Tmax = []
    for L in Ls:
        TT = run_case(L)
        Tmax.append(TT.max())
    plt.figure()
    plt.plot(Ls*1e3, np.array(Tmax) - 273.15, marker='o')
    plt.xlabel('Laser lead L (mm)')
    plt.ylabel('Peak T at z = 50 μm (°C)')
    plt.tight_layout()
    plt.savefig('A2_LAL_Tmax_vs_lead.png', dpi=200)


# ===============================================
# E2 — Force reduction via temperature softening
# ===============================================

# ====================================================
# E3 — Is my MQL flow enough? Film thickness estimate
# ====================================================
def e3_film_thickness_plot_and_table():
    """
    Back-of-the-envelope continuity estimate of an equivalent boundary-film
    thickness on the rake face as a function of oil flow:

      h ≈ Vdot_oil / (b * ℓ * U_slip)

    where b = width of cut, ℓ = tool-chip contact length, U_slip = chip sliding speed.
    This is an *upper bound* (assumes all oil reaches and spreads uniformly).
    """
    b = 2.0e-3      # m (width of cut)
    ell = 0.8e-3    # m (contact length)
    U_slip = 60.0   # m/s (chip sliding speed)
    flows_mL_h = np.array([5, 10, 15, 20, 30, 40, 50, 60])  # mL/h (oil only, no air)

    flows_m3_s = flows_mL_h * 1e-6 / 3600.0
    h = flows_m3_s / (b * ell * U_slip)          # m
    h_um = h * 1e6                                # micrometers

    # Save table
    df = pd.DataFrame({
        'Oil flow (mL/h)': flows_mL_h,
        'Estimated film thickness h (µm)': np.round(h_um, 3)
    })
    df.to_csv('C1_film_thickness_table.csv', index=False)

    # Plot
    plt.figure()
    plt.plot(flows_mL_h, h_um, marker='o')
    plt.xlabel('Oil flow rate (mL/h)')
    plt.ylabel('Estimated film thickness h (µm)')
    plt.tight_layout()
    plt.savefig('C1_MQL_film_vs_flow.png', dpi=200)


# =========================================================
# E4 — Will droplets reach the junction? ToF vs evaporation
# =========================================================
def e4_droplet_delivery_plot_and_table():
    """
    Simple 1D Stokes-drag model for droplet flight + D^2-law evaporation time.

    Dynamics under Stokes drag (low Re):
      v(t) = u0 * exp(-t/τ),  where τ = ρ_ℓ d^2 / (18 μ_air)
      s(t) = u0 τ (1 - exp(-t/τ)) -> invert for ToF:
      t_f = -τ ln(1 - s/(u0 τ))   (valid if s < u0 τ; otherwise asymptote)

    Evaporation (very crude upper bound):
      D^2(t) = D0^2 - K t  =>  t_evap = D0^2 / K

    We compare t_f and t_evap for d = 10–40 µm, s = 20 mm, u0 = 30 m/s.
    """
    mu_air = 1.8e-5   # Pa·s
    rho_oil = 900.0   # kg/m^3
    K_evap = 5e-12    # m^2/s  (conservative for oil; slow evaporation)
    s = 0.02          # m      (standoff)
    u0 = 30.0         # m/s    (nozzle exit speed)

    d_um_list = np.array([10, 20, 30, 40])
    rows = []
    tof_ms = []

    for d_um in d_um_list:
        d = d_um * 1e-6
        tau = rho_oil * d**2 / (18.0 * mu_air)              # s
        if s < u0 * tau:
            t_f = -tau * np.log(1.0 - s / (u0 * tau))       # s
        else:
            # If s >= u0 τ, the exponential approach is asymptotic.
            # Use 99% of the distance as a practical proxy.
            t_f = -tau * np.log(1.0 - 0.99)

        t_evap = d**2 / K_evap                               # s
        survives = (t_evap > t_f)
        rows.append({
            'Droplet dia (µm)': d_um,
            'Stokes tau (ms)': round(tau * 1e3, 3),
            'Time-of-flight (ms)': round(t_f * 1e3, 2),
            'Evaporation time (ms)': round(t_evap * 1e3, 0),
            'Survives to target?': 'Yes' if survives else 'No'
        })
        tof_ms.append(t_f * 1e3)

    # Save table
    df = pd.DataFrame(rows)
    df.to_csv('D1_droplet_feasibility.csv', index=False)

    # Plot
    plt.figure()
    plt.plot(d_um_list, tof_ms, marker='o')
    plt.xlabel('Droplet diameter (µm)')
    plt.ylabel('Time-of-flight (ms)')
    plt.tight_layout()
    plt.savefig('D1_tof_vs_diameter.png', dpi=200)


# ============
# Entrypoint
# ============
if __name__ == "__main__":
    e1_preheat_simulation()
    
    e3_film_thickness_plot_and_table()
    e4_droplet_delivery_plot_and_table()
    print("Done. Figures and CSV tables saved.")
