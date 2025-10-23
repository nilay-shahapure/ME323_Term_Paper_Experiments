

---

## Digital Demonstrations — What each “experiment” does

These are **simple, defensible** models that generate the figures/tables used in the report. They’re built for intuition (not high-fidelity prediction), using stainless-like properties.

> Run: `src.py`
> Outputs: PNG figures and CSV tables listed under each experiment.

---

### E1 — Lead distance vs preheat at the shear plane

**Goal:** Show why **laser lead distance** (L) matters.

**Model (1D transient conduction):** The workpiece is a semi-infinite solid in depth (z). Temperature (T(z,t)) obeys
[
\frac{\partial T}{\partial t}=\alpha\frac{\partial^2 T}{\partial z^2},\qquad
-k,\Bigl.\frac{\partial T}{\partial z}\Bigr|_{z=0}=q(t),\quad T(z,0)=T_0.
]
The surface heat flux (q(t)) is a **Gaussian pulse** representing the laser spot passing the point; its **FWHM** equals (L/v) (lead/relative speed). Peak flux is approximated by
[
q_0 \approx \frac{2,\eta P}{\pi r^2},
]
with laser power (P), absorptivity (\eta), and Gaussian radius (r).

**Why 1D is OK here:** On sub-ms times and shallow depths ((\sim 50,\mu)m), depth-wise conduction dominates; 1D captures the **trend**.

**What to look for in the plots:**

* `A1_LAL_T_vs_t.png` — (T) at (z=50,\mu)m vs time for two leads (e.g., 0.3 and 0.8 mm). Longer lead → broader/higher thermal pulse.
* `A2_LAL_Tmax_vs_lead.png` — Peak (T) vs (L); you’ll see a **sweet region** where the peak aligns with the shear plane.



**Limitations:** No convection, no lateral heat spreading, constant properties (so absolute values are upper bounds).

---

### E2 — Temperature softening → force reduction (illustrative)

**Goal:** Visualize how raising shear-zone temperature can **reduce cutting force**.

**Model (Johnson–Cook-style softening):**
[
\sigma(T)=(A+B\epsilon^n),\bigl(1+C\ln\dot{\epsilon}\bigr),\bigl[1-(T^*)^m\bigr],\qquad
T^*=\frac{T-T_r}{T_m-T_r}.
]
Holding geometry constant, assume (F_c \propto \sigma(T)) and plot the **relative** force vs (T).

**What to look for:**

* `B1_LAL_force_vs_T.png` — A 200–300 °C rise typically shows a **~10–25 %** drop in the **relative** force in this toy model; 400–600 °C gives **~27–45 %**. This matches the **direction/order** seen in LAL experiments.



**Limitations:** Generic parameters, constant chip geometry → **illustrative**, not predictive.

---

### E3 — Is my MQL flow enough? (film adequacy back-of-the-envelope)

**Goal:** Check whether the chosen **oil flow rate** can plausibly sustain a boundary film on the rake face.

**Model (volume continuity):**
[
h ;\approx; \frac{\dot V_{\text{oil}}}{,b,\ell,U_{\text{slip}},},
]
where (h) is an **equivalent** average film thickness, (b) = width of cut, (\ell) = tool–chip contact length, (U_{\text{slip}}) = chip sliding speed along the rake face.

**What to look for:**

* `C1_MQL_film_vs_flow.png` — Film thickness vs oil flow (5–60 mL/h).
* `C1_film_thickness_table.csv` — Same as a table.



**Limitations:** Uniform spreading is idealized; ignores spray losses and chip shielding.

---

### E4 — Will droplets reach the junction before evaporating?

**Goal:** Sanity-check that **MQL droplets** actually arrive at the tool–chip interface.

**Model (1D Stokes drag + (D^2) evaporation):**

* Kinematics under Stokes drag:
  [
  \tau=\frac{\rho_\ell d^2}{18,\mu_{\text{air}}},\quad
  v(t)=u_0,e^{-t/\tau},\quad
  s(t)=u_0\tau\bigl(1-e^{-t/\tau}\bigr),
  ]
  so the **time-of-flight** is
  [
  t_f=-\tau\ln!\left(1-\frac{s}{u_0\tau}\right)\quad (\text{if } s<u_0\tau).
  ]
* Crude evaporation time (upper bound): (t_{\text{evap}}=D_0^2/K).

**What to look for:**

* `D1_tof_vs_diameter.png` — ToF (ms) vs droplet diameter (10–40 µm) for (s=20) mm, (u_0=30) m/s.
* `D1_droplet_feasibility.csv` — ToF vs (t_{\text{evap}}) (“Survives?” yes/no).



**Limitations:** 1D, no cross-flow or breakup; used for feasibility, not detailed trajectories.


