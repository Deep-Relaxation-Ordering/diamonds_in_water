# Breakout Note — Numerical Pilot: Brownian Motion and Sedimentation of Diamond Particles in Aqueous Suspension as a Function of Particle Size

| | |
|---|---|
| **Status** | v0.2 (revised proposal) |
| **Date** | 2026-04-27 |
| **Stewardship** | U. Warring, AG Schätz, Physikalisches Institut Freiburg |
| **Endorsement Marker** | Local stewardship; pilot proposal, not a Coastline. |
| **Council stance of origin** | Guardian (clarity gate on scope and deliverables) |
| **Parent project** | FND suspensions for water-freezing experiments (FP-2 demonstration arm and complementary research arm) |
| **Scope class** | Numerical pilot. No experimental work in this breakout. No phase-transition physics in this breakout. |

*Terminology note: "Guardian council stance", "Coastline", and the "lock/key" framing in §11 are vocabulary from the T(h)reehouse +EC project governance framework. External readers can treat them as project-internal labels for, respectively, a clarity/scope-gate review role, an endorsed cross-project commitment, and the standard-physics primitives a breakout depends on.*

---

## 1. Purpose

A constrained numerical study of two coupled effects — Brownian diffusion and gravitational sedimentation — for diamond particles in water, as a function of particle radius. The purpose is to establish a **quantitative baseline model** against which later experiments and more elaborate simulations can be compared.

This breakout deliberately excludes:
- particle–particle interactions (DLVO, aggregation),
- the liquid–solid phase transition,
- NV-fluorescence spectroscopy,
- any reference to the classical or quantum Mpemba effect,
- temperature gradients within the sample.

These are deferred to follow-up breakouts. Mixing them into this pilot would defeat its purpose, which is to establish the simplest defensible baseline first.

## 2. Scientific question

For diamond particles in water at fixed temperature, in a sample geometry typical of our planned setups (Ibidi µ-slide, ~0.4–1 mm depth; Hellma 110-QS, ~10 mm depth), how does the equilibrium particle distribution and the time evolution toward equilibrium depend on particle radius across the range of practical interest (5 nm – 10 µm)?

In particular:

a. At what particle size does the equilibrium scale height ℓ_g become comparable to the sample depth h, marking the transition from "Brownian-dominated, effectively homogeneous" to "gravitationally stratified"?

b. What is the characteristic equilibration time as a function of radius, and how does it compare to typical experimental observation windows (seconds to hours)?

c. How sensitive are these results to temperature in the range +5 °C to +35 °C, where viscosity η(T) varies by a factor of approximately two?

d. What is the regime map (radius × sample depth × observation time × temperature) that distinguishes (i) homogeneous suspension, (ii) developing gravitational gradient, (iii) fully sedimented?

## 3. Physical model

Single-particle Langevin dynamics in one spatial dimension (vertical), with reflecting boundaries at z = 0 (sample bottom) and z = h (sample top):

```
dz/dt = -v_sed + sqrt(2*D) * ξ(t)
```

with

```
v_sed   = m_eff * g / γ
γ       = 6 * π * η(T) * r
m_eff   = (4/3) * π * r^3 * (ρ_p - ρ_f)
D       = k_B * T / γ
ξ(t)    : zero-mean unit-variance Gaussian white noise
```

Material parameters:
- ρ_p = 3.51 g/cm³ (diamond)
- ρ_f = ρ_water(T) (tabulated)
- η = η_water(T) (Kestin–Sengers correlation or IAPWS formulation)
- T fixed per run, scanned across runs

Reflecting boundaries are appropriate for the closed-bottom, open-top geometry of typical sample chambers (no flux through walls; no escape through air interface on the timescale of interest, since evaporation is slow). Adhesion of particles to the bottom (glass or polymer) or to the air–water interface (top) is a separate loss channel and is documented in §7f as a deferred limitation, not modelled here.

The model assumes:
- non-interacting particles (validated by dilute concentration, c < 10⁻⁴ volume fraction),
- spherical particles (acknowledged simplification — see Section 7),
- no temperature gradient in the sample,
- isothermal water properties throughout the run,
- gravity along z only.

Water density ρ_f(T) varies by only ~0.5 % across 278–308 K and enters the model only through the buoyancy term (ρ_p − ρ_f); the dominant temperature dependence in this study is η(T), which varies by ~2× over the same range. The tabulated ρ_f(T) is nevertheless used (rather than held constant at 4 °C) so that the parameter module exposes a single self-consistent water-properties interface.

These assumptions are explicit so that later breakouts can lift them one at a time.

## 4. Numerical approach

### 4.1 Methods
- **Method A (analytical)**: closed-form evaluation of equilibrium quantities and characteristic timescales — gravitational scale height ℓ_g = k_B T / (m_eff·g), settling velocity v_sed, diffusivity D, equilibration time t_eq ~ min(h, ℓ_g)² / D — for every (r, T, h) cell. These are independent of t_obs. The barometric equilibrium distribution and Einstein–Smoluchowski relations are the textbook expressions evaluated. Time-dependent profiles at finite t_obs are *not* part of Method A; those are produced by Method C.
- **Method B (stochastic ensemble)**: explicit time integration of the overdamped Langevin equation for N = 10⁴–10⁵ independent trajectories. Euler–Maruyama scheme with adaptive timestep dt = min(α · ℓ_g / v_sed, β · ℓ_g² / D), with α, β ~ 10⁻², chosen so that each step resolves both drift across the gravitational scale height and diffusion within it. In the diffusion-dominated regime where ℓ_g ≫ h (small particles), the diffusion bound is replaced by β · h² / D, since h becomes the limiting length scale.
- **Method C (Smoluchowski / overdamped Fokker–Planck PDE)**: numerical solution of the corresponding 1D Smoluchowski equation on the interval [0, h] with no-flux boundaries. Acts both as a cross-check for Method B convergence and as the primary engine for time-dependent quantities across the full t_obs grid in §5.

### 4.2 Implementation
- Language: Python 3.11+
- Core libraries: NumPy, SciPy, Matplotlib
- Vectorization across particles in Method B (no per-particle Python loops)
- Random number generator with explicit seeding for reproducibility
- All parameters in SI units in the source; convenience converters at I/O boundaries only
- Single configuration file in TOML (consistent with `pyproject.toml`) listing all parameter scans for one run

### 4.3 Repository structure
```
fnd-suspension-baseline/
├── README.md
├── LICENSE                        (MIT)
├── pyproject.toml
├── src/
│   ├── parameters.py              (SI constants, water properties)
│   ├── analytical.py              (Method A)
│   ├── langevin.py                (Method B)
│   ├── fokker_planck.py           (Method C)
│   └── regime_map.py              (high-level orchestration)
├── tests/
│   ├── test_water_properties.py
│   ├── test_einstein_relation.py
│   ├── test_equilibrium.py
│   └── test_method_consistency.py
├── notebooks/
│   ├── 01_baseline_validation.ipynb
│   ├── 02_radius_scan.ipynb
│   ├── 03_temperature_scan.ipynb
│   └── 04_regime_map.ipynb
└── docs/
    └── conventions.md
```

### 4.4 Validation strategy
- Method B at long times must reproduce Method A equilibrium to within statistical noise (≤ 2 % deviation in mean height for N = 10⁵).
- Method B and Method C must agree on time-dependent moments (mean, variance) to within numerical tolerance.
- Pure Brownian limit (gravity off): MSD = 2·D·t recovered exactly by Method B.
- Pure sedimentation limit (D → 0): for a uniform initial distribution on [0, h], the bottom population reaches 100 % at t = h / v_sed, with mean arrival time h / (2·v_sed). Particles starting at z₀ arrive at t(z₀) = z₀ / v_sed (Method B).
- Einstein–Smoluchowski recovered: D·γ = k_B·T to machine precision in the parameter module.

## 5. Parameter scan

| Parameter | Range | Step type | Notes |
|---|---|---|---|
| Particle radius r | 5 nm – 10 µm | log-spaced, 30 points | covers nano-FND to micro-diamond regimes |
| Temperature T | 278 K – 308 K | linear, 7 points | +5 °C to +35 °C |
| Sample depth h | 0.1, 0.5, 1, 2, 10 mm | discrete | matches Ibidi (0.4–1 mm), short-path Hellma (2 mm), standard Hellma (10 mm) |
| Observation time t_obs | 1 s, 10 s, 1 min, 1 h, 1 day, 1 week | discrete | spans single measurement to overnight settling |

Total grid: 30 × 7 × 5 × 6 = 6 300 cells. Method A delivers equilibrium quantities for every (r, T, h) cell — these are independent of t_obs and so collapse onto a 30 × 7 × 5 = 1 050 sub-grid. Within that sub-grid the dimensionality is mixed: ℓ_g, v_sed, and D depend only on (r, T) — a 30 × 7 = 210-cell 2-D map — while t_eq ~ min(h, ℓ_g)² / D additionally depends on h, occupying the full 1 050-cell 3-D sub-grid. Method C provides the time-dependent concentration profiles needed across the full 6 300-cell t_obs dimension. Method B is run on a representative subset (~100 cells covering corners and diagonals) as cross-validation against Method C.

### 5.1 Regime classification (used in §6 deliverable 3)

A cell (r, T, h, t_obs) is classified by the concentration profile attained at t_obs:

- **homogeneous**: top-to-bottom ratio c(z=h) / c(z=0) ≥ 0.95 (≤ 5 % gradient — same threshold as deliverable 5);
- **stratified**: 0.05 < c(h) / c(0) < 0.95 (an exponential-like profile with finite scale height);
- **sedimented**: c(h) / c(0) ≤ 0.05 (essentially all particles within ~3 ℓ_g of the bottom).

The 5 % threshold is applied consistently in deliverable 3 (regime map) and deliverable 5 (design table), so that both products tell the same story.

## 6. Deliverables

1. **Working repository** with all code, tests, configuration, and documentation as listed in Section 4.3. Public on GitHub under the T(h)reehouse +EC organization or equivalent.
2. **Validation report** (notebook 01) demonstrating that all four checks in Section 4.4 pass.
3. **Scientific output** (notebooks 02–04):
   - radius dependence of v_sed, D, ℓ_g at room temperature (figure: log–log plot of all three over r);
   - temperature dependence of equilibration time at fixed r (figure: t_eq vs T for representative radii);
   - regime map (figure: phase-diagram-style plot of homogeneous / stratified / sedimented as function of r and t_obs at fixed h, using the §5.1 thresholds).
4. **Short technical note** (this document, refined into v1.0) summarizing findings, suitable as input to later breakouts and as standalone reference for FP-2 students.
5. **Single quantitative table** for design use: "for sample depth h and observation time t, what is the largest r at which the suspension is effectively homogeneous (≤ 5 % concentration gradient between top and bottom)?"

The fifth deliverable is the most important practical output — it directly informs cuvette and protocol choices for any later experimental work.

## 7. Acknowledged limitations of this baseline

These are documented up front so that the model is not later misunderstood as more general than it is:

a. **Spherical-particle assumption.** Real diamond particles, especially in the µm range, are anisotropic and faceted. Stokes drag and Brownian diffusion become tensorial. Correction factors are typically ~1.1 for near-spherical FNDs (HPHT-derived nanocrystals after milling and oxidation, aspect ratio < 1.5) and ~1.3–1.5 for heavily faceted µm microdiamonds. The latter exceeds the 2 % numerical-validation target and so dominates the experimental-comparison error budget at the µm end of the scan — flagged here so that comparisons against the µm regime are read as order-of-magnitude rather than as precision tests.

b. **No particle–particle interactions.** Valid only at dilute concentrations. At fixed mass concentration c_m, the geometric volume fraction is φ = c_m / ρ_p, *independent of particle radius* — so at the Merck stock concentration of 1 mg/mL, φ ≈ 2.85·10⁻⁴ for diamond (ρ_p = 3.51 g/cm³) at any size in the scan range. What changes with r at fixed c_m is the number density (∝ 1/r³), not φ. Working concentrations for tracking experiments will be 10²–10⁴× lower, well within the dilute regime. Note: for r ≲ 10 nm a hydrodynamic hydration shell of ~1 nm contributes a non-negligible correction to the *effective* volume fraction (factor ~2 at r = 5 nm); this remains small relative to the 10²–10⁴× safety margin but is flagged at the small-radius end of the scan. Aggregation is a separate concern and is not modelled here.

c. **No phase transition.** Water remains liquid throughout. The model is invalid in any regime where ice nucleation, growth, or interface motion is relevant. Those phenomena are the subject of a separate, later breakout.

d. **No solute effects.** DI water assumed throughout. Realistic FND suspensions may contain residual surface chemistry or stabilizers; their effect on η or ρ_f is below the modelling threshold here but should be checked against vendor data sheets.

e. **No optical-trapping or external forces.** If experiments later use optical tweezers for active particle positioning, this model does not apply directly to that case.

f. **No particle–wall adhesion.** Reflecting boundaries assume no sticking at the bottom (glass or polymer) or at the air–water interface (top). In real samples — particularly for hydrophobic surface chemistries, low-stability pH, or ionic-strength conditions favouring secondary-minimum capture — both walls can act as loss channels on timescales not necessarily long compared to the equilibration times in this study. Adhesion is *not* modelled here; the regime-map predictions assume full reflection. (Follow-up: "Particle–wall interaction breakout".)

g. **Bulk-diamond density.** ρ_p = 3.51 g/cm³ is the bulk-diamond value. Real FNDs (HPHT-derived nanocrystals with surface graphitic carbon and irradiation-induced lattice defects) may be effectively 3.3–3.4 g/cm³ — a ~3–8 % reduction in (ρ_p − ρ_f) and hence in v_sed and m_eff·g. This sits below the 5 % regime-map threshold and is acceptable for a baseline; vendor-specific densities should be substituted when comparing against specific experimental data.

## 8. Out of scope (explicit deferrals)

The following items are *not* part of this breakout and must not be folded in. Each is listed here so it can be picked up cleanly as a follow-up.

- DLVO-mediated aggregation kinetics. (Follow-up: "Suspension stability and aggregation breakout".)
- Phase-field or sharp-interface modelling of an advancing solidification front. (Follow-up: "Particle–ice front interaction breakout".)
- NV-fluorescence temperature dependence and Brownian thermometry consistency check. (Follow-up: "Dual-mode FND thermometry breakout".)
- Comparison to experimental data. (Follow-up: depends on experimental setup status.)
- Connection to the classical or quantum Mpemba effect. (Follow-up: only after the three above are stable.)

## 9. Timeline and effort estimate

For a single researcher with Python and numerical-methods fluency:

| Phase | Effort | Cumulative |
|---|---|---|
| Repository scaffold, parameters module, water-properties module, tests | 1–2 days | day 2 |
| Method A implementation and notebook 01 (validation) | 1 day | day 3 |
| Method B implementation, vectorized, with tests | 2–3 days | day 6 |
| Method C implementation as cross-check | 2 days | day 8 |
| Parameter scan execution and notebooks 02–04 | 2 days | day 10 |
| Regime map figure and quantitative table (deliverable 5) | 1 day | day 11 |
| Documentation, README, this note refined to v1.0 | 1 day | day 12 |

Total: approximately two working weeks of focused effort. Realistic calendar time with other commitments: four to six weeks.

## 10. Stop conditions

The breakout terminates successfully when:
- all five validation checks (§4.4) pass,
- the parameter scan completes,
- deliverables 1 through 5 exist and are reproducible from a clean clone of the repository.

The breakout terminates unsuccessfully and triggers a re-deliberation if:
- the three numerical methods (A, B, C) cannot be reconciled to within 5 % on any tested cell of the parameter grid,
- the model produces predictions that contradict standard textbook expressions (e.g. Perrin's barometric distribution, Einstein–Smoluchowski),
- effort estimate is exceeded by more than a factor of two without clear scientific justification,
- vendor data sheets for the candidate FND product show η or ρ_f deviations greater than 10 % from pure-water values within the parameter range — in which case the water-properties module is no longer the appropriate baseline and the suspension medium needs to be modelled explicitly.

Either outcome is reportable. A failure mode is not a project failure — it is a finding about the difficulty of the baseline.

## 11. Lock–key declaration

This breakout uses three locks, all standard and broadly endorsed:

- Stokes drag in the low-Reynolds-number limit
- Einstein–Smoluchowski relation between diffusion and friction
- Boltzmann equilibrium distribution under conservative forces

It introduces no new locks. The keys (interpretive frames) are entirely local — they concern only the specific parameter regimes relevant to our experimental geometry.

The output is intended to be useful but unsurprising: it should reproduce known textbook physics in the regime where textbook physics applies, and it should make the radius × depth × time × temperature trade-offs *visible* and *quantitative* for the experimental design that follows.

---

## Change log

- **v0.1** (2026-04-23): initial proposal. Drafted in response to a request for a focused numerical pilot, separated from the broader FND-experiment design discussion. Scope deliberately narrow.
- **v0.2** (2026-04-27): revised in response to internal review (two independent passes). Substantive changes:
  - **Method B timestep (§4.1)**: replaced the undefined stiffness placeholder with two physically grounded constraints — drift across ℓ_g and diffusion within ℓ_g — with the diffusion bound switching to h² / D in the diffusion-dominated small-particle regime.
  - **Method A scope (§4.1, §5)**: clarified that Method A delivers equilibrium and characteristic-timescale quantities only (independent of t_obs); time-dependent grids are produced by Method C, which is now identified as the Smoluchowski (overdamped Fokker–Planck) equation rather than the more general FP form.
  - **Regime classification (§5.1, new)**: added explicit quantitative thresholds (c(h) / c(0) ratios of 0.95 and 0.05) for the homogeneous / stratified / sedimented boundaries, applied consistently across deliverables 3 and 5.
  - **Pure-sedimentation validation (§4.4)**: corrected the arrival-time statement (h / v_sed is the latest, not the universal, arrival; mean is h / (2·v_sed) for uniform initial conditions).
  - **Boundary-condition realism (§3, §7f)**: explicit acknowledgement that reflecting BCs ignore particle–wall adhesion at top (air–water) and bottom (glass), with a new §7f deferring this to a follow-up breakout.
  - **Radius range (§2)**: reconciled with §5; both now read 5 nm – 10 µm.
  - **Volume fraction (§7b)**: rewrote to make explicit that φ = c_m / ρ_p is size-independent at fixed mass concentration; added a small-r note about hydration-shell effects on effective hydrodynamic volume.
  - **Anisotropy correction (§7a)**: separated near-spherical FNDs (~1.1 correction) from heavily faceted µm microdiamonds (~1.3–1.5), and noted that the latter dominates the µm-range error budget.
  - **Density (§7g, new)**: noted that real FNDs may sit at 3.3–3.4 g/cm³ rather than the bulk 3.51 g/cm³ used in the model, with implications below the 5 % threshold.
  - **Water density temperature dependence (§3)**: clarified that ρ_f(T) is used (not held constant) but contributes ≪ η(T) to the temperature-dependent variation.
  - **Closed open decisions**: configuration format → TOML; license → MIT.
  - **Stop conditions (§10)**: added vendor-data-sheet check for η or ρ_f deviations beyond 10 %.
  - **Terminology pointer (top of note)**: brief gloss for "Guardian", "Coastline", "lock/key" framework vocabulary.
  - **Stop-conditions count (§10)**: corrected "all four validation checks" to "all five" to match the §4.4 list.
  - **Method A sub-grid dimensionality (§5)**: noted that ℓ_g, v_sed, and D depend on (r, T) only (a 2-D map) while t_eq additionally depends on h (a 3-D map within the same sub-grid).
