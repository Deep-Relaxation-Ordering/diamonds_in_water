# Breakout Note — Numerical Pilot: Brownian Motion and Sedimentation of Diamond Particles in Aqueous Suspension as a Function of Particle Size

| | |
|---|---|
| **Status** | v0.1 (proposal) |
| **Date** | 2026-04-23 |
| **Stewardship** | U. Warring, AG Schätz, Physikalisches Institut Freiburg |
| **Endorsement Marker** | Local stewardship; pilot proposal, not a Coastline. |
| **Council stance of origin** | Guardian (clarity gate on scope and deliverables) |
| **Parent project** | FND suspensions for water-freezing experiments (FP-2 demonstration arm and complementary research arm) |
| **Scope class** | Numerical pilot. No experimental work in this breakout. No phase-transition physics in this breakout. |

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

For diamond particles in water at fixed temperature, in a sample geometry typical of our planned setups (Ibidi µ-slide, ~0.4–1 mm depth; Hellma 110-QS, ~10 mm depth), how does the equilibrium particle distribution and the time evolution toward equilibrium depend on particle radius across the range of practical interest (10 nm – 10 µm)?

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

Reflecting boundaries are appropriate for the closed-bottom, open-top geometry of typical sample chambers (no flux through walls; no escape through air interface on the timescale of interest, since evaporation is slow).

The model assumes:
- non-interacting particles (validated by dilute concentration, c < 10⁻⁴ volume fraction),
- spherical particles (acknowledged simplification — see Section 7),
- no temperature gradient in the sample,
- isothermal water properties throughout the run,
- gravity along z only.

These assumptions are explicit so that later breakouts can lift them one at a time.

## 4. Numerical approach

### 4.1 Methods
- **Method A (analytical)**: closed-form equilibrium distribution (barometric exponential), mean first-passage times, and ensemble drift–diffusion comparison. No simulation; pure evaluation of textbook expressions across the parameter grid.
- **Method B (stochastic ensemble)**: explicit time integration of the Langevin equation for N = 10⁴–10⁵ independent trajectories. Euler–Maruyama scheme with adaptive timestep dt = min(0.01·γ/k, h²/(100·D)) where k is a stiffness placeholder.
- **Method C (Fokker–Planck / drift–diffusion PDE)**: numerical solution of the corresponding 1D Fokker–Planck equation on the interval [0, h] with no-flux boundaries. Used as cross-check for Method B convergence.

### 4.2 Implementation
- Language: Python 3.11+
- Core libraries: NumPy, SciPy, Matplotlib
- Vectorization across particles in Method B (no per-particle Python loops)
- Random number generator with explicit seeding for reproducibility
- All parameters in SI units in the source; convenience converters at I/O boundaries only
- Single configuration file (TOML or YAML) listing all parameter scans for one run

### 4.3 Repository structure
```
fnd-suspension-baseline/
├── README.md
├── LICENSE                        (MIT or BSD-3-Clause; to be decided)
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
- Pure sedimentation limit (D → 0): all particles arrive at z = 0 at t = h/v_sed (Method B).
- Einstein–Smoluchowski recovered: D·γ = k_B·T to machine precision in the parameter module.

## 5. Parameter scan

| Parameter | Range | Step type | Notes |
|---|---|---|---|
| Particle radius r | 5 nm – 10 µm | log-spaced, 30 points | covers nano-FND to micro-diamond regimes |
| Temperature T | 278 K – 308 K | linear, 7 points | +5 °C to +35 °C |
| Sample depth h | 0.1, 0.5, 1, 2, 10 mm | discrete | matches Ibidi (0.4–1 mm), short-path Hellma (2 mm), standard Hellma (10 mm) |
| Observation time t_obs | 1 s, 10 s, 1 min, 1 h, 1 day, 1 week | discrete | spans single measurement to overnight settling |

Total grid: 30 × 7 × 5 × 6 = 6 300 cells. Each cell evaluated by Method A (instantaneous), and a representative subset (~100 cells covering corners and diagonals) by Method B and C for cross-validation.

## 6. Deliverables

1. **Working repository** with all code, tests, configuration, and documentation as listed in Section 4.3. Public on GitHub under the T(h)reehouse +EC organization or equivalent.
2. **Validation report** (notebook 01) demonstrating that all four checks in Section 4.4 pass.
3. **Scientific output** (notebooks 02–04):
   - radius dependence of v_sed, D, ℓ_g at room temperature (figure: log–log plot of all three over r);
   - temperature dependence of equilibration time at fixed r (figure: t_eq vs T for representative radii);
   - regime map (figure: phase-diagram-style plot of homogeneous / stratified / sedimented as function of r and t_obs at fixed h).
4. **Short technical note** (this document, refined into v1.0) summarizing findings, suitable as input to later breakouts and as standalone reference for FP-2 students.
5. **Single quantitative table** for design use: "for sample depth h and observation time t, what is the largest r at which the suspension is effectively homogeneous (≤ 5 % concentration gradient between top and bottom)?"

The fifth deliverable is the most important practical output — it directly informs cuvette and protocol choices for any later experimental work.

## 7. Acknowledged limitations of this baseline

These are documented up front so that the model is not later misunderstood as more general than it is:

a. **Spherical-particle assumption.** Real diamond particles, especially in the µm range, are anisotropic and faceted. Stokes drag and Brownian diffusion become tensorial. Correction factors are typically 1.1–1.5 for mild anisotropy; this is within the systematic uncertainty band of the model and acceptable for a baseline.

b. **No particle–particle interactions.** Valid only at dilute concentrations. At Merck stock concentration (1 mg/mL), the volume fraction for 100-nm FNDs is ~3·10⁻⁴ — borderline; for 1 µm particles at the same mass concentration, ~3·10⁻⁴ as well, but with fewer larger objects. Working concentrations for tracking experiments will be 10²–10⁴× lower, well within the dilute regime. Aggregation is a separate concern and is not modelled here.

c. **No phase transition.** Water remains liquid throughout. The model is invalid in any regime where ice nucleation, growth, or interface motion is relevant. Those phenomena are the subject of a separate, later breakout.

d. **No solute effects.** DI water assumed throughout. Realistic FND suspensions may contain residual surface chemistry or stabilizers; their effect on η or ρ_f is below the modelling threshold here but should be checked against vendor data sheets.

e. **No optical-trapping or external forces.** If experiments later use optical tweezers for active particle positioning, this model does not apply directly to that case.

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
- all four validation checks pass,
- the parameter scan completes,
- deliverables 1 through 5 exist and are reproducible from a clean clone of the repository.

The breakout terminates unsuccessfully and triggers a re-deliberation if:
- the three numerical methods (A, B, C) cannot be reconciled to within 5 % on any tested cell of the parameter grid,
- the model produces predictions that contradict standard textbook expressions (e.g. Perrin's barometric distribution, Einstein–Smoluchowski),
- effort estimate is exceeded by more than a factor of two without clear scientific justification.

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
