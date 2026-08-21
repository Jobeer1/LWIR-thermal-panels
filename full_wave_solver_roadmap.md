# Full-Wave Solver Roadmap

## Purpose

This document describes how to extend the current Monte Carlo radiative-exchange app with a full electromagnetic solver for micro- and nano-structured cavities.

The goal is not to replace the entire application with an electromagnetic time-domain calculation. The practical target is a multiscale workflow:

1. Solve one representative 3-D cavity or periodic unit cell with Maxwell's equations.
2. Extract wavelength-, angle-, polarization-, and temperature-dependent optical data.
3. Convert those data into effective absorptivity, emissivity, modal transmission, and leakage quantities.
4. Feed the reduced results into the existing Monte Carlo and radiosity models.

Two solver families are considered:

- **FDTD**: general-purpose, broadband, time-domain Maxwell solver.
- **Coupled-Mode Theory (CMT)**: reduced-order frequency-domain model for a known set of guided, leaky, or resonant modes.

For this project, the recommended architecture is a **hybrid FDTD plus CMT workflow**. FDTD provides trustworthy reference data for selected geometries. CMT or a fitted modal model provides the fast model used during interactive parameter sweeps.

---

## 1. Why the Current Model Needs a Full-Wave Path

The current application is a physics-informed radiative model. Its main components are:

| Current component | Current responsibility | Full-wave replacement or input |
|---|---|---|
| `geometry.py` | Builds honeycomb, CNT, frustum, and rectangular cells | Creates electromagnetic meshes and periodic boundaries |
| `material_optics.py` | Supplies absorption depths and thin-film emissivity | Supplies complex `n(lambda, T)` and multilayer response |
| `ray_tracer.py` | Estimates internal escape and external absorption | Consumes modal transmission, angle-dependent absorption, and leakage tables |
| `simulator.py` | Applies unit-cell results to the macro enclosure | Interpolates full-wave data and solves radiosity/heat exchange |
| `sampling.py` | Samples Planck wavelengths | Supplies spectral weighting for full-wave output integration |

The existing approximations are useful in the geometric-optics regime, but they are insufficient when any of the following are important:

- cavity diameter, gap, wall thickness, or feature pitch is comparable to the wavelength;
- interference between front and back interfaces matters;
- the aperture supports only a small number of modes;
- evanescent fields tunnel across a short distance;
- resonances, guided modes, surface waves, or polarization conversion matter;
- the structure is periodic and diffracts into non-specular orders;
- local density of states (LDOS) or near-field heat transfer is required.

A full-wave solver should therefore replace the current sharp binary cutoff and empirical top-surface capture rule in the regimes where those effects control the result.

---

## 2. Recommended Overall Architecture

```mermaid
flowchart LR
    A[Geometry parameters] --> B[Unit-cell builder]
    B --> C[Material model n(lambda,T)]
    C --> D[FDTD reference solver]
    C --> E[CMT / modal solver]
    D --> F[Reference S parameters and fields]
    D --> G[Mode extraction and calibration]
    E --> H[Fast modal transmission model]
    F --> H
    H --> I[Planck and angular integration]
    I --> J[Effective alpha, epsilon, leakage, LDOS]
    J --> K[Monte Carlo / radiosity application]
```

### 2.1 Division of responsibilities

**Full-wave layer**

- solves Maxwell's equations in a finite 3-D unit cell;
- handles material dispersion, interfaces, diffraction, polarization, and evanescent fields;
- returns complex fields, power flux, S parameters, absorption, and modal content.

**Reduced-order layer**

- represents the full-wave result as a small set of modes or fitted response functions;
- evaluates transmission and absorption quickly for many temperatures and geometries;
- provides a stable interface to the existing application.

**Radiative layer**

- performs Planck weighting;
- maps unit-cell results to panel-scale properties;
- retains the current view-factor and radiosity network where far-field exchange is appropriate.

---

## 3. Choose the Solver Mode

### 3.1 FDTD: when to use it

Use FDTD when the geometry or physics is not yet understood well enough to choose a modal basis. FDTD is appropriate for:

- arbitrary 3-D cavity shapes;
- tapered CNT or pore geometry;
- broadband spectra from one impulse or short pulse;
- strong material contrast;
- multiple interacting modes;
- diffraction and polarization conversion;
- transient field buildup and decay;
- extracting reference data for CMT fitting.

FDTD is expensive. A realistic 3-D simulation may require millions to billions of Yee cells and substantial memory. It should not run inside every Flask request.

### 3.2 CMT: when to use it

Use CMT after the dominant modes and loss channels are known. CMT is appropriate for:

- a narrow spectral band;
- long, nearly uniform channels;
- weak taper or weak perturbations;
- a small number of propagating and evanescent modes;
- fast sweeps over height, loss, temperature, or aperture coupling;
- fitting FDTD results into an interactive engineering model.

CMT is not a substitute for discovering unknown modes. A bad or incomplete modal basis can produce confidently wrong results.

### 3.3 Recommended decision

Implement in this order:

1. **FDTD reference solver for one honeycomb cavity.**
2. **FDTD reference solver for one CNT periodic cell.**
3. **Mode extraction from the FDTD fields.**
4. **CMT transmission and loss model.**
5. **Use CMT for interactive scans; retain FDTD as an offline validation tool.**

Do not begin by attempting a general-purpose FDTD implementation from scratch in this repository. Use a validated electromagnetic solver library or an external solver and build an adapter around it.

---

## 4. Maxwell Equations and Conventions

For linear, isotropic, non-magnetic media, solve:

$$
\nabla \times \mathbf{E} = -\mu_0 \frac{\partial \mathbf{H}}{\partial t}
$$

$$
\nabla \times \mathbf{H} = \epsilon_0 \epsilon_r(\mathbf{r}, \omega) \frac{\partial \mathbf{E}}{\partial t} + \mathbf{J}
$$

For dispersive and lossy materials, use a complex refractive index:

$$
\tilde{n}(\lambda,T) = n(\lambda,T) + i k(\lambda,T)
$$

The sign of the imaginary part depends on the time convention. The implementation must choose one convention and use it consistently in:

- material files;
- FDTD dispersive equations;
- frequency-domain phasors;
- Poynting-vector calculations;
- TMM and CMT adapters.

The time-averaged power flux is:

$$
\langle \mathbf{S} \rangle = \frac{1}{2} \operatorname{Re}(\mathbf{E} \times \mathbf{H}^*)
$$

For a source-free volume, absorbed power density is:

$$
P_{abs}(\mathbf{r},\omega) = \frac{\omega \epsilon_0}{2} \operatorname{Im}(\epsilon_r) |\mathbf{E}|^2
$$

These definitions should be centralized in a new `wave_physics/conventions.py` module and tested against a simple homogeneous medium.

---

## 5. FDTD Implementation Plan

### 5.1 Proposed module layout

```text
wave_physics/
    __init__.py
    conventions.py
    materials.py
    mesh.py
    boundaries.py
    sources.py
    monitors.py
    fdtd_solver.py
    mode_decomposition.py
    s_parameters.py
    postprocess.py
    cache.py
    schemas.py
cmt/
    __init__.py
    modes.py
    coupled_mode_solver.py
    fitting.py
    reduced_model.py
wave_solver_service.py
wave_solver_cli.py
wave_solver_validation.py
```

### 5.2 Geometry-to-mesh adapter

Add a mesh adapter for each existing geometry class:

```python
class ElectromagneticCell:
    points: object
    material_id: object
    aperture_plane: object
    source_planes: object
    monitor_planes: object
    periodic_axes: tuple[str, ...]
```

Required adapters:

- `HoneycombCavityCell -> cylindrical pore in a periodic hexagonal or equivalent cell`;
- `CNTForestCell -> periodic square cell containing a tapered solid CNT`;
- `FrustumCavity3D -> isolated or periodic tapered pore`;
- `RectPit3D -> rectangular reference waveguide`.

The electromagnetic cell must explicitly represent:

- the actual solid wall thickness;
- the cavity void;
- the substrate/base;
- the exterior medium;
- periodic neighboring material;
- the top aperture and bottom termination.

Do not infer wall thickness from a ray-tracing boundary after meshing. A full-wave solver needs the material volume itself.

### 5.3 Spatial discretization

For a basic Yee-grid FDTD solver, choose the smallest wavelength inside the highest-index material:

$$
\Delta x,\Delta y,\Delta z \leq \frac{\lambda_0}{N n_{max}}
$$

where `N` is commonly 10-20 cells per wavelength for initial studies. Increase this for:

- high-index contrast;
- thin walls;
- curved surfaces;
- strongly dispersive metals;
- phase-accurate Q-factor calculations.

The mesh must resolve both the wavelength and the thinnest physical feature. If the wall is thinner than the required cell size, use subpixel smoothing, conformal meshing, or a frequency-domain finite-element method instead of pretending the wall is resolved.

Perform a mesh convergence study by repeating a reference case with progressively smaller cells and requiring stability in:

- total absorption;
- reflected and transmitted power;
- resonance wavelength;
- Q factor;
- aperture leakage.

### 5.4 Time step and stability

For a 3-D Cartesian Yee grid, obey the Courant condition:

$$
\Delta t \leq \frac{S}{c\sqrt{\Delta x^{-2}+\Delta y^{-2}+\Delta z^{-2}}}
$$

with safety factor `S < 1`, normally around `0.9` or lower depending on the implementation.

The solver must stop with a clear validation error if the requested grid and time step violate stability.

### 5.5 Boundary conditions

Use boundaries according to the physical experiment:

| Boundary | Use |
|---|---|
| Periodic / Bloch | Infinite periodic pore or CNT array |
| PML | Open exterior and absorbing axial ends |
| PEC / PMC | Symmetry planes only when the source polarization and geometry allow it |
| Bloch phase | Oblique incidence with transverse wavevector `k_parallel` |

For a periodic unit cell, the Bloch condition is:

$$
\mathbf{F}(\mathbf{r}+\mathbf{a}) = e^{i\mathbf{k}_{\parallel}\cdot\mathbf{a}}\mathbf{F}(\mathbf{r})
$$

A normal-incidence simulation with zero phase is not enough to characterize hemispherical absorption. At minimum, sample incident angle and both independent polarizations.

PML validation is mandatory. Confirm that reflected power from a homogeneous wave does not change the result when the PML thickness is increased.

### 5.6 Source and monitor setup

Each simulation should define:

- source wavelength or broadband pulse;
- incident angle and azimuth;
- polarization: TE, TM, or arbitrary;
- source plane above the aperture;
- reflected-power monitor above the source;
- transmitted-power monitor below or inside the cavity;
- absorption monitor throughout the walls and base;
- near-field monitors for mode extraction;
- optional thermal dipole or fluctuational source for LDOS studies.

For an incident plane wave, calculate:

$$
R = \frac{P_{ref}}{P_{inc}}, \qquad
T = \frac{P_{trans}}{P_{inc}}, \qquad
A = 1 - R - T
$$

The numerical residual `|1 - R - T - A|` must be reported. A result with poor energy balance must not be sent to the production model.

### 5.7 Broadband strategy

Use one of these approaches:

1. broadband temporal pulse plus Fourier transform;
2. independent frequency-domain runs;
3. a hybrid approach: broadband FDTD for smooth regions and narrow frequency scans near resonances.

For lossy and dispersive materials, independent frequency runs are often easier to validate. Broadband FDTD is efficient but requires careful source normalization, sampling, and material-pole treatment.

---

## 6. Material Model Requirements

The current absorption-depth tables are not sufficient for a quantitative Maxwell solver. Replace or extend them with measured or validated optical constants:

```json
{
  "material": "alumina",
  "temperature_K": 300,
  "wavelength_um": [2, 5, 8, 10, 12, 15, 20],
  "n": [1.7, 1.7, 1.7, 1.7, 1.7, 1.7, 1.7],
  "k": [0.01, 0.015, 0.02, 0.03, 0.04, 0.05, 0.06],
  "source": "measured dataset or publication identifier"
}
```

Requirements:

- preserve source and uncertainty metadata;
- interpolate `n` and `k` without producing negative absorption;
- enforce causality-compatible dispersion where possible;
- support temperature-dependent data;
- support anisotropic tensors for CNT forests when needed;
- support Drude, Lorentz, and Debye poles for time-domain FDTD;
- distinguish bulk material from effective-medium CNT forest material.

For a first implementation, use isotropic material data and explicitly mark CNT anisotropy as a later phase. Do not present scalar CNT constants as experimentally validated tensor data.

### 6.1 FDTD dispersive materials

FDTD cannot directly use an arbitrary complex value independently at every frequency. Fit the measured response to causal poles, for example:

- Drude model for free-carrier metals;
- Lorentz oscillators for resonances;
- Debye terms for low-frequency relaxation.

Validate the fitted time-domain material against the original `n(lambda)` and `k(lambda)` table before using it in a cavity simulation.

### 6.2 Thermal emission consistency

For a reciprocal structure at thermal equilibrium, absorption and emission must be generated from the same electromagnetic response. Do not combine an FDTD absorptivity curve with an unrelated empirical emission gate without documenting the approximation.

For far-field emission, use Kirchhoff's law per angle, polarization, and wavelength where reciprocity applies:

$$
\epsilon(\lambda, \theta, \phi, p) = \alpha(\lambda, \theta, \phi, p)
$$

For near-field emission and LDOS, use a fluctuational-electrodynamics or Green-function calculation. Far-field absorptivity alone does not determine near-field heat transfer.

---

## 7. CMT / Modal Solver Plan

### 7.1 Modal basis

For each uniform or slowly varying channel section, solve for transverse modes at each frequency:

$$
\mathbf{E}(x,y,z) = \sum_m a_m(z)\mathbf{e}_m(x,y)e^{i\beta_m z}
$$

The propagation constant is complex:

$$
\beta_m = \beta'_m + i\beta''_m
$$

Power attenuation in a uniform mode is approximately:

$$
P_m(z) = P_m(0)e^{-2\beta''_m z}
$$

This is the modal origin of the `exp(-2 kappa h)` tunneling factor for an evanescent mode.

### 7.2 Coupled-mode equations

For a slowly varying or perturbed structure, solve:

$$
\frac{d a_m}{dz} = i\beta_m a_m + \sum_n \kappa_{mn}(z)a_n - \gamma_m a_m + s_m(z)
$$

where:

- `a_m` is modal amplitude;
- `beta_m` is the propagation constant;
- `kappa_mn` is mode coupling;
- `gamma_m` represents material or radiation loss;
- `s_m` represents source coupling.

The implementation should retain both forward and backward waves when reflections or resonances matter:

$$
\frac{d}{dz}
\begin{bmatrix}
\mathbf{a}^{+} \\
\mathbf{a}^{-}
\end{bmatrix}
= \mathbf{M}(z,\omega)
\begin{bmatrix}
\mathbf{a}^{+} \\
\mathbf{a}^{-}
\end{bmatrix}
+ \mathbf{s}(z,\omega)
$$

### 7.3 Mode normalization

Normalize modes by power flux, not arbitrary field amplitude:

$$
P_m = \frac{1}{2}\operatorname{Re}\int_A
(\mathbf{e}_m \times \mathbf{h}_m^*)\cdot \hat{z}\,dA
$$

Evanescent modes carry no net power in a uniform region by themselves, but they are essential for aperture matching and near-field coupling. Keep their amplitude in the basis even when their standalone axial power is zero.

### 7.4 Aperture and interface matching

At each discontinuity, match tangential electric and magnetic fields. The resulting scattering matrix should include:

- propagating input/output modes;
- evanescent near-field modes;
- reflection at the aperture;
- absorption in the walls;
- mode conversion between polarizations and symmetry classes.

For a uniform channel of length `h`, propagation is represented by:

$$
\mathbf{P}(h) = \operatorname{diag}(e^{i\beta_m h})
$$

For evanescent modes where `beta_m = i kappa_m`, the field factor is `e^{-kappa_m h}` and the power-like transmission contribution scales as `e^{-2 kappa_m h}`.

Use S-matrix or transfer-matrix cascading. S-matrix cascading is usually more numerically stable for long, lossy, or strongly evanescent sections.

### 7.5 Extracting CMT parameters from FDTD

For each geometry and frequency:

1. Run a mode solve or FDTD excitation.
2. Project the field onto candidate transverse modes.
3. Fit `beta_m`, attenuation, and interface coupling from several axial planes.
4. Calculate the overlap matrix between adjacent sections.
5. Fit `kappa_mn` for taper or periodic perturbation.
6. Compare CMT S parameters against FDTD S parameters.
7. Increase the modal basis until the error target is met.

The reduced model should report its truncation error and the frequency range over which it is valid.

---

## 8. Extracting the Quantities the Current App Needs

The current simulator should receive a versioned full-wave result object rather than raw field arrays.

```python
@dataclass
class WaveResponse:
    wavelength_um: np.ndarray
    theta_rad: np.ndarray
    phi_rad: np.ndarray
    polarization: list[str]
    reflectance: np.ndarray
    transmittance: np.ndarray
    absorptance: np.ndarray
    aperture_leakage: np.ndarray
    modal_transmission: np.ndarray
    energy_balance_error: np.ndarray
    solver_kind: str
    geometry_hash: str
    material_hash: str
```

### 8.1 External absorptivity

For each incident wavelength, angle, azimuth, and polarization:

$$
\alpha_{eff}(\lambda,\theta,\phi,p) = A(\lambda,\theta,\phi,p)
$$

Integrate over the incident source distribution and Planck spectrum:

$$
\bar{\alpha}(T) =
\frac{\int M_\lambda(T)\alpha(\lambda,\theta,\phi,p)\,d\lambda\,d\Omega}
{\int M_\lambda(T)\,d\lambda\,d\Omega}
$$

The angular distribution must match the physical source. A normal-incidence result cannot be silently reused as a hemispherical absorptivity.

### 8.2 Far-field emissivity

For a reciprocal far-field calculation, integrate directional absorptivity using the same angular and polarization weighting:

$$
\bar{\epsilon}(T) =
\frac{\int M_\lambda(T)\epsilon(\lambda,\theta,\phi,p)\,d\lambda\,d\Omega}
{\int M_\lambda(T)\,d\lambda\,d\Omega}
$$

This replaces the current empirical `g_em` gate where full-wave data are available.

### 8.3 Internal leakage

For emission from a specified internal source region, calculate the fraction of source power crossing the aperture monitor:

$$
L(\lambda,T) =
\frac{P_{aperture}(\lambda,T)}{P_{source}(\lambda,T)}
$$

Use different source classes when needed:

- wall thermal noise sources;
- base thermal noise sources;
- localized dipoles at different depths;
- volume-distributed fluctuational sources.

A single source at the cavity center is not sufficient for a depth-dependent leakage model.

### 8.4 LDOS and near-field transfer

If the goal includes true sub-wavelength gap heat transfer, add a separate fluctuational-electrodynamics path. Required outputs may include:

- electric and magnetic Green tensors;
- electric LDOS;
- spectral heat-transfer coefficient;
- evanescent wavevector integration over `k_parallel`;
- surface-mode contributions.

Do not label the far-field FDTD aperture leakage as a near-field heat-transfer result.

---

## 9. How to Integrate with the Existing App

### 9.1 New service boundary

Add a service with two modes:

```python
response = wave_solver.solve(
    geometry=geometry,
    wavelengths_um=wavelength_grid,
    angles=angle_grid,
    polarization=('TE', 'TM'),
    solver='cached_cmt',
)
```

Suggested implementations:

- `FDTDWaveSolver`: offline or queued reference runs;
- `CMTWaveSolver`: fast local evaluation;
- `CachedWaveSolver`: loads validated response tables;
- `FallbackRaySolver`: current Monte Carlo approximation when no full-wave result exists.

### 9.2 Simulator behavior

Add a configuration option:

```python
wave_model = 'ray'       # existing approximation
wave_model = 'cmt'       # fast reduced-order model
wave_model = 'fdtd_cache' # validated precomputed data
```

The default should remain `ray` until the new solver passes validation. When a wave response is available:

1. sample or quadrature-integrate wavelength and angle;
2. compute `alpha_eff` from full-wave absorptance;
3. compute `epsilon_b` from reciprocal emission or explicit thermal-source leakage;
4. report modal and evanescent contributions;
5. send only effective properties to the radiosity layer.

The radiosity equations should not know whether the effective properties came from rays, CMT, or FDTD.

### 9.3 Cache key

Cache responses using a hash of:

- geometry type and all dimensions;
- mesh resolution;
- boundary conditions;
- material dataset and temperature;
- wavelength/angle grid;
- polarization;
- solver version;
- source and monitor configuration.

Never reuse a response when its material, geometry, or boundary-condition hash differs.

### 9.4 API and UI

Expose:

- selected wave model;
- solver status: `fallback`, `cached`, `cmt`, or `fdtd`;
- wavelength range and angular coverage;
- energy-balance error;
- modal truncation error;
- mesh convergence status;
- whether near-field effects are included;
- data provenance and material source.

Avoid presenting a CMT or cached FDTD result as an exact solution if it is only valid over a restricted band or angle range.

---

## 10. Validation Program

### 10.1 Unit tests

Implement analytic tests before testing complex geometries:

1. homogeneous medium: numerical wave speed and impedance;
2. normal-incidence Fresnel reflection;
3. lossless slab: `R + T = 1`;
4. lossy slab: `R + T + A = 1`;
5. PEC reflection;
6. PML reflection versus PML thickness;
7. periodic phase shift;
8. rectangular waveguide cutoff;
9. circular waveguide cutoff;
10. evanescent decay versus analytic `exp(-kappa h)` field amplitude.

### 10.2 Geometry benchmarks

Use the existing geometry classes and compare against known limits:

| Case | Expected behavior |
|---|---|
| Zero-height cavity | Approaches a flat interface or thin-film result |
| Very large aperture | Leakage approaches the open-surface response |
| Perfectly black walls | Absorption is unity after sufficient path length |
| Lossless walls | No material absorption; power exits or remains stored |
| Uniform rectangular channel | Matches analytic TE/TM mode propagation |
| Long sub-cutoff channel | Transmission decreases as `exp(-2 kappa h)` |
| Symmetric reciprocal cell | Forward and reverse S parameters obey reciprocity |

### 10.3 Conservation and reciprocity checks

Every production response must report:

$$
\delta_E = |1 - R - T - A|
$$

For a reciprocal passive structure, verify the appropriate S-matrix symmetry. For an isothermal reciprocal system, verify that the integrated net radiative exchange is zero within numerical tolerance.

### 10.4 FDTD-to-CMT acceptance criteria

For the chosen calibration geometries, require for example:

- absolute absorption error below 1-2% over the target band;
- resonance wavelength error below 1%;
- leakage error below 5% in the target aspect-ratio range;
- energy-balance residual below `1e-3` or a documented frequency-dependent threshold;
- stable result under mesh, time-step, PML, and modal-basis refinement.

These are engineering acceptance criteria, not universal physical guarantees. Set them per application and record them with the cached response.

---

## 11. Performance and Deployment

### FDTD

- run offline, in a worker process, or on a separate compute service;
- use MPI/GPU acceleration only after correctness is established;
- store fields selectively rather than retaining every time step;
- use checkpointing for long jobs;
- use job IDs and status polling from Flask;
- never block a web request for a large 3-D solve.

### CMT

- run synchronously for small modal systems;
- cache mode profiles and fitted parameters;
- use vectorized wavelength and angle evaluation;
- expose a validity interval and error estimate;
- fall back to FDTD cache or ray tracing outside the validity interval.

### Data storage

Store response data in a versioned format such as HDF5 or Zarr:

```text
response.h5
  /metadata
  /wavelength_um
  /angles/theta_rad
  /angles/phi_rad
  /polarization/TE/reflectance
  /polarization/TE/absorptance
  /polarization/TE/aperture_leakage
  /polarization/TM/...
  /modes/beta
  /validation/energy_balance_error
```

For lightweight deployments, store reduced CMT parameters and interpolate them in memory.

---

## 12. Phased Implementation Schedule

### Phase 0: Define contracts

- add `WaveResponse` schema;
- define unit, sign, and polarization conventions;
- add energy-balance and provenance fields;
- add `wave_model='ray'` configuration without changing behavior.

**Exit criterion:** existing application tests remain green and can consume a placeholder response object.

### Phase 1: Analytic mode and interface benchmarks

- implement Fresnel and rectangular/circular waveguide reference calculations;
- test propagation and evanescent decay;
- test S-parameter normalization;
- define material data schema.

**Exit criterion:** all analytic benchmark tests pass.

### Phase 2: First FDTD reference cell

- implement only the honeycomb cylindrical pore;
- use a lossless dielectric first;
- add lossy wall and substrate;
- add normal-incidence TE/TM sources;
- add aperture and absorption monitors;
- perform mesh, PML, and time-step convergence.

**Exit criterion:** energy conservation and Fresnel/open-channel limits pass.

### Phase 3: Angular, polarization, and periodic sweeps

- add Bloch phase boundaries;
- sample angle and azimuth;
- calculate both polarizations;
- produce angle-resolved `R`, `T`, `A`, and leakage tables.

**Exit criterion:** response tables are reproducible and have bounded interpolation error.

### Phase 4: CNT and tapered structures

- add periodic CNT cell;
- add taper and anisotropic material option;
- compare full-wave results with current cutoff and ray model;
- quantify where ray tracing fails.

**Exit criterion:** documented validity map for ray, CMT, and FDTD models.

### Phase 5: CMT reduction

- extract modes from FDTD;
- fit propagation constants and coupling matrices;
- implement stable S-matrix cascading;
- compare CMT against FDTD over the target parameter range.

**Exit criterion:** CMT meets the defined absorption and leakage error targets.

### Phase 6: Application integration

- add cached response loader;
- add CMT runtime path;
- preserve ray fallback;
- add solver status and validity diagnostics to the API/UI;
- update documentation and examples.

**Exit criterion:** the same input can be run with ray, cached FDTD, and CMT modes, with transparent provenance.

---

## 13. What Not to Do

- Do not call the current empirical cutoff rule a full Maxwell solution.
- Do not use a normal-incidence response as a hemispherical result without angular integration.
- Do not use absorption-depth data as if they were complete complex refractive-index data.
- Do not run a large FDTD mesh in the Flask request thread.
- Do not discard evanescent modes merely because they carry no net power in a uniform section.
- Do not fit CMT parameters from one field plane or one wavelength.
- Do not accept results without energy conservation and convergence checks.
- Do not mix FDTD absorptivity with an unrelated emission gate at thermal equilibrium.
- Do not claim near-field heat transfer from far-field aperture leakage alone.

---

## 14. First Concrete Milestone

The smallest useful full-wave milestone is:

1. Implement a validated rectangular waveguide solver or use a trusted external FDTD/FEM package.
2. Simulate a `10 um x 10 um x 450 um` rectangular cell with a lossy wall and aperture.
3. Sweep wavelengths around the predicted cutoff.
4. Measure reflected power, aperture transmission, wall absorption, and axial field decay.
5. Compare the measured sub-cutoff transmission against:

$$
T(h,\lambda) \propto e^{-2\kappa(\lambda)h}
$$

6. Store the response as a versioned cache.
7. Add a `cached_fdtd` path to `ray_tracer.py` or a new wave-response adapter.

This milestone directly tests the peer-review concern while keeping the scope manageable. After it works, move to the cylindrical honeycomb cell and then the CNT forest.

---

## 15. Definition of Done

The full-wave upgrade is complete when:

- the electromagnetic solver uses explicit complex material data;
- geometry, mesh, boundaries, and source conditions are reproducible;
- `R`, `T`, and `A` satisfy energy conservation;
- mesh, PML, time-step, and modal-basis convergence are documented;
- angular and polarization dependence are included;
- sub-cutoff transmission follows the full modal solution rather than a binary cutoff;
- CMT agrees with FDTD over its declared validity range;
- the existing radiosity model consumes wave-derived effective properties through a stable interface;
- ray tracing remains available as a clearly labeled fallback;
- the UI and API display solver provenance, validity, and residual errors;
- near-field results are clearly separated from far-field emissivity results.
