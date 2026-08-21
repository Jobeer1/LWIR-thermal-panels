# Monte Carlo Radiative Exchange Simulator

A Flask web app for modeling radiative heat exchange between two plates using Monte Carlo ray tracing with **thin-film optical depth corrections**. Specifically designed to prove and quantify anisotropic decoupling in micro/nano-structured thermal emitters.

## 🎯 Purpose & Physics

This simulator validates the claim that **anisotropic micro/nano-structures can achieve extreme decoupling** between:
- **Effective absorptivity (α_eff)** → 1.0 for incident radiation
- **Effective emissivity (ε_eff)** << α_eff for internal thermal emission

### Real Physics Modeled:
1. **Thin-Film Optical Depth** (Beer-Lambert law)
   - Walls thinner than absorption depth: ε_eff = ε_bulk × [1 - exp(-t/δ(λ))]
   - 100nm alumina at 300K: ε_eff ≈ 2.2% (vs. 80% bulk assumption)

2. **Waveguide Modal Cutoff** (sub-wavelength confinement)
   - Pores below λ_c cannot sustain propagating modes
   - Internal thermal photons become evanescent waves → decay before escaping
   - Reference: Narayanaswamy & Chen, PRB 2004

3. **Photonic Density of States** (LDOS suppression)
   - Deep cavity emission suppressed by geometry
   - Reference: Lin et al., PRB 2000

4. **Multi-Scale Physics**
   - Micro: Individual photon interactions in cavities
   - Meso: Unit cell periodic behavior
   - Macro: Two-plate radiative equilibrium

## 🧪 Key Features

### Geometry Modes
- **Honeycomb Cavities**: PAA (porous anodized alumina) hexagonal arrays
  - Adjustable: diameter, height, wall thickness, pitch
  - Calculates packing fraction & cavity enhancement factor
- **CNT Forest**: Vertically aligned carbon nanotube arrays
  - Tapered/graded forest geometry
  - Directional anisotropy built-in
- **Frustum Cavities**: Conical/tapered pores
- **Rectangular Pits**: Simple 2D reference geometry

### Material Physics
- **Wavelength-dependent absorption depth** (µm):
  - Alumina (Al₂O₃): δ ≈ 3.61µm at λ=10µm
  - CNT forests: δ ≈ 2.5µm (strong absorption)
  - Silver (Ag): δ ≈ 0.12µm (highly reflective)
- **Thin-film correction factor**: Applied per-photon based on wavelength
- **Spectral emissivity**: Planck-weighted integration over thermal spectrum

### Monte Carlo Engine
- **Photon-by-photon tracing** inside cavity geometry
- **Spectral sampling**: Wavelengths from Planck distribution at plate temperature
- **Lambertian reflection** at each bounce (realistic diffuse scattering)
- **Evanescent decay**: Photons above cutoff wavelength treated as non-escaping
- **Russian Roulette weighting**: Efficient treatment of deep cavities
- **Adaptive convergence**: Statistics monitored during run

### Radiative Exchange Model
- **Four-surface enclosure**: Plate A front, Plate A back, Plate B cavity, Surroundings
- **View factor calculation**: F_AB for two parallel rectangles
- **Radiosity network**: Solves coupled equations for heat exchange
- **Thermal equilibrium**: Verifies net flux → 0 at T_A = T_B

### Output Metrics
| Metric | Meaning | Typical Range |
|--------|---------|---------------|
| α_eff | Effective absorptivity (incident) | 0.8–0.99 |
| ε_B | Effective emissivity (emission) | 0.01–0.5 |
| P_esc | Escape probability | 0.001–0.1 |
| Decoupling ratio (α_eff / ε_B) | Anisotropy measure | 2–100 |
| Net flux Q | Heat from A to B | W/m² |
| T_B_stag | Adiabatic equilibrium temp | K |

## 🏗️ App Architecture

```
┌─────────────────────────────────────────────────────────┐
│  Browser UI (HTML/CSS/JS)                               │
│  ├─ Geometry selector (honeycomb, CNT, frustum, etc.)   │
│  ├─ Material inputs (temperature, emissivity, size)     │
│  ├─ MC settings (photon count, convergence)             │
│  └─ Results display (flux, emissivity, diagnostics)     │
└──────────────────────┬──────────────────────────────────┘
                       │ POST /api/simulate (JSON)
┌──────────────────────▼──────────────────────────────────┐
│  Flask REST API (app.py)                                │
│  ├─ Input validation & sanitization                     │
│  └─ Calls simulator.run_simulation()                    │
└──────────────────────┬──────────────────────────────────┘
                       │ Python simulation
┌──────────────────────▼──────────────────────────────────┐
│  Simulation Engine (simulator.py)                       │
│  ├─ Build cavity geometry (geometry.py)                 │
│  ├─ Run Monte Carlo ray tracer (ray_tracer.py)          │
│  ├─ Apply thin-film corrections (material_optics.py)    │
│  ├─ Solve radiosity enclosure (radiosity model)         │
│  └─ Compile results & statistics                        │
└──────────────────────┬──────────────────────────────────┘
                       │ JSON response
┌──────────────────────▼──────────────────────────────────┐
│  Results → Browser Display                              │
│  ├─ Key metrics (α_eff, ε_B, decoupling ratio)          │
│  ├─ 95% confidence intervals                            │
│  ├─ Diagnostics (cutoff wavelength, evanescent decay)   │
│  └─ Physical validation (Kirchhoff reciprocity check)   │
└─────────────────────────────────────────────────────────┘
```

## 📁 Project Structure

```
monte-carlo-ray-tracing/
├── README.md                           # This file
├── accuracy_improvement_plan.md        # Roadmap for accuracy enhancements
├── simulation_architecture.md          # Detailed physics documentation
├── peer_review.txt                     # Physics citations & validation
│
├── app.py                              # Flask server & /api/simulate route
├── simulator.py                        # Main orchestration layer
├── ray_tracer.py                       # Monte Carlo 3D ray tracing
├── geometry.py                         # Cavity geometry definitions
├── sampling.py                         # Spectral & directional sampling
├── spectral.py                         # Material spectral properties
├── material_optics.py                  # Thin-film physics (Beer-Lambert)
│
├── templates/
│   └── index.html                      # Web UI form & dashboard
├── static/
│   ├── app.js                          # Frontend API calls & rendering
│   └── style.css                       # UI styling (South African flag theme)
│
├── requirements.txt                    # Python dependencies (Flask, NumPy, etc.)
├── .gitignore                          # Git ignore patterns
│
├── test.py                             # Integration tests
├── test_thin_film.py                   # Thin-film physics validation
├── monte_carlo_radiation.py            # Legacy 2D reference model
├── rad_leakage.py                      # 2D educational model
└── LICENSE                             # MIT License
```

## 🚀 Quick Start

### 1) Set Up Python Environment

```bash
python -m venv .venv
```

**Windows (PowerShell):**
```powershell
.\.venv\Scripts\Activate.ps1
```

**macOS/Linux:**
```bash
source .venv/bin/activate
```

### 2) Install Dependencies

```bash
pip install -r requirements.txt
```

### 3) Run the App

```bash
python app.py
```

Open browser:
```
http://127.0.0.1:5000/
```

## 📊 Using the Web Interface

1. **Select Geometry**: Honeycomb, CNT Forest, Frustum, or Rectangular pit
2. **Enter Material Properties**:
   - Temperatures (Plate A hot, Plate B cool, Surroundings)
   - Emissivity/absorptivity values
   - Cavity dimensions & wall thickness
3. **Set MC Parameters**:
   - Photon count (default 20k, increase for accuracy)
   - Enable full-gap view factor check
4. **Click "Run Simulation"**
5. **Review Results**:
   - Effective emissivity & absorptivity with 95% CI
   - Net radiative flux
   - Decoupling ratio (α_eff / ε_B)
   - Cavity enhancement factor
   - Modal cutoff analysis

## 🔬 REST API

### Endpoint
```
POST /api/simulate
Content-Type: application/json
```

### Example Request
```json
{
  "geometry_mode": "honeycomb",
  "height": 20000,
  "cavity_diameter": 500,
  "wall_thickness": 50,
  "temp_a": 600,
  "temp_b": 300,
  "temp_surr": 300,
  "emissivity_a": 0.981,
  "emissivity_a_back": 0.051,
  "gap": 100.001,
  "n_photons": 20000,
  "full_gap_mc": false
}
```

### Example Response
```json
{
  "status": "ok",
  "results": {
    "alpha_eff": 0.9101,
    "epsilon_b": 0.7097,
    "epsilon_b_95ci": 0.1050,
    "p_esc": 0.0072,
    "p_esc_95ci": 0.0011,
    "decoupling_ratio": 1.2857,
    "net_flux_A_front": -6.6,
    "total_net_Q": 0.0,
    "T_B_stagnation_C": 11.6,
    "cutoff_wavelength_um": 853.0,
    "cavity_enhancement": 161.0,
    "view_factor_AB": 0.82699,
    "kirchhoff_error_percent": 0.0,
    "near_field_warning": false
  }
}
```

## 🧬 Physics Implementation

### Thin-Film Beer-Lambert Law
For walls thinner than optical penetration depth:
```
ε_eff(λ) = ε_bulk × [1 - exp(-t/δ(λ))]
```
where:
- t = wall thickness (µm)
- δ(λ) = wavelength-dependent absorption depth (µm)
- ε_bulk = bulk material emissivity

**Example**: 100nm alumina at λ=10µm (300K peak radiation)
- Bulk ε = 0.80
- Absorption depth δ ≈ 3.61µm
- Effective ε ≈ 0.80 × [1 - exp(-0.1/3.61)] ≈ 0.022 (2.2%)
- **Correction factor: 36×**

### Waveguide Cutoff (TE11 mode)
```
λ_c = 1.706 × diameter
```
Photons with λ > λ_c cannot propagate; evanescent decay:
```
δ_ev = (λ_c / 2π) / √(1 - (λ_c/λ)²)
```

### Planck-Weighted Emissivity
Effective emissivity averaged over thermal spectrum:
```
ε_eff(T) = ∫ ε(λ) × M_λ(λ,T) dλ
```
where M_λ is Planck spectral radiance.

## 📈 Expected Accuracy

| Metric | Current | Goal (Future) |
|--------|---------|---------------|
| Wall ε_eff error | < 5% | < 2% |
| Statistical CI | ±10% | ±2% |
| Net flux at equilibrium | Verified | < 0.1% σT⁴ |
| MC runtime | 5–30s | Adaptive, converged |

## 🔍 Validation & Limitations

### Strengths
- ✓ Accurate thin-film physics (Beer-Lambert law)
- ✓ Waveguide modal analysis (cutoff suppression)
- ✓ Spectral integration (Planck-weighted)
- ✓ Multi-scale physics (cavity → enclosure)
- ✓ Energy conservation (net flux → 0 at equilibrium)

### Limitations
- ✗ **Far-field only**: No near-field tunneling (gap < λ/2π triggers warning)
- ✗ **Geometric optics**: Ray-based, not full wave equation
- ✗ **Diffuse assumption**: Assumes Lambertian reflection (valid for rough surfaces)
- ✗ **Periodic idealization**: Real cavities have edge effects
- ✗ **Material data**: Absorption depth from literature (not measured)

### Recommended Validations
1. **Flat plate limit**: Set cavity depth → 0, compare to parallel-plate formula
2. **Blackbody limit**: Set all ε = 1.0, should get α_eff = ε_B = 1.0
3. **Thermal equilibrium**: Run at T_A = T_B, net flux should → 0
4. **Literature comparison**: Match published results for known geometries

## 📚 References

### Thin-Film Physics
- Born & Wolf, *Principles of Optics* (1999) — Beer-Lambert law
- Palik, *Handbook of Optical Constants of Solids* (1998) — Absorption depth data

### Anisotropic Decoupling
- Lin et al., *Phys. Rev. B*, vol. 62, pp. 3081–3084 (2000) — LDOS suppression
- Narayanaswamy & Chen, *Phys. Rev. B*, vol. 70, p. 125101 (2004) — Waveguide cutoff
- Sprafke et al., *Adv. Opt. Mater.*, vol. 1, pp. 527–535 (2013) — PAA light trapping
- Mizuno et al., *Proc. Natl. Acad. Sci. USA*, vol. 106, pp. 6044–6047 (2009) — CNT forests

### Monte Carlo Methods
- Howell et al., *Thermal Radiation Heat Transfer* (5th ed.) — View factors & ray tracing

## 💬 Feedback & Issues

Found a bug? Have a suggestion? Open an issue:

```
https://github.com/Jobeer1/LWIR-thermal-panels/issues
```

Include:
- Simulation parameters
- Expected vs. actual output
- Python version & OS
- Error logs (if any)

## 📄 License

MIT License — See [LICENSE](LICENSE) file for details.

---

**Last Updated**: August 2026  
**Physics Model**: Anisotropic radiative exchange with thin-film corrections  
**Status**: Active development with accuracy improvements in progress
