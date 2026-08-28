# Monte Carlo Radiative Heat Transfer — Architecture

This document explains how the simulator works in plain language, with visual
diagrams for quick understanding.

## 1. What the simulator answers

| Question                                         | Output                                     |
| ------------------------------------------------ | ------------------------------------------ |
| How easily do thermal photons escape the cavity? | **Escape Probability** P_esc         |
| How much incoming radiation is absorbed?         | **Effective Absorptivity** alpha_eff |
| How much does the structured surface emit?       | **Effective Emissivity** epsilon_B   |
| How much net heat flows between plates?          | **Net Radiative Flux** q             |
| What temperature does Plate B reach if isolated? | **Stagnation Temperature** T_stag    |

---

## 2. High-level data flow

```mermaid
flowchart TD
    %% Class Definitions
    classDef default fill:#1e293b,stroke:#475569,stroke-width:1.5px,color:#f8fafc,font-family:sans-serif;
    classDef inputNode fill:#0f172a,stroke:#38bdf8,stroke-width:1.5px,color:#e0f2fe;
    classDef orchestrator fill:#0f172a,stroke:#38bdf8,stroke-width:2.5px,color:#38bdf8,font-weight:bold;
    classDef rayMode fill:#064e3b,stroke:#34d399,stroke-width:1.5px,color:#ecfdf5;
    classDef fullWaveMode fill:#3730a3,stroke:#818cf8,stroke-width:1.5px,color:#e0e7ff;
    classDef emtMode fill:#701a75,stroke:#f0abfc,stroke-width:1.5px,color:#fdf4ff;
    classDef nearFieldMode fill:#7f1d1d,stroke:#f87171,stroke-width:1.5px,color:#fef2f2;
    classDef outputNode fill:#065f46,stroke:#10b981,stroke-width:2px,color:#ffffff,font-weight:bold;

    %% Subgraph Structure
    subgraph ENGINE["🚀 Simulator Pipeline Architecture"]
        style ENGINE fill:#0f172a,stroke:#0d9488,stroke-width:2px,color:#2dd4bf,font-weight:bold;
      
        A["User enters parameters<br/>(geometry, temperature, materials)"] --> B["Flask API<br/>app.py"]
        B --> C["Build cavity geometry<br/>geometry.py"]
        C --> D["Physics Orchestrator<br/>regime_selector.py"]

        %% Routed Engine Solvers
        D -->|RAY| E["Monte Carlo ray tracer<br/>ray_tracer.py"]
        D -->|FULL_WAVE| F["Waveguide modal solver<br/>waveguide_modes.py"]
        D -->|EMT| G["Maxwell-Garnett TMM slab<br/>material_optics.py"]
        D -->|NEAR_FIELD| H["Polder-Van Hove<br/>near_field_radiative_heat.py"]

        %% Aggregation & Exitance
        E --> I["Effective emissivity ε_B<br/>Effective absorptivity α_eff"]
        F --> I
        G --> I
        H --> I

        I --> J["4-surface radiosity<br/>Net flux, stagnation temperature"]
        J --> K["Results JSON<br/>+ orchestrator provenance"]
    end

    %% Apply Classes to Nodes
    class A,B,C inputNode;
    class D orchestrator;
    class E rayMode;
    class F fullWaveMode;
    class G emtMode;
    class H nearFieldMode;
    class I,J,K outputNode;

    %% Link Styling
    linkStyle 4 stroke:#34d399,stroke-width:2px,fill:none;
    linkStyle 5 stroke:#818cf8,stroke-width:2px,fill:none;
    linkStyle 6 stroke:#f0abfc,stroke-width:2px,fill:none;
    linkStyle 7 stroke:#f87171,stroke-width:2px,fill:none;
```

---

## 3. Physics Orchestrator: regime selection

Before any solver runs, the orchestrator computes dimensionless ratios:

| Ratio                        | Meaning                              | Threshold                                     |
| ---------------------------- | ------------------------------------ | --------------------------------------------- |
| **lambda/D**           | Wavelength vs. feature diameter      | < 0.2 → RAY, 0.2–5 → FULL_WAVE, > 5 → EMT |
| **lambda/P**           | Wavelength vs. pitch                 | EMT when both D, P ≪ lambda                  |
| **d_gap/(lambda/2pi)** | Gap vs. evanescent tunnelling length | < 1 → NEAR_FIELD                             |
| **t_wall/delta**       | Wall thickness vs. absorption depth  | Optically thin if < 0.5                       |

```mermaid
flowchart LR
    %% Class Definitions
    classDef default fill:#1e293b,stroke:#475569,stroke-width:1.5px,color:#f8fafc,font-family:sans-serif;
    classDef decision fill:#0f172a,stroke:#38bdf8,stroke-width:2px,color:#38bdf8,font-weight:bold;
    classDef rayMode fill:#064e3b,stroke:#34d399,stroke-width:1.5px,color:#ecfdf5;
    classDef fullWaveMode fill:#3730a3,stroke:#818cf8,stroke-width:1.5px,color:#e0e7ff;
    classDef emtMode fill:#701a75,stroke:#f0abfc,stroke-width:1.5px,color:#fdf4ff;
    classDef nearFieldMode fill:#7f1d1d,stroke:#f87171,stroke-width:1.5px,color:#fef2f2;

    %% Main Flowchart
    subgraph ORCH["⚙️ Physics Regime Orchestrator"]
        style ORCH fill:#0f172a,stroke:#0d9488,stroke-width:2px,color:#2dd4bf,font-weight:bold;
      
        T["Temperature T<br/>λ_peak = 2898 / T"] --> R["Compute Ratios<br/>λ/D, λ/P, d/λ, t/δ"]
        R --> D{"Which Regime?"}
      
        D -->|λ/D < 0.2| RAY["RAY<br/>Geometric optics<br/>Monte Carlo ray tracer"]
        D -->|0.2 ≤ λ/D ≤ 5| FW["FULL_WAVE<br/>Diffraction / resonance<br/>Modal solver or RCWA cache"]
        D -->|λ/D > 5| EMT["EMT<br/>Maxwell-Garnett<br/>Homogenised slab"]
        D -->|gap < λ / 2π| NF["NEAR_FIELD<br/>Polder-Van Hove<br/>Green-tensor LDOS"]
    end

    %% Apply Classes to Nodes
    class D decision;
    class RAY rayMode;
    class FW fullWaveMode;
    class EMT emtMode;
    class NF nearFieldMode;

    %% Link Styling
    linkStyle 2 stroke:#34d399,stroke-width:2px,fill:none;
    linkStyle 3 stroke:#818cf8,stroke-width:2px,fill:none;
    linkStyle 4 stroke:#f0abfc,stroke-width:2px,fill:none;
    linkStyle 5 stroke:#f87171,stroke-width:2px,fill:none;
```

The orchestrator also computes a **confidence score** (0–100%) from regime
penalties, temperature drift, and material extrapolation status. This is
surfaced in the UI as a color-coded badge with explicit warnings.

## 4. Monte Carlo cavity experiments

Two photon-scale experiments run inside the repeating unit cell:

```mermaid
flowchart LR
    %% Class Definitions
    classDef default fill:#1e293b,stroke:#475569,stroke-width:1.5px,color:#f8fafc,font-family:sans-serif;
    classDef decision fill:#0f172a,stroke:#38bdf8,stroke-width:2px,color:#38bdf8,font-weight:bold;
    classDef guided fill:#064e3b,stroke:#34d399,stroke-width:1.5px,color:#ecfdf5;
    classDef evanescent fill:#7f1d1d,stroke:#f87171,stroke-width:1.5px,color:#fef2f2;
    classDef emissionOut fill:#065f46,stroke:#10b981,stroke-width:2px,color:#ffffff,font-weight:bold;
    classDef incidenceOut fill:#3730a3,stroke:#818cf8,stroke-width:2px,color:#ffffff,font-weight:bold;

    %% Subgraph 1: Internal Emission
    subgraph E1["🔥 Internal Emission"]
        style E1 fill:#0f172a,stroke:#0d9488,stroke-width:2px,color:#2dd4bf,font-weight:bold;
      
        A1["Photons launched from walls + base"] --> B1["Planck-sampled wavelength"]
        B1 --> C1{"λ < λ_c ?"}
        C1 -->|Yes| D1["Guided mode<br/>T = T_ap × exp(-αL)"]
        C1 -->|No| D2["Evanescent<br/>exp(-2L/δ_ev)"]
        D1 --> OUT1["P_esc = escape fraction"]
        D2 --> OUT1
    end

    %% Subgraph 2: External Incidence
    subgraph E2["☀️ External Incidence"]
        style E2 fill:#0f172a,stroke:#6366f1,stroke-width:2px,color:#818cf8,font-weight:bold;
      
        A2["Photons enter from aperture"] --> B2["Planck-sampled wavelength"]
        B2 --> C2{"λ < λ_c ?"}
        C2 -->|Yes| D3["Trace into cavity<br/>TMM bounces"]
        C2 -->|No| D4["Diffract around rim"]
        D3 --> OUT2["α_eff = absorbed fraction"]
        D4 --> OUT2
    end

    %% Apply Classes to Nodes
    class C1,C2 decision;
    class D1,D3 guided;
    class D2,D4 evanescent;
    class OUT1 emissionOut;
    class OUT2 incidenceOut;

    %% Link Styling
    linkStyle 2,8 stroke:#34d399,stroke-width:2px,fill:none;
    linkStyle 3,9 stroke:#f87171,stroke-width:2px,fill:none;
```

- **P_esc** drives **epsilon_B** (emission side)
- **alpha_eff** drives **alpha_eff** (absorption side)
- Their ratio -> **decoupling ratio**

---

## 5. Thin-film and temperature-dependent optics

```mermaid
flowchart TD
    %% Class Definitions
    classDef default fill:#1e293b,stroke:#475569,stroke-width:1.5px,color:#f8fafc,font-family:sans-serif;
    classDef inputNode fill:#0f172a,stroke:#38bdf8,stroke-width:1.5px,color:#e0f2fe;
    classDef moduleHub fill:#0f172a,stroke:#38bdf8,stroke-width:2.5px,color:#38bdf8,font-weight:bold;
    classDef physicsStep fill:#3730a3,stroke:#818cf8,stroke-width:1.5px,color:#e0e7ff;
    classDef tmmStep fill:#701a75,stroke:#f0abfc,stroke-width:1.5px,color:#fdf4ff;
    classDef outputNode fill:#065f46,stroke:#10b981,stroke-width:2px,color:#ffffff,font-weight:bold;

    %% Main Flowchart Structure
    subgraph OPTICS["🧪 Material Optics Calculations"]
        style OPTICS fill:#0f172a,stroke:#0d9488,stroke-width:2px,color:#2dd4bf,font-weight:bold;
      
        L["Wavelength λ, Temperature T"] --> MO["material_optics.py"]
      
        MO --> T1["Tabulated n(λ), k(λ) at 300 K"]
        MO --> T2["Drude-Lorentz drift<br/>γ(T) damping"]
        MO --> T3["Non-local hydrodynamic<br/>correction"]
      
        T1 --> T4["TMM thin-film<br/>R, T, A"]
        T2 --> T4
        T3 --> T4
      
        T4 --> OUT["Effective n(λ,T), k(λ,T)"]
    end

    %% Apply Classes to Nodes
    class L inputNode;
    class MO moduleHub;
    class T1,T2,T3 physicsStep;
    class T4 tmmStep;
    class OUT outputNode;

    %% Corrected Link Styling (Indices 0 through 8)
    linkStyle 0 stroke:#38bdf8,stroke-width:2px,fill:none;
    linkStyle 1 stroke:#818cf8,stroke-width:2px,fill:none;
    linkStyle 2 stroke:#818cf8,stroke-width:2px,fill:none;
    linkStyle 3 stroke:#818cf8,stroke-width:2px,fill:none;
    linkStyle 4 stroke:#f0abfc,stroke-width:2px,fill:none;
    linkStyle 5 stroke:#f0abfc,stroke-width:2px,fill:none;
    linkStyle 6 stroke:#f0abfc,stroke-width:2px,fill:none;
    linkStyle 7 stroke:#34d399,stroke-width:2px,fill:none;
```

---

## 6. Waveguide modal physics

```mermaid
flowchart LR
    %% Class Definitions
    classDef default fill:#1e293b,stroke:#475569,stroke-width:1.5px,color:#f8fafc,font-family:sans-serif;
    classDef inputNode fill:#0f172a,stroke:#38bdf8,stroke-width:1.5px,color:#e0f2fe;
    classDef moduleHub fill:#0f172a,stroke:#38bdf8,stroke-width:2.5px,color:#38bdf8,font-weight:bold;
    classDef calcStep fill:#3730a3,stroke:#818cf8,stroke-width:1.5px,color:#e0e7ff;
    classDef propStep fill:#701a75,stroke:#f0abfc,stroke-width:1.5px,color:#fdf4ff;
    classDef outputNode fill:#065f46,stroke:#10b981,stroke-width:2px,color:#ffffff,font-weight:bold;

    %% Subgraph Structure
    subgraph WAVEGUIDE["🌊 Waveguide Modal Solver"]
        style WAVEGUIDE fill:#0f172a,stroke:#0d9488,stroke-width:2px,color:#2dd4bf,font-weight:bold;

        D["Cavity diameter D, Height H"] --> WG["waveguide_modes.py"]

        WG --> LC["λ_c = λ_c_PEC / n_real"]
        WG --> BETA["Complex β = β_real + iα"]
        WG --> RAP["R_ap, T_ap aperture mismatch"]
        WG --> ATTN["exp(-αL) attenuation"]

        LC --> PROP["f_prop = Planck power below λ_c"]
        BETA --> PROP

        RAP --> TT["T_total(λ) multi-mode sum"]
        ATTN --> TT
    end

    %% Apply Classes to Nodes
    class D inputNode;
    class WG moduleHub;
    class LC,BETA,RAP,ATTN calcStep;
    class PROP propStep;
    class TT outputNode;

    %% Link Styling (Explicitly mapped 0..8 to avoid render errors)
    linkStyle 0 stroke:#38bdf8,stroke-width:2px,fill:none;
    linkStyle 1,2,3,4 stroke:#818cf8,stroke-width:2px,fill:none;
    linkStyle 5,6 stroke:#f0abfc,stroke-width:2px,fill:none;
    linkStyle 7,8 stroke:#34d399,stroke-width:2px,fill:none;
```

---

## 7. Radiosity enclosure model

```mermaid
flowchart LR
    %% Class Definitions
    classDef default fill:#1e293b,stroke:#475569,stroke-width:1.5px,color:#f8fafc,font-family:sans-serif;
    classDef plateNode fill:#0f172a,stroke:#38bdf8,stroke-width:1.5px,color:#e0f2fe;
    classDef surrNode fill:#7f1d1d,stroke:#f87171,stroke-width:1.5px,color:#fef2f2;
    classDef matrixHub fill:#0f172a,stroke:#38bdf8,stroke-width:2.5px,color:#38bdf8,font-weight:bold;
    classDef outputNode fill:#065f46,stroke:#10b981,stroke-width:2px,color:#ffffff,font-weight:bold;

    %% Subgraph Structure
    subgraph RADIOSITY["🔥 4-Surface Radiosity Network"]
        style RADIOSITY fill:#0f172a,stroke:#0d9488,stroke-width:2px,color:#2dd4bf,font-weight:bold;

        %% Surface Exchange Links
        AF["Plate A front<br/>T_A, ε_A"] -->|F_AF→B| BF["Plate B front"]
        AB["Plate A back"] -->|F_AB→S| S["Surroundings<br/>T_surr"]
        BF -->|F_B→S| S
        AF -->|F_AF→S| S

        %% Matrix Calculations
        J["Radiosity matrix J"] --> Q["q_net = F × (J_AF - J_B)"]
        J --> TS["T_stag = (G_B / ε_B σ)^0.25"]
    end

    %% Apply Classes to Nodes
    class AF,BF,AB plateNode;
    class S surrNode;
    class J matrixHub;
    class Q,TS outputNode;

    %% Link Styling (Explicitly mapped 0..5 to avoid render errors)
    linkStyle 0 stroke:#38bdf8,stroke-width:2px,fill:none;
    linkStyle 1,2,3 stroke:#f87171,stroke-width:2px,fill:none;
    linkStyle 4,5 stroke:#34d399,stroke-width:2px,fill:none;
```

---

## 8. Interactive UI features

| Feature              | What it shows                                  |
| -------------------- | ---------------------------------------------- |
| Collapsible sections | Config panels collapse/expand                  |
| Quick presets        | One-click fill: 200K, 3000K, 12000K, room temp |
| Executive summary    | Plain-English interpretation                   |
| Physics status strip | Confidence dot + regime badge                  |
| Regime gauge         | Temperature scale with color zones             |
| Energy flow bars     | Proportional bars: Emitted / Net / Lost        |
| Radiation diagram    | Animated canvas: photon flow A to B            |
| Modal cutoff visual  | Planck curve with lambda markers               |
| Orchestrator banner  | Full regime breakdown + dimensionless ratios   |
| Download report      | Human-readable text with all values            |

---

## 9. One-paragraph summary

The simulator studies a micro-structured surface as a repeating cavity, asks how many photons escape from inside it and how many incoming photons get trapped, and applies thin-film, temperature-dependent complex-optical, roughness-scattering, and modal-cutoff corrections along the way. Those escape and absorption probabilities become effective emissivity and absorptivity values. The code then chooses between near-field and far-field physics and solves a four-surface enclosure to get the net radiative heat flow. The reason this matters is that a structured surface can absorb incident radiation very efficiently while emitting much less from inside its cavities - the anisotropic decoupling the model is built to capture.
