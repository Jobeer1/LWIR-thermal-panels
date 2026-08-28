/* =====================================================================
   app.js — Monte Carlo Radiative Exchange Simulator
   Handles: geometry mode toggle, form submission, result display,
            CI labels, near-field warnings, spectral cards.
   ===================================================================== */

'use strict';

// ---- Geometry mode toggle -----------------------------------------------

const modeBtns = document.querySelectorAll('.mode-btn');
const modeInput = document.getElementById('geometry_mode');
const fieldsHoneycomb = document.getElementById('fields-honeycomb');
const fieldsForest  = document.getElementById('fields-forest');
let lastSimulation = null;

const downloadReportBtn = document.getElementById('download-report-btn');
if (downloadReportBtn) {
    downloadReportBtn.addEventListener('click', () => {
        if (!lastSimulation) return;

        const report = [
            'Monte Carlo Radiative Exchange Simulation Report',
            `Generated: ${new Date().toLocaleString()}`,
            '',
            'INPUT VARIABLES',
            '==============',
            JSON.stringify(lastSimulation.inputs, null, 2),
            '',
            'OUTPUT VARIABLES',
            '================',
            JSON.stringify(lastSimulation.outputs, null, 2),
            '',
        ].join('\n');

        const blob = new Blob([report], { type: 'text/plain;charset=utf-8' });
        const url = URL.createObjectURL(blob);
        const link = document.createElement('a');
        link.href = url;
        link.download = `radiative-simulation-${new Date().toISOString().slice(0, 10)}.txt`;
        document.body.appendChild(link);
        link.click();
        link.remove();
        URL.revokeObjectURL(url);
    });
}

function switchMode(mode) {
    modeInput.value = mode;
    modeBtns.forEach(b => b.classList.toggle('active', b.dataset.mode === mode));
    
    if (fieldsHoneycomb) fieldsHoneycomb.style.display = mode === 'honeycomb' ? '' : 'none';
    if (fieldsForest) fieldsForest.style.display = mode === 'cnt_forest' ? '' : 'none';
}

modeBtns.forEach(btn => {
    btn.addEventListener('click', () => switchMode(btn.dataset.mode));
});

// ---- Full-gap MC toggle -------------------------------------------------

const fullGapCb  = document.getElementById('full_gap_mc');
const gapMcRow   = document.getElementById('gap-mc-row');
fullGapCb.addEventListener('change', () => {
    gapMcRow.style.display = fullGapCb.checked ? '' : 'none';
});


// ---- Collapsible handlers -------------------------------------------------
document.querySelectorAll(".collapsible-header").forEach(function(hdr){
    hdr.addEventListener("click",function(){var p=hdr.closest(".collapsible");if(p)p.classList.toggle("open");});
});

// ---- Quick presets ---------------------------------------------------------
var PRESETS={
    subwave:{temp_a:200,temp_b:200,temp_surr:3,emissivity_a:0.981,emissivity_a_back:0.081,height:200.001,gap:200.001,cavity_diameter:4,wall_thickness:0.051,emissivity_cnt_hc:0.8,emissivity_base_hc:0.05,eps_flat_wall:0.8,n_photons:2000},
    propagating:{temp_a:3000,temp_b:3000,temp_surr:3,emissivity_a:0.981,emissivity_a_back:0.081,height:200.001,gap:200.001,cavity_diameter:4,wall_thickness:0.051,emissivity_cnt_hc:0.8,emissivity_base_hc:0.05,eps_flat_wall:0.8,n_photons:2000},
    shortwave:{temp_a:12000,temp_b:12000,temp_surr:3,emissivity_a:0.981,emissivity_a_back:0.081,height:200.001,gap:200.001,cavity_diameter:4,wall_thickness:0.051,emissivity_cnt_hc:0.8,emissivity_base_hc:0.05,eps_flat_wall:0.8,n_photons:2000},
    room:{temp_a:300,temp_b:300,temp_surr:300,emissivity_a:0.981,emissivity_a_back:0.081,height:200.001,gap:200.001,cavity_diameter:4,wall_thickness:0.051,emissivity_cnt_hc:0.8,emissivity_base_hc:0.05,eps_flat_wall:0.8,n_photons:2000},
};
document.querySelectorAll("[data-preset]").forEach(function(btn){
    btn.addEventListener("click",function(){
        var key=btn.dataset.preset;
        window.switchMode("honeycomb");
    if(!PRESETS[key])return;
        for(var k in PRESETS[key]){var id=document.getElementById(k);if(id)id.value=PRESETS[key][k];}
        document.querySelectorAll(".collapsible").forEach(function(c){c.classList.add("open");});
        setTimeout(function(){document.getElementById("sim-form").requestSubmit();},100);
    });
});

// ---- Form submit ---------------------------------------------------------

document.getElementById('sim-form').addEventListener('submit', async (e) => {
    e.preventDefault();

    const runBtn  = document.getElementById('run-btn');
    const btnText = runBtn.querySelector('.btn-text');
    const loader  = runBtn.querySelector('.loader');

        const resultIds = [
        'res-leak','res-epsilon-b','res-emitted-flux','res-enhancement','res-esc','res-abs',
        'res-kirchhoff','res-view-factor','res-cutoff','res-fprop','res-confine',
        'res-decoup','res-omega','res-total-heat','res-coupling',
        'res-aback','res-flux-afront','res-eps-a-spec','res-spec-corr',
    ];
    const loaderHtml = '<div class="loader" style="width:22px;height:22px;margin:0 auto;"></div>';

    // Loading state
    runBtn.disabled = true;
    btnText.classList.add('hidden');
    loader.classList.remove('hidden');
    resultIds.forEach(id => {
        const el = document.getElementById(id);
        if (el) el.innerHTML = loaderHtml;
    });

        // Clear CI labels and banners
    ['ci-pesc','ci-alpha','ci-epsb'].forEach(id => {
        const el = document.getElementById(id);
        if (el) el.textContent = '';
    });
    _hideBanner('near-field-banner');

    // Gather form data (all named inputs)
    const formData = new FormData(e.target);
    const data = Object.fromEntries(formData.entries());

    // Checkbox: ensure it sends its value if checked
    if (!data.full_gap_mc) data.full_gap_mc = 'false';

    try {
        const response = await fetch('/api/simulate', {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify(data),
        });

        const result = await response.json();

        if (result.status === 'success') {
            lastSimulation = { inputs: data, outputs: result.results };
            if (downloadReportBtn) downloadReportBtn.disabled = false;
            _renderResults(result.results);
        } else {
            alert('Simulation error:\n' + result.message +
                  (result.traceback ? '\n\n' + result.traceback : ''));
            resultIds.forEach(id => {
                const el = document.getElementById(id);
                if (el) el.innerHTML = '—';
            });
        }
    } catch (err) {
        console.error(err);
        alert('Failed to reach the Flask server. Is it running on port 5000?');
    } finally {
        runBtn.disabled = false;
        btnText.classList.remove('hidden');
        loader.classList.add('hidden');
    }
});

// ---- Result rendering ---------------------------------------------------

function _renderResults(r) {
    // Executive summary
    _renderExecutiveSummary(r);

    // Solver-mode badge (Phase 6): ray fallback vs cached full-wave
    _updateSolverBadge(r.solver_mode, r.wave_model, r.wave_response_info);
    
    // Physics regime badge (Phase 3): near-field vs far-field
    if (r.physics_regime !== undefined && r.gap_ratio !== undefined) {
        _updatePhysicsRegimeBadge(r.physics_regime, r.gap_ratio);
    }

    // MC results
    const pEsc     = (r.p_esc != null) ? (r.p_esc     * 100).toFixed(3) : '—';
    const aEff     = (r.alpha_eff * 100).toFixed(3);
    const epsilonB = (r.epsilon_b * 100).toFixed(3);
    const kErr     = (Number(r.kirchhoff_error) * 100).toFixed(2);
    const enh      = Number(r.cavity_enhancement).toFixed(2);
    const vf       = r.view_factor_A_B != null
        ? Number(r.view_factor_A_B).toFixed(5)
        : 'n/a';

    // Flux values
    const flux    = _fmtW(r.net_flux_A_front);
    const qNet    = _fmtW(r.total_leakage,   6);
    const qAtoB   = _fmtW(r.q_af_to_b,      6);
    const qAback  = _fmtW(r.q_ab_net,        6);
    const fluxB   = _fmtW(r.flux_emitted_b);

    // One-way physical emission fluxes
    const qEmitAtoB   = (r.q_emit_a_to_b_one_way != null) ? _fmtW(r.q_emit_a_to_b_one_way, 1) : null;
    const qEmitBtoA   = (r.q_emit_b_to_a_one_way != null) ? _fmtW(r.q_emit_b_to_a_one_way, 1) : null;
    const qNetPhys    = (r.q_net_a_to_b_physical != null) ? _fmtW(r.q_net_a_to_b_physical, 1) : null;

        // Animate numeric results
    _animateValue('res-esc',       0, pEsc,     900, '%');
    _animateValue('res-abs',       0, aEff,     900, '%');
    _animateValue('res-epsilon-b', 0, epsilonB, 900, '%');
    _animateValue('res-kirchhoff', 0, kErr,     900, '%');
    _setText('res-emitted-flux', fluxB);
    _setText('res-leak',       flux);
    _setText('res-enhancement', enh + '×');
    _setText("res-view-factor", vf);
    _setText("res-total-heat", qNet);
    _setText("res-coupling", _fmtW(r.q_af_to_b, 6));
    _setText("res-aback", _fmtW(r.q_ab_net, 6));
    // Physical net flux A->B (accounts for eps_A != eps_B)
    _setText("res-flux-afront", _fmtW(r.q_net_a_to_b_physical != null ? r.q_net_a_to_b_physical : r.net_flux_A_front, 1));
    // Stagnation temperature (adiabatic wall temperature of Plate B)
    if (r.T_B_stag != null) {
        _setText("res-t-stag", (Number(r.T_B_stag) - 273.15).toFixed(1) + String.fromCharCode(176) + "C  (" + Number(r.T_B_stag).toFixed(1) + " K)");
    }
    // Bug #13 - wave-optics decoupling diagnostics
    _setText('res-cutoff',
        r.cutoff_wavelength_um != null
            ? Number(r.cutoff_wavelength_um).toFixed(4) + ' µm' : '∞ (no cutoff)');
    if (r.propagating_fraction != null)
        _animateValue('res-fprop', 0, Number(r.propagating_fraction) * 100, 700, '%');
    if (r.confinement_pct != null)
        _animateValue('res-confine', 0, Number(r.confinement_pct), 700, '%');
    if (r.decoupling_ratio != null)
        _setText('res-decoup', Number(r.decoupling_ratio).toFixed(3) + '×');
    if (r.escape_solid_angle_sr != null)
        _setText('res-omega', Number(r.escape_solid_angle_sr).toExponential(2) + ' sr');

    // Phase 3 near-field diagnostics
    if (r.physics_regime === 'near-field') {
        if (r.evanescent_fraction !== undefined && r.evanescent_fraction > 0) {
            const evan_pct = (100 * Number(r.evanescent_fraction)).toFixed(1);
            const evan_flux = _fmtW(r.evanescent_flux_W_m2, 2);
            const prop_flux = _fmtW(r.propagating_flux_W_m2, 2);
            _setText('near-field-info',
                `<strong>Evanescent waves:</strong> ${evan_pct}% of total flux (${evan_flux} W/m²)<br/>` +
                `<strong>Propagating:</strong> ${prop_flux} W/m²`);
        }
    }

    // 95% Confidence intervals
    const isCached = r.solver_mode === 'cached';
    if (r.p_esc_ci95 != null) {
        _setText('ci-pesc', isCached
            ? 'deterministic (cached table)'
            : `±${(r.p_esc_ci95 * 100).toFixed(4)}% (95% CI)`);
    }
    if (r.alpha_eff_ci95 != null) {
        _setText('ci-alpha', isCached
            ? 'deterministic (cached table)'
            : `±${(r.alpha_eff_ci95 * 100).toFixed(4)}% (95% CI)`);
    }
    if (r.epsilon_b_ci95 != null) {
        _setText('ci-epsb', isCached || Number(r.epsilon_b_ci95) <= 0
            ? 'deterministic macro operator'
            : `±${(r.epsilon_b_ci95 * 100).toFixed(4)}% (95% CI)`);
    }
        // ε_B < α_eff (anisotropic decoupling — the correct physics)
    if (r.epsilon_b_cavity_part != null || r.epsilon_b_flat_part != null) {
        const cav  = r.epsilon_b_cavity_part != null ? (r.epsilon_b_cavity_part * 100).toFixed(2) + '%' : '—';
        const flat = r.epsilon_b_flat_part   != null ? (r.epsilon_b_flat_part   * 100).toFixed(2) + '%' : '—';
        _setText('res-epsb-breakdown', `cavity ${cav} + flat top ${flat}`);
    }

    // Spectral
    if (r.eps_a_spectral != null) {
        _setText('res-eps-a-spec', (r.eps_a_spectral * 100).toFixed(2) + '%');
    }
    if (r.spectral_correction_pct != null) {
        _setText('res-spec-corr', r.spectral_correction_pct.toFixed(2) + '%');
    }

        // Geometry label
    _setText('geometry-label', r.geometry_label || '—');

    // Wall thickness — dynamic footnote value (replaces hardcoded "0.3 µm")
    _setText('wall-thickness-val', r.wall_thickness_um != null
        ? r.wall_thickness_um.toFixed(3) : '—');

    // Kirchhoff card colour: reflects the anisotropic decoupling ratio
    const kCard = document.getElementById('kirchhoff-card');
    if (kCard) {
        kCard.classList.remove('kirchhoff-good', 'kirchhoff-bad', 'kirchhoff-decoupled');
        if (r.decoupling_ratio != null) {
            const ratio = r.decoupling_ratio;
            if (ratio > 0.95) {
                // Above-cutoff gray cavity: reciprocity holds
                kCard.classList.add('kirchhoff-good');
            } else {
                // Sub-wavelength: anisotropic decoupling active (α >> ε)
                kCard.classList.add('kirchhoff-decoupled');
            }
        } else {
            kCard.classList.add('kirchhoff-good');
        }
    }

    // Near-field warning
    if (r.near_field_warning) {
        _showBanner('near-field-banner', r.near_field_warning);
    }
}

// ---- Solver-mode badge (Phase 6 UI indicator) ------------------------------

function _updateSolverBadge(solverMode, waveModel, waveInfo) {
    const badge = document.getElementById('solver-badge');
    if (!badge) return;

    const mode = (solverMode || waveModel || 'ray');
    if (mode === 'cached') {
        badge.className = 'solver-badge solver-cached';
        let label = 'Solver Mode: Cached Full-Wave';
        if (waveInfo && waveInfo.energy_conservation_ok) {
            label += ' ✓';
        }
        badge.title = waveInfo
            ? `Cache: ${waveInfo.source} · grid ${JSON.stringify(waveInfo.response_shape)} · ` +
              `max δE ${Number(waveInfo.energy_balance_error).toExponential(1)} · ` +
              `λ ${Number(waveInfo.wavelength_range_um[0]).toFixed(2)}–${Number(waveInfo.wavelength_range_um[1]).toFixed(1)} µm`
            : 'Pre-computed full-wave response table';
        badge.textContent = label;
    } else {
        badge.className = 'solver-badge solver-ray';
        badge.textContent = 'Solver Mode: Ray Tracing (Fallback)';
        badge.title = 'Monte Carlo 3-D ray tracer';
    }
}

// ---- Physics regime badge (Phase 3 near-field indicator) --------------------

function _updatePhysicsRegimeBadge(physicsRegime, gapRatio) {
    const badge = document.getElementById('physics-regime-badge');
    if (!badge) return;
    
    if (physicsRegime === 'near-field') {
        badge.className = 'physics-regime-badge physics-regime-near-field';
        badge.innerHTML = `⚡ NEAR-FIELD MODE (Gap Ratio: ${Number(gapRatio).toFixed(2)})`;
        badge.title = 'Polder-Van Hove near-field radiative transfer active (Phase 3)';
    } else if (physicsRegime && physicsRegime.includes('far-field')) {
        badge.className = 'physics-regime-badge physics-regime-far-field';
        badge.innerHTML = `📡 FAR-FIELD MODE (Gap Ratio: ${Number(gapRatio).toFixed(2)})`;
        badge.title = 'Classical far-field radiosity calculation (Phase 0/6)';
    } else {
        badge.className = 'physics-regime-badge physics-regime-unknown';
        badge.innerHTML = '? Unknown Regime';
        badge.title = 'Physics regime not determined';
    }
}

// ---- Helpers -------------------------------------------------------------

function _setText(id, html) {
    const el = document.getElementById(id);
    if (el) el.innerHTML = html;
}

function _fmtW(val, dec = 1) {
    if (val == null) return '—';
    return Number(val).toLocaleString('en-US', {
        maximumFractionDigits: dec,
        minimumFractionDigits: dec,
    });
}

function _showBanner(id, msg) {
    const el = document.getElementById(id);
    if (!el) return;
    el.innerHTML = msg;
    el.classList.remove('hidden');
}

function _hideBanner(id) {
    const el = document.getElementById(id);
    if (el) { el.classList.add('hidden'); el.innerHTML = ''; }
}

function _animateValue(id, start, end, duration, suffix) {
    const el = document.getElementById(id);
    if (!el) return;
    const endNum = parseFloat(end);
    let startTs = null;
    // Executive summary (plain English)
    _renderExecutiveSummary(r);


    const step = (ts) => {
        if (!startTs) startTs = ts;
        const progress = Math.min((ts - startTs) / duration, 1);
        // Ease-out quart
        const eased = 1 - Math.pow(1 - progress, 4);
        el.innerHTML = (eased * (endNum - start) + start).toFixed(2) + suffix;
        if (progress < 1) requestAnimationFrame(step);
        else el.innerHTML = endNum.toFixed(2) + suffix;
    };
    requestAnimationFrame(step);
}


// ---- Physics status strip -----------------------------------------------
function _renderStatusStrip(r) {
    var el=document.getElementById("physics-status-strip");
    if(!el)return; if(!r){el.innerHTML="";return;}
    var conf=Number(r.physics_confidence||100);
    var dc="safe",st="Normal confidence";
    if(conf<60){dc="critical";st="Low confidence";}
    else if(conf<85){dc="warning";st="Moderate confidence";}
    var rg=(r.physics_orchestrator||{}).geometry_regime||"?";
    el.innerHTML="<span class=\"status-dot "+dc+"\"></span> <span>"+st+"</span> <span class=\"status-regime\">"+rg+"</span>";
}

// ---- Regime gauge -------------------------------------------------------
function _renderRegimeGauge(orch,T){
    var el=document.getElementById("regime-gauge");
    if(!el||!orch||!T){if(el)el.innerHTML="";return;}
    var pct=Math.min(100,Math.max(0,(Number(T)/16000)*100));
    el.innerHTML="<div class=\"gauge-track\"><div class=\"gauge-needle\" style=\"left:"+pct+"%;\"></div></div><div class=\"gauge-labels\"><span class=\"gauge-label\">RAY</span><span class=\"gauge-label\">FULLWAVE</span><span class=\"gauge-label\">EMT</span></div><div style=\"text-align:center;font-size:0.72rem;color:var(--text-muted);\">T="+Number(T).toLocaleString()+"K "+(orch.geometry_regime||"?")+"</div>";
}

// ---- Energy flow bars ---------------------------------------------------
function _renderEnergyFlow(r){
    var el=document.getElementById("energy-flow");
    if(!el||!r||!r.flux_emitted_b){if(el)el.innerHTML="";return;}
    var qE=Number(r.flux_emitted_b||0),qN=Number(r.net_flux_A_front||0),qB=Number(r.q_ab_net||r.net_back_loss||0),mv=Math.max(qE,qN,qB,1);
    var rows="";
    var arr=[["Emitted",qE,"emiss"],["Net A->B",qN,"net"],["Lost back",qB,"loss"]];
    for(var i=0;i<3;i++){var d=arr[i];rows=rows+"<div class=\"flow-row\"><span class=\"flow-label\">"+d[0]+"</span><div class=\"flow-bar-track\"><div class=\"flow-bar-fill "+d[2]+"\" style=\"width:"+Math.max(2,d[1]/mv*100)+"%;\"></div></div><span class=\"flow-value\">"+_fmtW(d[1],0)+"</span></div>";}
    el.innerHTML=rows+"<div style=\"font-size:0.68rem;color:var(--text-muted);text-align:right;\">W/m2</div>";
}

// ---- Executive summary --------------------------------------------------
function _renderExecutiveSummary(r){
    var el=document.getElementById("exec-summary");
    if(!el)return;
    if(!r||!r.epsilon_b){el.classList.remove("visible");el.innerHTML="";return;}
    var eB=(r.epsilon_b*100).toFixed(1),aE=(r.alpha_eff*100).toFixed(1),dR=Number(r.decoupling_ratio||0),qN=Number(r.net_flux_A_front||0),tS=Number(r.T_B_stag||0);
    var l=[];
    if(dR>3)l.push("Cavity traps light: "+aE+"% absorbed, only <b>"+eB+"%</b> re-emitted &mdash; <b>"+dR.toFixed(1)+"x</b> decoupling.");
    else l.push("Near-blackbody: epsB="+eB+"%, alpha_eff="+aE+"%.");
    if(qN>1){var z=qN>=1e6?(qN/1e6).toFixed(1)+" MW/m2":qN>=1e3?(qN/1e3).toFixed(1)+" kW/m2":qN.toFixed(1)+" W/m2";l.push("Net flux <b>"+z+"</b>.");}
    if(tS>0)l.push("Stagnation: <b>"+(tS-273.15).toFixed(1)+"C ("+tS.toFixed(1)+" K)</b>.");
    el.innerHTML=l.join(" ");el.classList.add("visible");
}


// ---- Modal cutoff visual ----------------------------------------------------
function _renderModalVisual(r){
    var c=document.getElementById("modal-canvas");if(!c)return;
    var ctx=c.getContext("2d"),dpr=window.devicePixelRatio||1;
    c.width=720*dpr;c.height=120*dpr;c.style.width="720px";c.style.height="120px";
    ctx.scale(dpr,dpr);ctx.clearRect(0,0,720,120);
    var lC=Number(r.cutoff_wavelength_um||0),tB=Number(r.temp_b||r.temperature_B||300);
    var lP=r.operator_wavelength_um||(2898/Math.max(tB,1));
    var fP=Number(r.propagating_fraction||0),fE=Math.max(0,1-fP);
    if(!lC||!lP){ctx.fillStyle="#6b8174";ctx.font="12px DM Sans";ctx.textAlign="center";ctx.fillText("No data",360,62);return;}
    var x0=30,xM=690,yB=102,yT=10;
    var ls=function(lam){return x0+660*Math.log(lam/0.1)/Math.log(1000);};
    var cX=ls(lC);
    ctx.fillStyle="rgba(0,122,77,0.10)";ctx.fillRect(x0,yT,cX-x0,yB-yT);
    ctx.fillStyle="rgba(200,47,47,0.08)";ctx.fillRect(cX,yT,xM-cX,yB-yT);
    ctx.beginPath();
    for(var i=0;i<=200;i++){var lam=0.1*Math.exp(i/200*Math.log(1000));var x=lam/100;var s=Math.pow(x,-5)/(Math.exp(1/Math.max(x*120,0.001))-1);var px=ls(lam),py=yB-s*3400;i===0?ctx.moveTo(px,Math.min(py,yB)):ctx.lineTo(px,Math.min(py,yB));}
    ctx.strokeStyle="#17231d";ctx.lineWidth=2;ctx.stroke();
    function dl(pos,color){ctx.beginPath();ctx.moveTo(pos,yT);ctx.lineTo(pos,yB);ctx.strokeStyle=color;ctx.lineWidth=1.5;ctx.setLineDash([3,3]);ctx.stroke();ctx.setLineDash([]);}
    dl(ls(lP),"#c65b1a");dl(cX,"#002395");
    ctx.fillStyle="#c65b1a";ctx.font="bold 9px DM Sans";ctx.textAlign="center";ctx.fillText("peak="+Number(lP).toFixed(2)+"um",ls(lP),yT-4);
    ctx.fillStyle="#002395";ctx.fillText("cutoff="+lC.toFixed(2)+"um",cX,yB+14);
    ctx.fillStyle="#007a4d";ctx.font="9px DM Sans";ctx.textAlign="left";ctx.fillText("Prop:"+(fP*100).toFixed(3)+"%",x0+4,yT+12);
    ctx.fillStyle="#c82f2f";ctx.fillText("Evan:"+(fE*100).toFixed(3)+"%",x0+4,yT+26);
    if(cX-x0>40){ctx.fillStyle="#007a4d";ctx.fillRect(x0,yT+34,cX-x0,6);ctx.fillStyle="#c82f2f";ctx.fillRect(cX,yT+34,xM-cX,6);}
}

// ---- Radiation diagram (animated) -------------------------------------------
var _radAnimRunning=false;
function _renderRadiationDiagram(r){
    var c=document.getElementById("rad-canvas");if(!c)return;
    var ctx=c.getContext("2d"),dpr=window.devicePixelRatio||1;
    if(!c._init){c._W=720;c._H=240;c.width=c._W*dpr;c.height=c._H*dpr;c.style.width=c._W+"px";c.style.height=c._H+"px";c._init=true;}
    ctx.setTransform(dpr,0,0,dpr,0,0);ctx.clearRect(0,0,720,240);
    ctx.fillStyle="#c65b1a";ctx.fillRect(56,60,8,120);
    ctx.fillStyle="#002395";ctx.fillRect(656,60,8,120);
    ctx.fillStyle="#17231d";ctx.font="bold 11px DM Sans";ctx.textAlign="center";
    ctx.fillText("A",60,175);ctx.fillText("B",660,175);
    var now=Date.now()/1000,nP=Math.min(20,Math.max(5,600/24));
    for(var i=0;i<nP;i++){var t=((now*0.3+i/nP)%1);
        ctx.beginPath();ctx.arc(60+t*600,120+Math.sin(t*Math.PI*3+i)*18,2+3*(1-Math.abs(t-0.5)*2),0,7);
        ctx.fillStyle="rgba(197,87,35,"+((0.4+0.6*(1-Math.abs(t-0.5)*2))*0.7)+")";ctx.fill();}
    ctx.fillStyle="#6b8174";ctx.font="9px DM Sans";ctx.textAlign="center";
    var q=r.net_flux_A_front;if(q){var qs=q>=1e6?(q/1e6).toFixed(1)+" MW/m2":q>=1e3?(q/1e3).toFixed(1)+" kW/m2":q.toFixed(1)+" W/m2";ctx.fillText("q="+qs,360,18);}
    if(!_radAnimRunning){_radAnimRunning=true;(function loop(){if(!_radAnimRunning)return;window._radAnimId=requestAnimationFrame(loop);_renderRadiationDiagram(r);})();}
}
