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
    // MC results
    const pEsc     = (r.p_esc     * 100).toFixed(3);
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

    // 95% Confidence intervals
    if (r.p_esc_ci95 != null) {
        _setText('ci-pesc',  `±${(r.p_esc_ci95 * 100).toFixed(4)}% (95% CI)`);
    }
    if (r.alpha_eff_ci95 != null) {
        _setText('ci-alpha', `±${(r.alpha_eff_ci95 * 100).toFixed(4)}% (95% CI)`);
    }
    if (r.epsilon_b_ci95 != null) {
        _setText('ci-epsb',  `±${(r.epsilon_b_ci95 * 100).toFixed(4)}% (95% CI)`);
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
        ? r.wall_thickness_um.toFixed(2) : '—');

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
