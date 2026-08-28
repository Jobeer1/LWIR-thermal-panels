/* =====================================================================
   near_field_ui.js — Track B3 UI Controls for Near-Field Physics
   
   Features:
     1. Near-field toggle switches
     2. Gap distance slider with dynamic thresholds
     3. LDOS heatmap visualization
     4. Distance decay curves (q_net vs gap)
     5. Surface-mode indicator badges
   
   Author: Kiro
   Date: August 2026
   ===================================================================== */

'use strict';

class NearFieldUIManager {
    constructor() {
        this.nearFieldEnabled = true;
        this.gapThreshold = 5.0;
        this.ldosChart = null;
        this.distanceChart = null;
        this.lastResults = null;
    }

    /**
     * Initialize near-field UI controls.
     */
    initializeControls() {
        console.log('[NearFieldUI] Initializing near-field controls...');
        
        // Create control panel if not exists
        const controlPanel = this.createNearFieldControlPanel();
        const configPanel = document.querySelector('.config-panel');
        if (configPanel && !document.getElementById('nf-control-panel')) {
            configPanel.appendChild(controlPanel);
        }

        // Bind event listeners
        this.bindEventListeners();
    }

    /**
     * Create near-field control panel HTML.
     */
    createNearFieldControlPanel() {
        const panel = document.createElement('div');
        panel.id = 'nf-control-panel';
        panel.className = 'near-field-controls';
        panel.innerHTML = `
            <div class="section-label">⚡ Near-Field Physics</div>
            
            <div class="input-row">
                <div class="input-group">
                    <label>
                        <input type="checkbox" id="enable-near-field" checked>
                        Enable Near-Field Model
                    </label>
                </div>
            </div>

            <div class="input-row">
                <div class="input-group">
                    <label for="nf-threshold">
                        Activation Threshold (gap ratio)
                    </label>
                    <input type="number" id="nf-threshold" 
                           value="5.0" step="0.5" min="1.0" max="20.0">
                </div>
                <div class="input-group">
                    <label for="nf-gap-distance">
                        Gap Distance (nm)
                    </label>
                    <input type="range" id="nf-gap-distance" 
                           value="100" min="10" max="10000" step="10">
                    <span id="nf-gap-display">100 nm</span>
                </div>
            </div>

            <div class="input-row">
                <div class="input-group">
                    <label for="nf-quadrature-omega">
                        Frequency Quadrature (n_ω)
                    </label>
                    <input type="number" id="nf-quadrature-omega" 
                           value="80" min="20" max="200" step="10">
                </div>
                <div class="input-group">
                    <label for="nf-quadrature-kp">
                        Parallel Wavevector (n_k∥)
                    </label>
                    <input type="number" id="nf-quadrature-kp" 
                           value="50" min="20" max="100" step="10">
                </div>
            </div>

            <div id="nf-regime-indicator" class="regime-badge" style="display:none;"></div>
        `;
        
        // Add styling
        const style = document.createElement('style');
        style.textContent = `
            .near-field-controls {
                margin-top: 1.5rem;
                padding: 1rem;
                border: 1px solid rgba(100, 200, 255, 0.3);
                border-radius: 8px;
                background: linear-gradient(135deg, rgba(100, 200, 255, 0.05), transparent);
            }

            .regime-badge {
                margin-top: 1rem;
                padding: 0.75rem;
                border-radius: 6px;
                font-weight: 600;
                text-align: center;
                font-size: 0.9rem;
            }

            .regime-badge.near-field {
                background: #fff3cd;
                border: 1px solid #ffc107;
                color: #856404;
            }

            .regime-badge.far-field {
                background: #d1ecf1;
                border: 1px solid #17a2b8;
                color: #0c5460;
            }

            #nf-gap-display {
                display: inline-block;
                min-width: 80px;
                text-align: right;
                font-family: 'Space Mono', monospace;
                font-size: 0.9rem;
                color: #666;
            }
        `;
        document.head.appendChild(style);
        
        return panel;
    }

    /**
     * Bind event listeners to UI controls.
     */
    bindEventListeners() {
        const gapSlider = document.getElementById('nf-gap-distance');
        const gapDisplay = document.getElementById('nf-gap-display');
        
        if (gapSlider && gapDisplay) {
            gapSlider.addEventListener('input', (e) => {
                const nm = parseInt(e.target.value);
                gapDisplay.textContent = `${nm} nm`;
                this.updateRegimeIndicator(nm);
            });
        }

        const enableNF = document.getElementById('enable-near-field');
        if (enableNF) {
            enableNF.addEventListener('change', (e) => {
                this.nearFieldEnabled = e.target.checked;
                console.log(`[NearFieldUI] Near-field: ${this.nearFieldEnabled}`);
            });
        }
    }

    /**
     * Update regime indicator (near-field vs far-field).
     */
    updateRegimeIndicator(gapNm) {
        const T_peak = 600.0;  // Assume hot plate temp
        const lambda_peak_um = 2898.0 / T_peak;  // Wien's displacement
        const lambda_peak_m = lambda_peak_um * 1e-6;
        
        const gap_m = gapNm * 1e-9;
        const threshold_m = lambda_peak_m / (2.0 * Math.PI);
        const gap_ratio = gap_m / threshold_m;
        
        const indicator = document.getElementById('nf-regime-indicator');
        if (!indicator) return;
        
        indicator.style.display = 'block';
        
        if (gap_ratio < 5.0) {
            indicator.className = 'regime-badge near-field';
            indicator.innerHTML = `⚡ <strong>NEAR-FIELD MODE</strong><br/>
                                   Gap ratio: ${gap_ratio.toFixed(2)} (< 5.0)`;
        } else {
            indicator.className = 'regime-badge far-field';
            indicator.innerHTML = `📡 <strong>FAR-FIELD MODE</strong><br/>
                                   Gap ratio: ${gap_ratio.toFixed(2)} (> 5.0)`;
        }
    }

    /**
     * Collect near-field parameters from form.
     */
    getNearFieldParams() {
        return {
            enable_near_field: document.getElementById('enable-near-field')?.checked ?? true,
            near_field_threshold: parseFloat(
                document.getElementById('nf-threshold')?.value ?? '5.0'
            ),
            near_field_n_omega: parseInt(
                document.getElementById('nf-quadrature-omega')?.value ?? '80'
            ),
            near_field_n_kparallel: parseInt(
                document.getElementById('nf-quadrature-kp')?.value ?? '50'
            ),
            gap: parseInt(
                document.getElementById('nf-gap-distance')?.value ?? '100'
            ) * 1e-9,  // Convert nm to m
        };
    }

    /**
     * Display LDOS heatmap from simulation results.
     */
    displayLDOSHeatmap(results) {
        if (!results.ldos_peak_ratio) {
            console.log('[NearFieldUI] No LDOS data in results');
            return;
        }

        const container = document.getElementById('ldos-visualization');
        if (!container) {
            console.log('[NearFieldUI] No LDOS visualization container');
            return;
        }

        const ldosRatio = results.ldos_peak_ratio || 1.0;
        const dominantWavelength = results.dominant_wavelength_um || 10.0;

        container.innerHTML = `
            <div class="ldos-card">
                <h4>🌊 LDOS Enhancement</h4>
                <div class="ldos-metric">
                    <span class="metric-label">Peak Ratio:</span>
                    <span class="metric-value">${ldosRatio.toFixed(2)}×</span>
                </div>
                <div class="ldos-metric">
                    <span class="metric-label">Dominant λ:</span>
                    <span class="metric-value">${dominantWavelength.toFixed(1)} µm</span>
                </div>
                <div class="ldos-bar">
                    <div class="ldos-fill" style="width: ${Math.min(ldosRatio * 10, 100)}%;"></div>
                </div>
            </div>
        `;

        // Add styling
        if (!document.getElementById('ldos-styles')) {
            const style = document.createElement('style');
            style.id = 'ldos-styles';
            style.textContent = `
                .ldos-card {
                    background: linear-gradient(135deg, rgba(150, 100, 255, 0.1), transparent);
                    border: 1px solid rgba(150, 100, 255, 0.3);
                    border-radius: 8px;
                    padding: 1rem;
                    margin: 1rem 0;
                }

                .ldos-metric {
                    display: flex;
                    justify-content: space-between;
                    align-items: center;
                    padding: 0.5rem 0;
                    font-size: 0.95rem;
                }

                .metric-label {
                    font-weight: 600;
                    color: #444;
                }

                .metric-value {
                    font-family: 'Space Mono', monospace;
                    font-weight: 700;
                    color: #9964ff;
                    font-size: 1.1rem;
                }

                .ldos-bar {
                    width: 100%;
                    height: 20px;
                    background: #f0f0f0;
                    border-radius: 4px;
                    overflow: hidden;
                    margin-top: 0.5rem;
                }

                .ldos-fill {
                    height: 100%;
                    background: linear-gradient(90deg, #9964ff, #ff64b0);
                    transition: width 0.3s ease;
                }
            `;
            document.head.appendChild(style);
        }

        this.lastResults = results;
    }

    /**
     * Display distance-decay curve (q_net vs gap).
     */
    displayDistanceDecay(gapData, fluxData) {
        const container = document.getElementById('distance-decay-chart');
        if (!container) return;

        // Simple ASCII chart fallback (no Chart.js dependency)
        const html = this.renderDistanceDecayChart(gapData, fluxData);
        container.innerHTML = html;
    }

    /**
     * Render simple ASCII distance-decay chart.
     */
    renderDistanceDecayChart(gaps, fluxes) {
        if (!gaps || gaps.length === 0) return '';

        const maxFlux = Math.max(...fluxes);
        const minFlux = Math.min(...fluxes);
        const range = maxFlux - minFlux || 1.0;

        let html = '<div class="decay-chart"><pre>';
        
        // Header
        html += `Gap (nm)  │ Flux (W/m²)\n`;
        html += '─────────┼─────────────────────────────────\n';

        // Data points
        for (let i = 0; i < gaps.length; i++) {
            const gapNm = gaps[i] * 1e9;
            const flux = fluxes[i];
            const normalized = (flux - minFlux) / range;
            const barLength = Math.round(normalized * 30);
            const bar = '█'.repeat(barLength);
            
            html += `${gapNm.toFixed(0).padStart(7)} │ ${bar}\n`;
        }

        html += '</pre></div>';
        return html;
    }

    /**
     * Update UI with simulation results.
     */
    updateResults(results) {
        console.log('[NearFieldUI] Updating results...', results);

        // Update LDOS visualization
        if (results.ldos_peak_ratio !== undefined) {
            this.displayLDOSHeatmap(results);
        }

        // Update physics regime badge
        if (results.physics_regime) {
            const badge = document.getElementById('regime-badge');
            if (badge) {
                badge.textContent = results.physics_regime.toUpperCase();
                badge.className = `badge ${results.physics_regime}`;
            }
        }

        // Update evanescent fraction
        if (results.evanescent_fraction !== undefined) {
            const evanEl = document.getElementById('evanescent-fraction');
            if (evanEl) {
                evanEl.textContent = `${(results.evanescent_fraction * 100).toFixed(1)}%`;
            }
        }

        this.lastResults = results;
    }
}

// Initialize on DOM ready
document.addEventListener('DOMContentLoaded', () => {
    window.nearFieldUI = new NearFieldUIManager();
    window.nearFieldUI.initializeControls();
    
    console.log('[NearFieldUI] Initialized');
});
