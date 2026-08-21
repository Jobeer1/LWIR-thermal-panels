import math

from flask import Flask, render_template, request, jsonify
from simulator import run_simulation

app = Flask(__name__)


@app.route('/')
def index():
    return render_template('index.html')


@app.route('/api/simulate', methods=['POST'])
def simulate():
    data = request.json
    try:
        # ---- Geometry mode -----------------------------------------------
        geometry_mode = str(data.get('geometry_mode', 'honeycomb'))

        # ---- Shared Plate B cavity parameters ----------------------------
        height     = float(data.get('height', 450.0))        # µm
        eps_flat_wall = float(data.get('eps_flat_wall', 0.1))

        if geometry_mode == 'honeycomb':
            alpha_cnt  = float(data.get('emissivity_cnt_hc', 0.95))
            alpha_ag   = float(data.get('emissivity_base_hc', 0.95))
        else:
            alpha_cnt  = float(data.get('emissivity_cnt_forest', 0.98))
            alpha_ag   = float(data.get('emissivity_base_forest', 0.40))

        # ---- Honeycomb-specific ------------------------------------------
        cavity_diameter = float(data.get('cavity_diameter', 20.0))  # µm
        wall_thickness  = float(data.get('wall_thickness', 1.0))    # µm
        if not math.isfinite(cavity_diameter) or cavity_diameter < 0.001:
            raise ValueError('Honeycomb cavity diameter must be at least 0.001 µm (1 nm).')
        if not math.isfinite(wall_thickness) or wall_thickness < 0.051:
            raise ValueError('Honeycomb wall thickness must be at least 0.051 µm.')

        # Calculate packing fraction for a hexagonal close-packed array of holes
        pitch_hc = cavity_diameter + wall_thickness
        if pitch_hc > 0:
            packing_fraction = (math.pi / (2 * math.sqrt(3))) * (cavity_diameter / pitch_hc)**2
        else:
            packing_fraction = 0.9069

        # ---- CNT Forest-specific (pitch in µm, diameters in nm) ----------
        cnt_pitch    = float(data.get('cnt_pitch', 0.05))     # µm
        cnt_dia_base = float(data.get('cnt_diameter_base', 10.0))  # nm
        cnt_dia_top  = float(data.get('cnt_diameter_top',  5.0))   # nm

        # ---- Legacy rect pit params (kept for compatibility) --------------
        width = float(data.get('pitch', 10.0))
        depth = float(data.get('trench_depth', 10.0))

        # ---- Plate A parameters ------------------------------------------
        temp_a          = float(data.get('temp_a', 600.0))
        emissivity_a    = float(data.get('emissivity_a', 1.0))
        emissivity_a_back = float(data.get('emissivity_a_back', 0.1))
        width_a         = float(data.get('width_a', 1000.0))  # µm
        depth_a         = float(data.get('depth_a', 1000.0))  # µm
        material_a      = str(data.get('material_a', ''))

        # ---- Gap / surroundings ------------------------------------------
        gap      = float(data.get('gap', 100.0))
        temp_b   = float(data.get('temp_b', 300.0))
        temp_surr= float(data.get('temp_surr', 300.0))
        material_b = str(data.get('material_b', ''))

        # ---- MC settings --------------------------------------------------
        n_photons = int(data.get('n_photons', 1000))
        if n_photons > 200000:
            n_photons = 200000

        full_gap_mc  = str(data.get('full_gap_mc', '')).lower() in ('true','1','on','yes')
        n_gap_photons = int(data.get('n_gap_photons', 1000))

        # ---- Wave-model solver selector (Phase 0/6) ------------------------
        # 'ray' (default fallback, Monte Carlo) or 'cached' (pre-computed full-wave).
        wave_model = str(data.get('wave_model', 'ray')).strip().lower()
        if wave_model not in ('ray', 'cached'):
            wave_model = 'ray'
        cache_path = str(data.get('cache_path', '')).strip()

        results = run_simulation(
            geometry_mode    = geometry_mode,
            cavity_diameter  = cavity_diameter,
            wall_thickness   = wall_thickness,
            packing_fraction = packing_fraction,
            height           = height,
            alpha_cnt        = alpha_cnt,
            alpha_ag         = alpha_ag,
            eps_flat_wall    = eps_flat_wall,
            width            = width,
            depth            = depth,
            cnt_dia_base     = cnt_dia_base,
            cnt_dia_top      = cnt_dia_top,
            cnt_pitch        = cnt_pitch,
            temp_a           = temp_a,
            emissivity_a     = emissivity_a,
            emissivity_a_back= emissivity_a_back,
            width_a          = width_a,
            depth_a          = depth_a,
            gap              = gap,
            temp_b           = temp_b,
            temp_surr        = temp_surr,
            n_photons        = n_photons,
            full_gap_mc      = full_gap_mc,
            n_gap_photons    = n_gap_photons,
            material_a       = material_a,
            material_b       = material_b,
            wave_model       = wave_model,
            cache_path       = cache_path,
        )
        return jsonify({'status': 'success', 'results': results})


    except Exception as e:
        import traceback
        return jsonify({'status': 'error', 'message': str(e),
                        'traceback': traceback.format_exc()}), 400


if __name__ == '__main__':
    app.run(debug=True, port=5000)
