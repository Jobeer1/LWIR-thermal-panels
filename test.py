from simulator import run_simulation
r2 = run_simulation(geometry_mode='cnt_forest', height=10, alpha_cnt=0.98, alpha_ag=0.02,
                    cnt_dia_base=10, cnt_dia_top=10, cnt_pitch=50,
                    temp_a=600, emissivity_a=1.0, gap=100, temp_b=300,
                    n_photons=5000)
print(f"Kirchhoff error = {r2['kirchhoff_error']:.2f}%")
print(f"p_esc = {r2['p_esc']:.6f}")
print(f"alpha_eff = {r2['alpha_eff']:.4f}")
print(f"eps_b_raw = {r2['epsilon_b_raw']:.4f}")
print("--- Bug #13 wave-optics diagnostics ---")
print(f"cutoff lambda_c = {r2['cutoff_wavelength_um']:.4g} um")
print(f"propagating fraction (emit, analytic) = {r2['propagating_fraction']:.5f}")
print(f"deep-emission confinement = {r2['confinement_pct']:.2f} %")
print(f"escape solid angle = {r2['escape_solid_angle_sr']:.3e} sr")
print(f"decoupling ratio eps_B/alpha_eff = {r2['decoupling_ratio']:.4f}")
print(f"alpha_top (top surface) = {r2['alpha_top_surface']:.3f}")
print(f"LDOS emission gate G_em = {r2['confinement_gate']:.4f}")
print(f"evanescent rim fraction = {r2['rim_fraction']:.4f}")
print(f"Planck-averaged delta_ev = {r2['evanescent_decay_avg_um']:.3g} um")
print(f"epsB breakdown: cavity {r2['epsilon_b_cavity_part']:.4f} + flat {r2['epsilon_b_flat_part']:.4f}")

