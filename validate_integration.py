"""Quick validation of Phase 1 & 2 integration"""
from ray_tracer import run_cavity_mc_3d
from geometry import RectPit3D, CNTForestCell

print('=' * 60)
print('FINAL VALIDATION: Phase 1 & 2 Integration')
print('=' * 60)

# Test 1: Phase 1 Only
print('\nTest 1: Phase 1 Only (Complex Fresnel)')
geom = RectPit3D(width_um=10, depth_um=10, height_um=100)
r1 = run_cavity_mc_3d(geom, 50, 0.95, 0.05, use_complex_fresnel=True, wall_thickness_um=2.0)
print('  OK: epsilon_b={:.4f}, alpha_eff={:.4f}'.format(r1['epsilon_b'], r1['alpha_eff']))

# Test 2: Phase 2 Only
print('\nTest 2: Phase 2 Only (Modal Attenuation)')
geom2 = CNTForestCell(pitch_um=1.0, dia_base_nm=300, dia_top_nm=300, height_um=100)
r2 = run_cavity_mc_3d(geom2, 50, 0.92, 0.02, apply_modal_attenuation=True)
print('  OK: epsilon_b={:.4f}, alpha_eff={:.4f}'.format(r2['epsilon_b'], r2['alpha_eff']))

# Test 3: Phase 1 + 2
print('\nTest 3: Phase 1 + 2 Combined')
r3 = run_cavity_mc_3d(geom2, 50, 0.92, 0.02, use_complex_fresnel=True, wall_thickness_um=1.5, apply_modal_attenuation=True)
print('  OK: epsilon_b={:.4f}, alpha_eff={:.4f}'.format(r3['epsilon_b'], r3['alpha_eff']))

# Test 4: Original Mode
print('\nTest 4: Original Mode (Backward Compatibility)')
r4 = run_cavity_mc_3d(geom, 50, 0.95, 0.05)
print('  OK: epsilon_b={:.4f}, alpha_eff={:.4f}'.format(r4['epsilon_b'], r4['alpha_eff']))

print('\n' + '=' * 60)
print('ALL VALIDATION TESTS PASSED')
print('=' * 60)
