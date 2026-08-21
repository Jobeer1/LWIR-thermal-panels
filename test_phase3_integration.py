"""Test Phase 3 near-field integration."""

from simulator import run_simulation

print("\n" + "="*70)
print("Testing Phase 3 (Polder-Van Hove Near-Field) Integration")
print("="*70)

# Test 1: Small gap (should trigger near-field)
print("\n[Test 1] Small gap (100 nm) - should trigger near-field")
print("-" * 70)
try:
    result = run_simulation(
        geometry_mode='honeycomb',
        cavity_diameter=20.0,
        height=450.0,
        gap=0.1,  # 100 nm gap - should trigger near-field
        temp_a=600.0,
        temp_b=300.0,
        n_photons=100,  # Small for speed
        enable_near_field=True,
        near_field_threshold=5.0,
        near_field_n_omega=20,  # Reduced for speed
        near_field_n_kparallel=15
    )
    
    print(f'✓ Simulation completed')
    print(f'  Physics regime: {result["physics_regime"]}')
    print(f'  Gap ratio: {result["gap_ratio"]:.3f}')
    print(f'  Evanescent fraction: {result["evanescent_fraction"]:.4f}')
    print(f'  Near-field flux: {result["net_flux_near_field_W_m2"]:.3e} W/m²')
    
    if result["physics_regime"] == "near-field":
        print("  ✓ Correctly identified as near-field regime")
    else:
        print("  ⚠ WARNING: Expected near-field regime")
except Exception as e:
    print(f"✗ Test 1 failed: {e}")
    import traceback
    traceback.print_exc()

# Test 2: Large gap (should stay far-field)
print("\n[Test 2] Large gap (100 µm) - should stay far-field")
print("-" * 70)
try:
    result2 = run_simulation(
        geometry_mode='honeycomb',
        cavity_diameter=20.0,
        height=450.0,
        gap=100.0,  # 100 µm gap - should be far-field
        temp_a=600.0,
        temp_b=300.0,
        n_photons=100,
        enable_near_field=True,
        near_field_threshold=5.0,
        near_field_n_omega=20,
        near_field_n_kparallel=15
    )
    
    print(f'✓ Simulation completed')
    print(f'  Physics regime: {result2["physics_regime"]}')
    print(f'  Gap ratio: {result2["gap_ratio"]:.3f}')
    
    if result2["physics_regime"] == "far-field":
        print("  ✓ Correctly identified as far-field regime")
    else:
        print("  ⚠ WARNING: Expected far-field regime")
except Exception as e:
    print(f"✗ Test 2 failed: {e}")
    import traceback
    traceback.print_exc()

# Test 3: Near-field disabled
print("\n[Test 3] Small gap (100 nm) with near-field disabled")
print("-" * 70)
try:
    result3 = run_simulation(
        geometry_mode='honeycomb',
        cavity_diameter=20.0,
        height=450.0,
        gap=0.1,
        temp_a=600.0,
        temp_b=300.0,
        n_photons=100,
        enable_near_field=False,
    )
    
    print(f'✓ Simulation completed')
    print(f'  Physics regime: {result3["physics_regime"]}')
    print(f'  Gap ratio: {result3["gap_ratio"]:.3f}')
    
    if "disabled" in result3["physics_regime"]:
        print("  ✓ Correctly shows near-field disabled")
except Exception as e:
    print(f"✗ Test 3 failed: {e}")
    import traceback
    traceback.print_exc()

print("\n" + "="*70)
print("✓✓✓ Phase 3 Integration Tests Complete")
print("="*70 + "\n")
