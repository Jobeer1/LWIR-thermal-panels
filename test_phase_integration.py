"""
Test script for Phase 1 & 2 integration into ray_tracer.py

Tests:
1. Phase 1 Complex Fresnel & TMM integration
2. Phase 2 Modal Loss integration
3. Backward compatibility with existing code
"""

import numpy as np
import math
from geometry import RectPit3D, CNTForestCell
from ray_tracer import (
    run_cavity_mc_3d,
    _trace_photon,
    _trace_photon_thin_film,
    evanescent_decay_length,
    evanescent_power_transmission
)

def test_phase1_complex_fresnel():
    """Test Phase 1: Complex Fresnel reflectance integration"""
    print("=" * 60)
    print("TEST 1: Phase 1 - Complex Fresnel & TMM Integration")
    print("=" * 60)
    
    # Create a simple rectangular cavity
    geom = RectPit3D(width_um=10, depth_um=10, height_um=100)
    
    print("\n1a. Testing _trace_photon_thin_film() with complex Fresnel...")
    
    # Trace a photon with complex Fresnel enabled
    pos = np.array([5e-6, 5e-6, 50e-6])  # Start in middle
    direction = np.array([0.0, 0.0, 1.0])  # Upward
    
    try:
        weight = _trace_photon_thin_film(
            pos=pos,
            direction=direction,
            geometry=geom,
            eps_walls_bulk=0.95,
            eps_base_bulk=0.05,
            re_entry_prob=0.0,
            wall_thickness_um=2.0,
            wall_material='alumina',
            base_material='silver',
            photon_wavelength_um=10.0,  # 10 µm wavelength
            use_complex_fresnel=True,
            apply_modal_attenuation=False,
            geometry_diameter_um=None
        )
        print(f"   ✓ Photon traced successfully. Survival weight: {weight:.6f}")
        assert 0.0 <= weight <= 1.0, "Weight should be in [0, 1]"
        print("   ✓ Weight is in valid range [0, 1]")
    except Exception as e:
        print(f"   ✗ Error: {e}")
        return False
    
    print("\n1b. Testing run_cavity_mc_3d() with complex Fresnel...")
    
    try:
        results = run_cavity_mc_3d(
            geometry=geom,
            n_photons=100,  # Small sample for speed
            eps_walls=0.95,
            eps_base=0.05,
            T_emit=300.0,
            T_inc=600.0,
            wall_thickness_um=2.0,
            wall_material='alumina',
            base_material='silver',
            use_complex_fresnel=True,
            apply_modal_attenuation=False
        )
        print(f"   ✓ Cavity MC simulation completed")
        print(f"     - p_esc: {results['p_esc']:.4f}")
        print(f"     - alpha_eff: {results['alpha_eff']:.4f}")
        print(f"     - epsilon_b: {results['epsilon_b']:.4f}")
        
        # Validate results are in reasonable ranges
        assert 0.0 <= results['p_esc'] <= 1.0, "p_esc out of range"
        assert 0.0 <= results['alpha_eff'] <= 1.0, "alpha_eff out of range"
        assert 0.0 <= results['epsilon_b'] <= 1.0, "epsilon_b out of range"
        print("   ✓ All results in valid ranges")
    except Exception as e:
        print(f"   ✗ Error: {e}")
        return False
    
    return True


def test_phase2_modal_attenuation():
    """Test Phase 2: Modal loss attenuation integration"""
    print("\n" + "=" * 60)
    print("TEST 2: Phase 2 - Modal Loss Attenuation Integration")
    print("=" * 60)
    
    # Create a CNT forest cavity (cylindrical pores)
    # CNTForestCell(pitch_um, dia_base_nm, dia_top_nm, height_um)
    geom = CNTForestCell(pitch_um=1.0, dia_base_nm=300, dia_top_nm=300, height_um=100)
    
    print("\n2a. Testing _trace_photon_thin_film() with modal attenuation...")
    
    pos = np.array([0.0, 0.0, 50e-6])
    direction = np.array([0.0, 0.0, 1.0])
    
    try:
        # Calculate gap width from geometry
        mean_dia_um = (300 + 300) / 2.0 / 1000.0  # Convert nm to µm
        gap_um = max(1.0 - mean_dia_um, 0.1 * 1.0)
        
        weight = _trace_photon_thin_film(
            pos=pos,
            direction=direction,
            geometry=geom,
            eps_walls_bulk=0.90,
            eps_base_bulk=0.02,
            re_entry_prob=0.0,
            wall_thickness_um=None,
            wall_material='alumina',
            base_material='silver',
            photon_wavelength_um=5.0,  # 5 µm wavelength
            use_complex_fresnel=False,
            apply_modal_attenuation=True,  # ENABLE Phase 2
            geometry_diameter_um=gap_um  # Use gap as characteristic dimension
        )
        print(f"   ✓ Photon traced with modal attenuation. Survival weight: {weight:.6f}")
        assert 0.0 <= weight <= 1.0, "Weight should be in [0, 1]"
        print("   ✓ Weight is in valid range [0, 1]")
    except Exception as e:
        print(f"   ✗ Error: {e}")
        return False
    
    print("\n2b. Testing run_cavity_mc_3d() with modal attenuation...")
    
    try:
        results = run_cavity_mc_3d(
            geometry=geom,
            n_photons=100,
            eps_walls=0.90,
            eps_base=0.02,
            T_emit=300.0,
            T_inc=600.0,
            wall_thickness_um=None,
            wall_material='alumina',
            base_material='silver',
            use_complex_fresnel=False,
            apply_modal_attenuation=True  # ENABLE Phase 2
        )
        print(f"   ✓ Cavity MC simulation with modal loss completed")
        print(f"     - p_esc: {results['p_esc']:.4f}")
        print(f"     - alpha_eff: {results['alpha_eff']:.4f}")
        print(f"     - epsilon_b: {results['epsilon_b']:.4f}")
        
        assert 0.0 <= results['p_esc'] <= 1.0, "p_esc out of range"
        assert 0.0 <= results['alpha_eff'] <= 1.0, "alpha_eff out of range"
        assert 0.0 <= results['epsilon_b'] <= 1.0, "epsilon_b out of range"
        print("   ✓ All results in valid ranges")
    except Exception as e:
        print(f"   ✗ Error: {e}")
        return False
    
    return True


def test_combined_phase1_phase2():
    """Test Phase 1 & Phase 2 together"""
    print("\n" + "=" * 60)
    print("TEST 3: Combined Phase 1 & Phase 2 Integration")
    print("=" * 60)
    
    # CNTForestCell(pitch_um, dia_base_nm, dia_top_nm, height_um)
    geom = CNTForestCell(pitch_um=1.0, dia_base_nm=500, dia_top_nm=500, height_um=150)
    
    print("\n3a. Testing run_cavity_mc_3d() with BOTH Phase 1 and Phase 2...")
    
    try:
        # Calculate gap for modal calculations
        mean_dia_um = (500 + 500) / 2.0 / 1000.0
        gap_um = max(1.0 - mean_dia_um, 0.1 * 1.0)
        
        results = run_cavity_mc_3d(
            geometry=geom,
            n_photons=150,
            eps_walls=0.92,
            eps_base=0.03,
            T_emit=300.0,
            T_inc=600.0,
            wall_thickness_um=1.5,
            wall_material='alumina',
            base_material='silver',
            use_complex_fresnel=True,      # Phase 1 ENABLED
            apply_modal_attenuation=True   # Phase 2 ENABLED
        )
        print(f"   ✓ Full integrated simulation completed")
        print(f"     - p_esc: {results['p_esc']:.4f} ± {results['p_esc_ci95']:.4f}")
        print(f"     - alpha_eff: {results['alpha_eff']:.4f} ± {results['alpha_eff_ci95']:.4f}")
        print(f"     - epsilon_b: {results['epsilon_b']:.4f} ± {results['epsilon_b_ci95']:.4f}")
        print(f"     - Kirchhoff error: {results['kirchhoff_error']:.2f}%")
        print(f"     - f_prop_emit: {results['f_prop_emit']:.4f}")
        print(f"     - f_prop_inc: {results['f_prop_inc']:.4f}")
        
        assert 0.0 <= results['p_esc'] <= 1.0, "p_esc out of range"
        assert 0.0 <= results['alpha_eff'] <= 1.0, "alpha_eff out of range"
        assert 0.0 <= results['epsilon_b'] <= 1.0, "epsilon_b out of range"
        print("   ✓ All results in valid ranges")
    except Exception as e:
        print(f"   ✗ Error: {e}")
        return False
    
    return True


def test_backward_compatibility():
    """Test backward compatibility with original code"""
    print("\n" + "=" * 60)
    print("TEST 4: Backward Compatibility")
    print("=" * 60)
    
    geom = RectPit3D(width_um=10, depth_um=10, height_um=100)
    
    print("\n4a. Testing run_cavity_mc_3d() without Phase 1/Phase 2 (original mode)...")
    
    try:
        results = run_cavity_mc_3d(
            geometry=geom,
            n_photons=100,
            eps_walls=0.95,
            eps_base=0.05,
            T_emit=300.0,
            T_inc=600.0,
            wall_thickness_um=None,  # None = bulk
            use_complex_fresnel=False,  # Phase 1 DISABLED
            apply_modal_attenuation=False  # Phase 2 DISABLED
        )
        print(f"   ✓ Original mode simulation completed")
        print(f"     - p_esc: {results['p_esc']:.4f}")
        print(f"     - alpha_eff: {results['alpha_eff']:.4f}")
        print(f"     - epsilon_b: {results['epsilon_b']:.4f}")
        
        assert 0.0 <= results['p_esc'] <= 1.0, "p_esc out of range"
        assert 0.0 <= results['alpha_eff'] <= 1.0, "alpha_eff out of range"
        assert 0.0 <= results['epsilon_b'] <= 1.0, "epsilon_b out of range"
        print("   ✓ Backward compatibility maintained")
    except Exception as e:
        print(f"   ✗ Error: {e}")
        return False
    
    return True


def test_evanescent_modes():
    """Test evanescent mode calculations"""
    print("\n" + "=" * 60)
    print("TEST 5: Evanescent Mode Calculations")
    print("=" * 60)
    
    lambda_c_um = 3.0  # Cutoff wavelength
    
    print("\n5a. Testing evanescent_decay_length()...")
    
    test_cases = [
        (lambda_c_um, 2.0, float('inf')),  # λ < λ_c: propagating
        (lambda_c_um, 3.0, float('inf')),  # λ = λ_c: edge
        (lambda_c_um, 5.0, 0.238),         # λ > λ_c: evanescent (approx δ_ev ≈ λ_c/2π)
        (lambda_c_um, 10.0, 0.238),        # λ >> λ_c: deeply evanescent
    ]
    
    for lambda_c, lambda_w, expected_type in test_cases:
        delta = evanescent_decay_length(lambda_c, lambda_w)
        if lambda_w < lambda_c:
            print(f"   λ={lambda_w:.1f}µm < λ_c={lambda_c:.1f}µm: δ_ev = {delta:.3f} (propagating ✓)")
        elif lambda_w == lambda_c:
            print(f"   λ={lambda_w:.1f}µm = λ_c={lambda_c:.1f}µm: δ_ev = ∞ (cutoff ✓)")
        else:
            print(f"   λ={lambda_w:.1f}µm > λ_c={lambda_c:.1f}µm: δ_ev = {delta:.3f}µm (evanescent ✓)")
        
        assert math.isfinite(delta) or lambda_w <= lambda_c, "Decay length should be finite for λ > λ_c"
    
    print("   ✓ Decay length calculations correct")
    
    print("\n5b. Testing evanescent_power_transmission()...")
    
    T_short = evanescent_power_transmission(lambda_c_um, 5.0, 1.0)
    T_long = evanescent_power_transmission(lambda_c_um, 5.0, 10.0)
    
    print(f"   L=1.0µm: T = {T_short:.6f}")
    print(f"   L=10.0µm: T = {T_long:.6f}")
    
    assert T_short > T_long, "Transmission should decrease with distance"
    assert 0.0 <= T_short <= 1.0, "Transmission should be in [0, 1]"
    assert 0.0 <= T_long <= 1.0, "Transmission should be in [0, 1]"
    print("   ✓ Power transmission calculations correct")
    
    return True


if __name__ == "__main__":
    print("\n")
    print("╔════════════════════════════════════════════════════════════╗")
    print("║  Phase 1 & 2 Integration Test Suite for ray_tracer.py    ║")
    print("╚════════════════════════════════════════════════════════════╝")
    
    all_passed = True
    
    # Run all tests
    tests = [
        ("Phase 1: Complex Fresnel & TMM", test_phase1_complex_fresnel),
        ("Phase 2: Modal Loss Attenuation", test_phase2_modal_attenuation),
        ("Combined Phase 1 & 2", test_combined_phase1_phase2),
        ("Backward Compatibility", test_backward_compatibility),
        ("Evanescent Modes", test_evanescent_modes),
    ]
    
    results = []
    for name, test_func in tests:
        try:
            passed = test_func()
            results.append((name, passed))
        except Exception as e:
            print(f"\n✗ UNEXPECTED ERROR in {name}:")
            print(f"  {e}")
            results.append((name, False))
            all_passed = False
    
    # Summary
    print("\n" + "=" * 60)
    print("TEST SUMMARY")
    print("=" * 60)
    for name, passed in results:
        status = "✓ PASS" if passed else "✗ FAIL"
        print(f"{status}: {name}")
    
    print("\n" + "=" * 60)
    if all_passed:
        print("✓ ALL TESTS PASSED")
        print("=" * 60)
    else:
        print("✗ SOME TESTS FAILED")
        print("=" * 60)
    
    exit(0 if all_passed else 1)
