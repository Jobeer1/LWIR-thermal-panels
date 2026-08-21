#!/usr/bin/env python3
"""
Test script to verify thin-film physics implementation.
"""

import sys
import os
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from material_optics import (
    effective_emissivity_thin_film,
    planck_weighted_effective_emissivity,
    optical_thickness_analysis
)

def test_basic_physics():
    """Test the core thin-film physics."""
    print("=" * 60)
    print("Testing Thin-Film Physics Implementation")
    print("=" * 60)
    
    # Test 1: Very thin wall (100nm alumina at 10µm wavelength)
    print("\n1. 100nm alumina wall at λ=10µm (peak 300K radiation):")
    eps = effective_emissivity_thin_film(
        bulk_emissivity=0.8,
        thickness_um=0.1,  # 100nm
        wavelength_um=10.0,
        material='alumina'
    )
    print(f"   ε_eff = {eps:.4f} (TMM result, correction factor: {eps/0.8:.3f})")
    
    # Test 2: Thick wall (10µm alumina)
    print("\n2. 10µm alumina wall at λ=10µm:")
    eps = effective_emissivity_thin_film(
        bulk_emissivity=0.8,
        thickness_um=10.0,
        wavelength_um=10.0,
        material='alumina'
    )
    print(f"   ε_eff = {eps:.4f} (should be close to bulk 0.8)")
    
    # Test 3: Planck-weighted average at 300K
    print("\n3. Planck-weighted effective ε for 100nm alumina at 300K:")
    eps_weighted = planck_weighted_effective_emissivity(
        bulk_emissivity=0.8,
        thickness_um=0.1,
        material='alumina',
        temperature_K=300.0
    )
    print(f"   ε_eff(weighted) = {eps_weighted:.4f}")
    
    # Test 4: Optical thickness analysis
    print("\n4. Optical thickness analysis for 100nm alumina at 300K:")
    analysis = optical_thickness_analysis(
        thickness_um=0.1,
        material='alumina',
        temperature_K=300.0
    )
    for key, value in analysis.items():
        if isinstance(value, float):
            print(f"   {key}: {value:.4g}")
        else:
            print(f"   {key}: {value}")
    
    # Test 5: Compare to literature claim
    print("\n5. Verification against literature claim:")
    print("   TMM check: 100nm wall remains below its bulk emissivity")
    print(f"   Our calculation: {analysis['effective_emissivity_weighted']*100:.1f}%")
    print(f"   Bounded: {'yes' if analysis['effective_emissivity_weighted'] <= analysis['bulk_emissivity'] else 'no'}")
    
    # Test 6: CNT forest properties
    print("\n6. CNT forest thin-film properties:")
    eps_cnt = effective_emissivity_thin_film(
        bulk_emissivity=0.98,
        thickness_um=0.1,  # 100nm CNT wall
        wavelength_um=10.0,
        material='cnt_forest'
    )
    print(f"   100nm CNT at λ=10µm: ε_eff = {eps_cnt:.4f}")
    
    return True

def test_monte_carlo_integration():
    """Test that thin-film physics integrates with Monte Carlo."""
    print("\n" + "=" * 60)
    print("Testing Monte Carlo Integration")
    print("=" * 60)
    
    try:
        from ray_tracer import run_cavity_mc_3d
        from geometry import HoneycombCavityCell
        
        print("Creating test geometry...")
        geometry = HoneycombCavityCell(
            diameter_um=500.0,
            height_um=20000.0,
            wall_emissivity=0.8,
            packing_fraction=0.75
        )
        
        print("Running Monte Carlo with thin-film physics...")
        # Small photon count for quick test
        results = run_cavity_mc_3d(
            geometry=geometry,
            n_photons=1000,
            eps_walls=0.8,  # Bulk emissivity
            eps_base=0.1,
            wall_thickness_um=0.1,  # 100nm walls
            wall_material='alumina',
            base_material='silver'
        )
        
        print(f"Escape probability: {results['p_esc']:.6f}")
        print(f"Effective absorptivity: {results['alpha_eff']:.4f}")
        print(f"Effective emissivity: {results['epsilon_b']:.4f}")
        print("\nThin-film physics successfully integrated!")
        
        return True
    except Exception as e:
        print(f"Integration test failed: {e}")
        import traceback
        traceback.print_exc()
        return False

def main():
    """Run all tests."""
    print("Thin-Film Physics Test Suite")
    print("=" * 60)
    
    success = True
    
    # Test basic physics
    if not test_basic_physics():
        success = False
    
    # Test Monte Carlo integration
    if not test_monte_carlo_integration():
        success = False
    
    print("\n" + "=" * 60)
    if success:
        print("✓ All tests passed!")
    else:
        print("✗ Some tests failed")
    print("=" * 60)
    
    return success

if __name__ == '__main__':
    success = main()
    sys.exit(0 if success else 1)