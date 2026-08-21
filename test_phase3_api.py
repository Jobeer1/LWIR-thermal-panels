"""Test Phase 3 integration via Flask API."""

import json
import sys
import io
from app import app

# Fix encoding for Windows console
if sys.platform == 'win32':
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

print("\n" + "="*70)
print("Testing Phase 3 Integration via Flask API")
print("="*70)

with app.test_client() as client:
    
    # Test 1: API with near-field enabled for small gap
    print("\n[Test 1] API Request: Small gap (100 nm) with near-field enabled")
    print("-" * 70)
    
    request_data = {
        'geometry_mode': 'honeycomb',
        'cavity_diameter': 20.0,
        'height': 450.0,
        'gap': 0.1,  # 100 nm
        'temp_a': 600.0,
        'temp_b': 300.0,
        'n_photons': 50,
        'enable_near_field': 'true',
        'near_field_threshold': 5.0,
        'near_field_n_omega': 15,
        'near_field_n_kparallel': 10,
    }
    
    try:
        response = client.post('/api/simulate',
                              data=json.dumps(request_data),
                              content_type='application/json')
        
        if response.status_code == 200:
            result = response.get_json()
            if result['status'] == 'success':
                r = result['results']
                print(f"✓ API request successful")
                print(f"  Status: {result['status']}")
                print(f"  Physics regime: {r.get('physics_regime', 'N/A')}")
                print(f"  Gap ratio: {r.get('gap_ratio', 'N/A')}")
                print(f"  Evanescent fraction: {r.get('evanescent_fraction', 'N/A')}")
                print(f"  Near-field flux: {r.get('net_flux_near_field_W_m2', 'N/A')} W/m²")
            else:
                print(f"✗ API error: {result.get('message', 'Unknown error')}")
        else:
            print(f"✗ HTTP error: {response.status_code}")
            print(response.get_json())
    except Exception as e:
        print(f"✗ Exception: {e}")
        import traceback
        traceback.print_exc()
    
    # Test 2: API with near-field disabled
    print("\n[Test 2] API Request: Near-field disabled")
    print("-" * 70)
    
    request_data['enable_near_field'] = 'false'
    
    try:
        response = client.post('/api/simulate',
                              data=json.dumps(request_data),
                              content_type='application/json')
        
        if response.status_code == 200:
            result = response.get_json()
            if result['status'] == 'success':
                r = result['results']
                print(f"✓ API request successful")
                print(f"  Physics regime: {r.get('physics_regime', 'N/A')}")
                print(f"  Gap ratio: {r.get('gap_ratio', 'N/A')}")
                if 'disabled' in str(r.get('physics_regime', '')):
                    print("  ✓ Near-field correctly disabled")
            else:
                print(f"✗ API error: {result.get('message', 'Unknown error')}")
        else:
            print(f"✗ HTTP error: {response.status_code}")
    except Exception as e:
        print(f"✗ Exception: {e}")

print("\n" + "="*70)
print("✓✓✓ Phase 3 API Tests Complete")
print("="*70 + "\n")
