p = r'static/app.js'
L = open(p, encoding='utf-8', errors='replace').readlines()

# Fix lines 138-148 (0-indexed 137-147) with clean content
# (verified against current app.js line numbers)
L[137] = '    _setText("res-view-factor", vf);\n'
L[138] = '    _setText("res-total-heat", qNet);\n'
L[139] = '    _setText("res-coupling", _fmtW(r.q_af_to_b, 6));\n'
L[140] = '    _setText("res-aback", _fmtW(r.q_ab_net, 6));\n'
L[141] = '    // Physical net flux A->B (accounts for eps_A != eps_B)\n'
L[142] = '    _setText("res-flux-afront", _fmtW(r.q_net_a_to_b_physical != null ? r.q_net_a_to_b_physical : r.net_flux_A_front, 1));\n'
L[143] = '    // Stagnation temperature (adiabatic wall temperature of Plate B)\n'
L[144] = '    if (r.T_B_stag != null) {\n'
L[145] = '        _setText("res-t-stag", (Number(r.T_B_stag) - 273.15).toFixed(1) + String.fromCharCode(176) + "C  (" + Number(r.T_B_stag).toFixed(1) + " K)");\n'
L[146] = '    }\n'
L[147] = '    // Bug #13 - wave-optics decoupling diagnostics\n'

open(p, 'w', encoding='utf-8').writelines(L)
print('Fixed JS lines 137-147')
for i in range(136, 148):
    print(f'{i+1:4d}: {L[i]}', end='')
