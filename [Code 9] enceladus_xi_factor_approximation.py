# COUTANT Aubin L3 Internship: Modeling physical properties of Enceladus' icy surface
# Numerical approximation of the geometric factor Xi (Ξ) for thermal conductivity.
# Extrapolation of reference data from Chan & Tien (1973) to high porosity 
# using a power law model.

# Start

import numpy as np
import matplotlib.pyplot as plt

# ==============================================================================
# 1. REFERENCE POINTS CALCULATION
# ==============================================================================

def calcul_exact_Xi_ref(phi, r=1.0):
    """
    Computes exact Xi values using Eq. 9 from Gundlach (2012) 
    and geometric data from Table 2 of Chan and Tien (1973).
    """
    Na, Nl, S_val = 0.0, 0.0, 0.0
    
    # Reference packing configurations (SC, BCC, FCC)
    if phi == 0.524: # Simple Cubic
        Na = 1.0 / (4.0 * r**2)
        Nl = 1.0 / (2.0 * r)
        S_val = 1.0
    elif phi == 0.68: # Body-Centered Cubic
        Na = 3.0 / (16.0 * r**2)
        Nl = np.sqrt(3.0) / (2.0 * r)
        S_val = 1.0 / 4.0
    elif phi == 0.74: # Face-Centered Cubic
        Na = 1.0 / (2.0 * np.sqrt(3.0) * r**2)
        Nl = np.sqrt(3.0 / 8.0) / r 
        S_val = 1.0 / 3.0
    else:
        raise ValueError("Unknown Phi value")

    constante = 0.5315
    # Exact Formula
    Xi = (1.0 / (constante * S_val)) * (Na / Nl)
    return Xi

# Generate the 3 base reference points
phi_refs = np.array([0.524, 0.68, 0.74])
xi_refs = np.array([calcul_exact_Xi_ref(p) for p in phi_refs])

# ==============================================================================
# 2. POWER LAW APPROXIMATION
# ==============================================================================
# We search for a model of the form: Xi = A * phi^B

log_phi = np.log(phi_refs)
log_xi  = np.log(xi_refs)

# Linear regression on logarithms
coeffs_log = np.polyfit(log_phi, log_xi, 1)

B = coeffs_log[0]           # The slope corresponds to the exponent
A = np.exp(coeffs_log[1])   # The intercept corresponds to the pre-factor A

def modele_puissance(phi):
    return A * (phi ** B)

# ==============================================================================
# 3. TARGET POINTS CALCULATION
# ==============================================================================
phi_cibles = np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6])
xi_cibles = modele_puissance(phi_cibles)

# ==============================================================================
# 4. RESULTS DISPLAY
# ==============================================================================
print("1. REFERENCE POINTS")
for p, x in zip(phi_refs, xi_refs):
    print(f"Phi: {p:.2f} -> Xi: {x:.4f}")

print(f"\n 2. FITTED PARAMETERS ")
print(f"A = {A:.4f}")
print(f"B = {B:.4f}")
print(f"Formula : Xi = {A:.4f} * phi^({B:.4f})")

print(f"\n 3. APPROXIMATED POINTS (Targets) ")
print(f"{'Phi':<10} | {'Xi (Approximated)':<15}")
print("-" * 30)
for p, x in zip(phi_cibles, xi_cibles):
    print(f"{p:<10.1f} | {x:<15.4f}")

# ==============================================================================
# 5. PLOTTING
# ==============================================================================
plt.figure(figsize=(10, 6))

# Smooth curve for the model
x_smooth = np.linspace(0.05, 0.8, 100)
y_smooth = modele_puissance(x_smooth)

plt.plot(x_smooth, y_smooth, 'b-', label=f"Power Law: $\Xi = {A:.2f} \phi^{{{B:.2f}}}$")

# Plot Reference Points
plt.plot(phi_refs, xi_refs, 'ro', markersize=10, label='Ref Points (Calculated)')

# Plot Target Points
plt.plot(phi_cibles, xi_cibles, 'gx', markersize=8, markeredgewidth=2, label='Target Points (0.1-0.6)')

plt.xlabel(r'$\phi$ (Solid Volume Fraction)')
plt.ylabel(r'Geometric Factor $\Xi$')
plt.title(r'Approximation via Power Law ($\Xi = A \cdot \phi^B$)')
plt.legend()
plt.grid(True, linestyle='--', alpha=0.6)
plt.show()

# End
