# COUTANT Aubin L3 Internship: Modeling physical properties of Enceladus' icy surface
# Dynamic modeling of Thermal Conductivity evolution driven by sintering.
# Uses the geometric factor Xi calculated in [Code 9] for high porosity media (phi=0.1).
# Application to "Hot Ice" regions (T=120K).

# Start

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

# ==============================================================================
# 1. PHYSICAL PARAMETERS & CONSTANTS
# ==============================================================================

# Thermodynamic Parameters (HOT ICE)
T = 120.0                   # Surface Temperature [K]
R_gas = 8.314               # Ideal gas constant [J/(mol K)]
Q = 4600.0                  # Activation energy [J/mol]

# T-dependent properties calculation
gamma_val = 0.17 * np.exp(-Q / (R_gas * T)) # Surface Energy [J/m^2]
g1 = 567.0                                  # Klinger's Constant
lambda_bulk = g1 / T                        # Bulk Ice Thermal Conductivity [W/mK]

# Material Parameters (Ice)
E = 10.5e9                  # Young's Modulus [Pa]
nu = 0.31                   # Poisson's Ratio
phi = 0.1                   # Filling volume factor

# Geometric Factor Xi for phi = 0.1
# Note: Value determined by numerical approximation [Code 9]
Xi_val = 0.0088

# Derived Elastic Constant K = 2/3 * E / (1 - nu^2)
K = (2/3) * E / (1 - nu**2)

# Sintering Constant Beta (Eq 12 Gundlach 2018)
# Beta = sqrt(3/2 * pi * gamma * K)
beta = np.sqrt(1.5 * np.pi * gamma_val * K)

# ==============================================================================
# 2. PHYSICAL MODEL (Gundlach 2018)
# ==============================================================================

def JKR_theorie(rp):
    """
    Calculates initial contact radius (JKR Neck) for a grain of radius rp.
    Formula: a0^3 = (6 * pi * gamma * (rp/2)^2) / K
    """
    term = (6 * np.pi * gamma_val * (rp/2)**2) / K
    return term**(1/3)

def conductivity_sintering(rn, rp):
    """
    Calculates Thermal Conductivity Lambda as a function of neck radius rn.
    Based on Equation 20 from Gundlach et al. (2018).
    
    Lambda = Lambda_Bulk * [ (3/4) * ((1-nu^2)/E) * Beta ]^(1/3) * Xi * (rn^(1/2) / rp^(2/3))
    """
    # Hertzian term (constant bracket to power 1/3)
    hertz_term_bracket = (3.0 / 4.0) * ((1 - nu**2) / E) * beta
    hertz_factor = hertz_term_bracket**(1.0/3.0)
    
    # Full scaling law
    # Lambda scales as rn^0.5 / rp^0.66
    lam = lambda_bulk * hertz_factor * Xi_val * (np.sqrt(rn) / (rp**(2.0/3.0)))
    
    return lam

# ==============================================================================
# 3. SIMULATION CONFIGURATION
# ==============================================================================

# Grain radii to plot
Grain_Sizes = [1e-6, 5e-6, 10e-6, 20e-6, 50e-6, 100e-6]
Labels = ['1 µm', '5 µm', '10 µm', '20 µm', '50 µm', '100 µm']

# Distinct color palette
colors = ['#1f77b4', '#d62728', '#2ca02c', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f']

# ==============================================================================
# 4. PLOTTING
# ==============================================================================

plt.figure(figsize=(11, 7))

for rp, color, label in zip(Grain_Sizes, colors, Labels):
    
    # 1. Define Sintering Domain
    rn_start = JKR_theorie(rp)  # Physical starting point (JKR)
    rn_end = rp                 # Physical limit (fully sintered)

    if rn_start >= rn_end: rn_end = rp # Safety check

    # 2. Compute Values
    rn_vec = np.geomspace(rn_start, rn_end, 300)
    Y_vec = conductivity_sintering(rn_vec, rp)
    
    # 3. Plot Curve
    plt.loglog(rn_vec, Y_vec, linewidth=2.5, color=color, label=f"$r_p =$ {label}")
    
    # 4. Mark Initial Point (State at 120K without sintering)
    plt.scatter(rn_start, Y_vec[0], color=color, s=50, marker='o', edgecolors='black', zorder=5)

# ==============================================================================
# 5. FORMATTING & LEGEND
# ==============================================================================

plt.xlabel("Sintering Neck Radius $r_n$ [m]", fontsize=14)
plt.ylabel("Thermal Conductivity $\lambda$ [W m$^{-1}$ K$^{-1}$]", fontsize=14)
plt.title(f"Thermal Conductivity Evolution with Neck Growth\n Effect of Grain Size ($T={T}$ K, $\phi={phi}$)", fontsize=16)

plt.grid(True, which="major", linestyle='-', alpha=0.5, color='gray')
plt.grid(True, which="minor", linestyle=':', alpha=0.3, color='gray')

plt.xlim(1e-9, 2e-4)

# Parameter Info Box
info_text = (f"PARAMETERS:\n"
             f"$T = {T}$ K\n"
             f"$\phi = {phi}$\n"
             f"$\\lambda_{{bulk}} \\approx {lambda_bulk:.2f}$ W/mK\n"
             f"$\\gamma \\approx {gamma_val:.2e}$ J/m²")

plt.text(0.02, 0.95, info_text, transform=plt.gca().transAxes, fontsize=15, va='top',
         bbox=dict(facecolor='white', alpha=0.9, edgecolor='black', boxstyle='round,pad=0.5'))

ax = plt.gca()
ax.tick_params(which='both', direction='in', top=True, right=True, labelsize=12)
ax.tick_params(which='major', length=8, width=1.5)
ax.tick_params(which='minor', length=4, width=1)
loc_majeurs_y = ticker.LogLocator(base=10.0, numticks=10)
ax.yaxis.set_major_locator(loc_majeurs_y)

plt.legend(title="Grain Size $r_p$", fontsize=11, loc='upper left', bbox_to_anchor=(1, 1))
plt.tight_layout()
plt.show()

# End
