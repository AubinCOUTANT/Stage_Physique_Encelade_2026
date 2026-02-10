# COUTANT Aubin L3 Internship: Modeling physical properties of Enceladus' icy surface
# Synthesis of static Tensile Strength (Y) as a function of grain size and temperature.
# Based on the cohesion model without sintering (Gundlach et al. 2018), applicable to cold plains.
# Comparison of two porosity scenarios: High porosity (phi=0.1) vs more compacted (phi=0.3).

# Start

import numpy as np
import matplotlib.pyplot as plt

# ==============================================================================
# 1. PHYSICAL PARAMETERS & MODELS
# ==============================================================================

R = 8.314                   # Ideal gas constant [J/(mol K)]
Q = 4600.0                  # Activation energy [J/mol] 

def get_gamma(T):
    """ Effective Surface Energy as a function of Temperature (Jabaud et al. 2024) """
    return 0.17 * np.exp(-Q / (R * T))

def tensile_strength(rp, T, phi_val):
    """ Formula: Y = (9 * phi * gamma(T)) / (2 * rp) """
    gamma_T = get_gamma(T)
    return (9 * phi_val * gamma_T) / (2 * rp)

# ==============================================================================
# 2. DATA SETUP
# ==============================================================================

rp_array = np.geomspace(0.5e-6, 300e-6, 1000)
Temps = [80, 120, 140, 160]
Colors = ['#00BFFF', '#FFC125', '#FF7F24', '#CD3333'] 
Labels = ['80 K (Plaines)', '120 K (Marges)', '140 K (Transition)', '160 K (Cœur Actif)']

Marker_Sizes = [1e-6, 10e-6, 100e-6] 
phi_list = [0.1, 0.3]  # The two cases to compare

# ==============================================================================
# 3. PLOTTING (Multi-panel)
# ==============================================================================

fig, axes = plt.subplots(2, 1, figsize=(11, 12), sharex=True)
plt.subplots_adjust(hspace=0.15)

for idx, phi_val in enumerate(phi_list):
    ax = axes[idx]
    
    for T_val, col, lab in zip(Temps, Colors, Labels):
        # 1. Plot the continuous line
        Y_vals = tensile_strength(rp_array, T_val, phi_val)
        ax.loglog(rp_array, Y_vals, color=col, linewidth=3.0, label=lab if idx==0 else "")
        
        # 2. Plot the markers
        for r_mark in Marker_Sizes:
            y_mark = tensile_strength(r_mark, T_val, phi_val)
            ax.scatter(r_mark, y_mark, color=col, edgecolor='black', 
                       s=120, marker='o', zorder=10, linewidth=1.5)

    # Subplot Formatting
    ax.set_title(f"Influence of Grain Size and Temperature ($\phi={phi_val}$)", fontsize=16, pad=10)
    ax.set_ylabel("Tensile Strength $\sigma$ [Pa]", fontsize=13, fontweight='bold')
    ax.grid(True, which='major', linestyle='-', alpha=0.5, color='gray')
    ax.grid(True, which='minor', linestyle=':', alpha=0.3, color='gray')
    ax.tick_params(which='both', direction='in', top=True, right=True, labelsize=11)
    
    # Text Annotations 
    y_ref_1um = tensile_strength(1e-6, 80, phi_val)
    ax.text(1.2e-6, y_ref_1um * 1.5, "Petits grains\n(Plus cohésifs)", color='#00BFFF', fontsize=10, fontweight='bold')
    
    y_ref_100um = tensile_strength(100e-6, 80, phi_val)
    ax.text(1.2e-4, y_ref_100um * 1.5, "Gros grains\n(Poudre libre)", color='#00BFFF', fontsize=10, fontweight='bold')

# Final formatting
axes[1].set_xlabel("Particle Radius $r_p$ [m]", fontsize=14, fontweight='bold')
axes[0].legend(loc='upper right', fontsize=10, framealpha=0.9, edgecolor='black', title="Régimes Thermiques")

plt.tight_layout()
plt.savefig('Fig_Tensile_Strength_Comparison.png', dpi=300)
plt.show()

# End
