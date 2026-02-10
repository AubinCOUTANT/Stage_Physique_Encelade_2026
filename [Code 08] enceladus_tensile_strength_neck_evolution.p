# COUTANT Aubin L3 Internship: Modeling physical properties of Enceladus' icy surface
# Dynamic modeling of Tensile Strength (Y) evolution driven by sintering neck growth.
# Focus on "Hot Ice" regions (Active Zones: 120 K - 160 K).


# Start

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec # Module for custom subplot grid layout

# ==============================================================================
# 1. PHYSICAL PARAMETERS & CONSTANTS
# ==============================================================================

R_gas = 8.314               # Ideal gas constant [J/(mol K)]
Q = 4600.0                  # Activation energy [J/mol] (Value for surface diffusion)

# Material Parameters 
E = 10.5e9                  # Young's Modulus [Pa]
nu = 0.31                   # Poisson's ratio
xi = 0.1                    # Geometric scaling factor (Calibration parameter)
phi = 0.1                   # Filling volume factor

# Derived Elastic Constant K
K = (2/3) * E / (1 - nu**2)

# ==============================================================================
# 2. PHYSICAL MODELS
# ==============================================================================

def get_gamma(T_val):
    """ Effective Surface Energy as a function of Temperature [J/m^2] """
    return 0.17 * np.exp(-Q / (R_gas * T_val))

def get_rn_JKR(rp, g_val):
    """ 
    Calculates initial adhesion contact radius (JKR Theory).
    This is the starting point before any sintering occurs.
    """
    term = (6 * np.pi * g_val * (rp/2)**2) / K
    return term**(1/3)

def res_tract(rn, rp, g_val):
    """ 
    Calculates Tensile Strength Y [Pa] as a function of neck radius rn.
    Formula based on Gundlach et al. (2018).
    """
    beta = np.sqrt(1.5 * np.pi * g_val * K)
    return xi * (3 * phi / (2 * np.pi)) * beta * (rn**1.5 / rp**2)

# ==============================================================================
# 3. PLOTTING CONFIGURATION 
# ==============================================================================

Grain_Sizes = [1e-6, 5e-6, 10e-6, 20e-6, 50e-6, 100e-6]
Labels = ['1 µm', '5 µm', '10 µm', '20 µm', '50 µm', '100 µm']
Temps_panel = [120.0, 140.0, 160.0]
colors = plt.cm.viridis(np.linspace(0, 0.9, len(Grain_Sizes)))

fig = plt.figure(figsize=(14, 11)) 
# GridSpec 2 rows x 4 columns to center the bottom plot
gs = gridspec.GridSpec(2, 4, figure=fig, hspace=0.35, wspace=0.4)

ax1 = fig.add_subplot(gs[0, 0:2]) # Top Left
ax2 = fig.add_subplot(gs[0, 2:4]) # Top Right
ax3 = fig.add_subplot(gs[1, 1:3]) # Bottom Center

axes_list = [ax1, ax2, ax3]

# Loop over temperatures (panels)
for idx, (ax, T) in enumerate(zip(axes_list, Temps_panel)):
    current_gamma = get_gamma(T)
    
    # Loop over grain sizes (curves)
    for rp, color, label in zip(Grain_Sizes, colors, Labels):
        rn_start = get_rn_JKR(rp, current_gamma)
        rn_end = rp 
        
        # Logarithmic space for neck radius evolution
        rn_vec = np.geomspace(rn_start, rn_end, 300)
        Y_vec = res_tract(rn_vec, rp, current_gamma)
        
        # Legend labels attached only to the bottom graph (idx == 2)
        lbl = f"$r_p =$ {label}" if idx == 2 else ""
        
        ax.loglog(rn_vec, Y_vec, linewidth=2.5, color=color, label=lbl)
        
        # Mark the initial state (JKR Adhesion)
        ax.scatter(rn_start, Y_vec[0], color=color, s=60, marker='o', edgecolors='black', zorder=5)

    # Axis Labels & Title
    ax.set_ylabel("Tensile Strength $Y$ [Pa]", fontsize=12, fontweight='bold')
    ax.set_xlabel("Sintering Pont Radius $r_n$ [m]", fontsize=12, fontweight='bold')
    ax.set_title(f"Regime: $T = {int(T)}$ K ( $\\phi={phi}$ )", fontsize=14, fontweight='bold', loc='center')
    
    # Grid & Ticks
    ax.grid(True, which="major", linestyle='-', alpha=0.4, color='gray')
    ax.grid(True, which="minor", linestyle=':', alpha=0.2, color='gray')
    
    # Info Box (Parameters)
    info_box = (f"PARAMETERS:\n$\\Phi= {phi}$\n$T = {T}$ K\n$\\gamma \\approx {current_gamma:.2e}$ J/m²")
    ax.text(0.02, 0.95, info_box, transform=ax.transAxes, fontsize=10, va='top', 
            bbox=dict(facecolor='white', alpha=0.9, edgecolor='black', boxstyle='round,pad=0.5'))

    ax.tick_params(which='both', direction='in', top=True, right=True, labelsize=11)
    ax.set_ylim(1e-2, 1e7)
    ax.set_xlim(1e-9, 1e-4)

# Place Legend to the right of the bottom plot 
axes_list[2].legend(title="Grain Size $r_p$", fontsize=11, loc='center left', bbox_to_anchor=(1.05, 0.5), borderaxespad=0.)
plt.show()

# End
