# COUTANT Aubin L3 Internship: Modeling physical properties of Enceladus' icy surface
# Validation: Saturation Vapor Pressure
# Model: Andreas (2007) / Murphy & Koop (2005)

# Start

import numpy as np
import matplotlib.pyplot as plt

# ==========================================
# 1. PHYSICAL LAW DEFINITION
# ==========================================

def p_sat_Andreas(T):
    """
    Computes the saturation vapor pressure of water ice (sublimation).
    Based on the formula from Andreas (2007), derived from Murphy & Koop (2005).
    
    Arguments:
        T (float or np.array): Temperature in Kelvin [K].
        
    Returns:
        P_sat (float or np.array): Saturation pressure in Pascals [Pa].
        
    Note on units:
    The original equation in the paper gives the pressure in hPa (with a 0.01 pre-factor).
    I removed this factor to get the result directly in Pascals (SI units),
    since 1 hPa = 100 Pa.
    """
    return np.exp(9.550426 - 5723.265/T + 3.53068*np.log(T) - 0.00728332*T)

# ==========================================
# 2. DATA GENERATION
# ==========================================

# Temperature range: from 60 K to 180 K 
T = np.linspace(60, 180, 500)
P = p_sat_Andreas(T)

# ==========================================
# 3. PLOTTING
# ==========================================

# Create figure
plt.figure(figsize=(8, 6))

# Main curve 
plt.semilogy(T, P, color='#00646E', linewidth=2.5, label='Andreas (2007) / Murphy & Koop Law')

# Highlight Enceladus Temperature Range (80 K - 160 K)
plt.axvspan(80, 160, color='grey', alpha=0.15, label='Enceladus T° Range')

# Annotations for key points (Plains, Tiger Stripes)
P_80 = p_sat_Andreas(80)
P_160 = p_sat_Andreas(160)

# Point at 80 K
plt.plot(80, P_80, 'ko')
plt.text(85, P_80*0.5, f'Plains (80 K)\n$P = {P_80:.2e} Pa$', fontsize=10, va='top')

# Point at 160 K
plt.plot(160, P_160, 'ko')
plt.text(155, P_160*2, f'Tiger Stripes (160 K)\n$P = {P_160:.2e} Pa$', fontsize=10, ha='right', va='bottom')

# Layout and Labels 
plt.xlabel("Temperature [K]", fontsize=12, fontweight='bold')
plt.ylabel("Saturation Vapor Pressure [Pa]", fontsize=12, fontweight='bold')
plt.title("Ice Thermodynamics (Andreas 2007 Model)", fontsize=13)
plt.grid(True, which="both", linestyle='--', alpha=0.6)
plt.legend(fontsize=11, loc='upper left')

# Set limits to focus on relevant data
plt.xlim(60, 180)
plt.ylim(1e-35, 1e0)

plt.tight_layout()

# Show
plt.show()

# End
