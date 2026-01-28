# COUTANT Aubin L3 Internship: Modeling physical properties of Enceladus' icy surface
# Figure 8 Reproduction: Sintering vs Sublimation Competition
# Model: Gundlach et al. (2018)

# Start

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

# ==========================================
# 1. NUMERICAL SOLVER 
# ==========================================

def runge_kutta_4_vec(f, y0, a, b, N):
    """
    Vectorized version of the Runge-Kutta 4.
    
    NOTE: While the standard RK4 was studied in the L3 Numerical Analysis course,
    this version is adapted to handle vector states Y = [r_p, r_n].
    It allows solving the coupled evolution of particle radius and neck radius simultaneously.
    """
    h = (b - a) / N
    T = np.linspace(a, b, N + 1)
    y0 = np.asarray(y0)
    
    # Initialization of the results matrix 
    Y = np.zeros((N + 1,) + y0.shape)
    Y[0] = y0
    
    for i in range(N):
        k1 = h * f(T[i], Y[i])
        k2 = h * f(T[i] + h/2, Y[i] + 0.5 * k1)
        k3 = h * f(T[i] + h/2, Y[i] + 0.5 * k2)
        k4 = h * f(T[i] + h, Y[i] + k3)
        
        Y[i+1] = Y[i] + (1/6) * (k1 + 2*k2 + 2*k3 + k4)
        
    return T, Y

# ==========================================
# 2. PARAMETERS AND CONSTANTS 
# ==========================================

R = 8.31447        # Ideal gas constant [J/(mol·K)]
mu = 18.015e-3     # Molar mass of water [kg/mol]
gamma = 0.2        # Specific surface energy [J/m^2]
Omega = 2.0e-5     # Molar volume [m^3/mol]
alpha = np.pi/2    # Geometric packing factor related to porosity

# Pressure condition in the experimental chamber
p_gas = 1e-5       # Partial pressure of water [Pa] (Adjusted to fit data)
T_gas = 300.0      # Chamber wall temperature [K]

def rho_ice(T):
    """ Ice density as a function of Temperature [kg/m^3] """
    return 916 - 0.175 * T - 5.0e-4 * T**2

def p_sat(T):
    """
    Computes saturation vapor pressure (Andreas 2007 / Murphy & Koop 2005).
    """
    return np.exp(9.550426 - 5723.265/T + 3.53068*np.log(T) - 0.00728332*T)

# ===================================
# 3. PHYSICAL MODEL 
# ===================================

def hertz_knudsen(T_surface, p_gas_val, T_gas_val):
    """ 
    Calculates the net MASS flux Z (Eq. 8 in Gundlach et al., 2018).
    Includes the dimensional correction (multiplication by mu inside sqrt).
    """
    # Outgoing flux (Sublimation) 
    flux_out = p_sat(T_surface) * np.sqrt(mu / (2 * np.pi * R * T_surface))
    
    # Incoming flux (Condensation from background gas)
    flux_in = p_gas_val * np.sqrt(mu / (2 * np.pi * R * T_gas_val))
    
    return (flux_out - flux_in)

def system_func(T_pa, T_mt, T_su):
    """
    Defines the system of coupled differential equations:
    dy/dt = f(t, y) where y = [r_p, r_n]
    
    Parameters:
    - T_pa: Particle Temperature (controls particle sublimation)
    - T_mt: Mass Transport Temperature (controls diffusion/sintering)
    - T_su: Neck Surface Temperature (controls neck sublimation)
    """
    rho = rho_ice(T_pa) 
    
    # Calculate fluxes based on Hertz-Knudsen
    Z_particule = hertz_knudsen(T_pa, p_gas, T_gas) # For particle erosion
    Z_neck = hertz_knudsen(T_su, p_gas, T_gas)      # For neck erosion
    
    # Sintering efficiency factor Zeta (Eq. 3 in Gundlach)
    P_mt = p_sat(T_mt)
    zeta = (Omega**2 * gamma * P_mt) / (R * T_mt * np.sqrt(2 * np.pi * mu * R * T_mt))
    
    def f(t, Y):
        rp = Y[0] # Current particle radius
        rn = Y[1] # Current neck radius
        
        # 1. Particle Evolution
        # Critical radius 
        rc_val = 2 * mu * gamma / (rho * R * T_pa)
        
        # Evolution of particle radius (Eq. 7 derivative)
        facteur_courbure = rp / (rp - rc_val) 
        drp = - (Z_particule / rho) * facteur_courbure 
        
        # 2. Neck Evolution 
        
        # Geometry factors for Sintering (Eq. 4, 5, 6)
        denom = 2 * (rp - rn)
        if denom <= 0: denom = 1e-15 # Avoid division by zero
        
        Delta = rn**2 / denom
        phi = np.arcsin(rn / rp)
        
        d_source = rp * (alpha/2 + phi - np.pi/2)
        d_sink = Delta * np.arctan(rp / (rn + Delta))
        S = d_source / (d_source - d_sink)
        
        # Growth term: Sintering (Eq. 2, first term)
        # Driven by curvature difference 
        terme_crochet = (2/rp) + (1/Delta) - (1/rn)
        dv_sinter = zeta * S * terme_crochet
        
        # Loss term: Sublimation of the neck (Eq. 9)
        dv_sub_neck = - Z_neck / rho 
        
        # Master Equation (Eq. 1): Total neck evolution
        drn = dv_sinter + dv_sub_neck
        
        # Physical Constraint: Neck cannot exceed Particle size
        if rn >= rp: 
            drn = drp 
            
        return np.array([drp, drn])
    return f

# ==========================================
# 4. SIMULATION
# ==========================================

# Initial conditions [r_p0, r_n0]
Y0 = [1.5e-6, 7.5e-8]  

# Time parameters
t_max = 10580     
N_pas = 200000      

# Temperature Setup
# NOTE: T_su is adjusted (-9.5 K) to account for self-shadowing/concavity 
# protection of the neck
T_exp = 160.0
T_pa = T_exp        # Particle T
T_mt = T_exp        # Diffusion T
T_su = T_exp - 9.5  # Neck T (Cooler)

# Run Simulation
fonction_systeme = system_func(T_pa, T_mt, T_su)
Temps, Resultats = runge_kutta_4_vec(fonction_systeme, Y0, 0, t_max, N_pas)

# Extract results
Rayon_P = Resultats[:, 0]
Rayon_N = Resultats[:, 1]

# ==========================================
# 5. PLOTTING
# ==========================================

plt.figure(figsize=(10, 7))

# Plotting curves
plt.loglog(Temps, Rayon_N, 'b-', linewidth=2.5, label='Sinter Neck Radius ($r_n$)')
plt.loglog(Temps, Rayon_P, color='tab:blue', linestyle='--', linewidth=2, label='Particle Radius ($r_p$)')

# Limit lines
plt.axhline(y=1e-7, color='black', linestyle='-.', linewidth=1)
plt.axvline(x=9e3, color='gray', linestyle=':', alpha=0.5)

# Annotations (matching Gundlach Fig 8 style)
plt.text(2.0e3, 4.2e-7, "Neck\nEvolution\nStage", fontsize=13, ha='center', color='black')
plt.text(2.5e4, 4.2e-7, "Ice\nSolidification\nStage", fontsize=13, ha='center', color='black')
plt.text(100, 2.5e-6, "T = 160.0 K", fontsize=20, ha='center', color='black')
plt.text(4, 2.5e-6, "$r_p$ = 1.5 µm", fontsize=20, ha='center', color='black')
plt.text(1.1, 1.1e-7, "Resolution Limit of the Cryo-SEM", fontsize=13, ha='left', color='black')

# Styling
plt.title("Reproduction of Theoretical Models (Fig. 8 Gundlach 2018)", fontsize=13)
plt.xlabel("Time [s]", fontsize=12, fontweight='bold')
plt.ylabel("Sinter Neck Radius, Particle Radius [m]", fontsize=12, fontweight='bold')

ax = plt.gca()
ax.tick_params(which='both', direction='in', top=True, right=True, labelsize=12)
ax.tick_params(which='major', length=8, width=1.5)
ax.tick_params(which='minor', length=4, width=1)
loc_majeurs_y = ticker.LogLocator(base=10.0, numticks=10)
ax.yaxis.set_major_locator(loc_majeurs_y)

plt.xlim(1, 1e5)
plt.ylim(5e-8, 4e-6)
plt.legend(loc='upper right', fontsize=11)
plt.grid(True, which="both", alpha=0.2)

plt.show()

# End
