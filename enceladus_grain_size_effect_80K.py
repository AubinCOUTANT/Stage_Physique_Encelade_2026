# COUTANT Aubin L3 Internship: Modeling physical properties of Enceladus' icy surface
# Study of the influence of grain size (granulometry) on sintering kinetics.
# Includes dynamic time-stepping per grain size to avoid numerical artifacts.

# Start

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

# ==============================================================================
# 1. NUMERICAL SOLVER
# ==============================================================================

def runge_kutta_4_vec(f, y0, T_array):
    """ Vectorized Runge-Kutta 4 solver """
    N = len(T_array) - 1
    y0 = np.asarray(y0)
    Y = np.zeros((N + 1,) + y0.shape)
    Y[0] = y0
    
    for i in range(N):
        t_curr = T_array[i]
        h = T_array[i+1] - T_array[i]
        y_curr = Y[i]
        
        k1 = h * f(t_curr, y_curr)
        k2 = h * f(t_curr + h/2, y_curr + 0.5 * k1)
        k3 = h * f(t_curr + h/2, y_curr + 0.5 * k2)
        k4 = h * f(t_curr + h, y_curr + k3)
        
        Y[i+1] = y_curr + (1/6) * (k1 + 2*k2 + 2*k3 + k4)
        
        # Numerical safety
        if Y[i+1][0] < 1e-10: Y[i+1][0] = 1e-10
        if Y[i+1][1] < 1e-10: Y[i+1][1] = 1e-10
        
        # Check for full solidification 
        # If neck >= particle, we can theoretically stop or clamp
        if Y[i+1][1] >= 0.999 * Y[i+1][0]:
            Y[i+1][1] = Y[i+1][0]
            
    return Y

# ==============================================================================
# 2. PHYSICAL PARAMETERS
# ==============================================================================

# Universal Constants
R = 8.31447       # Ideal gas constant [J/(mol·K)]

# Water Ice Properties
mu      = 18.015e-3     # Molar mass [kg/mol]
gamma   = 0.2           # Surface energy [J/m^2]
omega   = 2.0e-5        # Molar volume [m^3/mol]
alpha   = np.pi/2       # Packing factor (Porosity ~0.5)

# JKR Elastic Properties
E = 10.5e9         
nu = 0.31          
K = (2/3) * E / (1 - nu**2)

# Scaling Factors
SCALING_PARTICULE = 0.9
SCALING_NECK      = 0.01

# ==============================================================================
# 3. PHYSICS & THERMODYNAMICS
# ==============================================================================

def rho_ice(T):
    """ 
    Calculates the density of water ice [kg/m^3] as a function of Temperature.
    Based on standard empirical fits for low-temperature ice evolution.
    """
    return 916 - 0.175 * T - 5.0e-4 * T**2

def p_sat(T):
    """ 
    Computes the Saturation Vapor Pressure of ice [Pa].
    Source: Andreas (2007), derived from Murphy & Koop (2005)..
    """
    return np.exp(9.550426 - 5723.265/T + 3.53068*np.log(T) - 0.00728332*T)

def JKR_theory(rp):
    """ 
    Calculates the initial neck radius (r_n0) formed immediately upon contact due to adhesion.
    Based on the JKR (Johnson-Kendall-Roberts, 1971) theory for elastic contact.
    This defines the starting point of cohesion before any sintering occurs.
    """
    return (6 * np.pi * gamma * (rp/2)**2 / K)**(1/3)

def hertz_knudsen_flux(T, P_sat):
    """ 
    Calculates the theoretical maximum mass flux Z [kg/(m^2 s)].
    Based on the Hertz-Knudsen-Langmuir equation.
    It quantifies the kinetic exchange of water molecules between the solid surface and the vacuum.
    """
    return P_sat * np.sqrt(mu / (2 * np.pi * R * T))

# ==============================================================================
# 4. DIFFERENTIAL MODEL
# ==============================================================================

def system_derivatives(T_env):
    rho = rho_ice(T_env)
    P_sat = p_sat(T_env)
    Z_theory = hertz_knudsen_flux(T_env, P_sat)
    zeta = (omega**2 * gamma * P_sat) / (R * T_env * np.sqrt(2 * np.pi * mu * R * T_env))
    
    def f(t, Y):
        rp = Y[0]
        rn = Y[1]
        
        drp_dt = - (Z_theory / rho) * SCALING_PARTICULE
        
        denom = 2 * (rp - rn)
        if denom <= 1e-20: denom = 1e-20
        Delta = rn**2 / denom
        
        ratio = rn/rp
        if ratio > 1.0: ratio = 1.0
        phi = np.arcsin(ratio)
        
        d_source = rp * (alpha/2 + phi - np.pi/2)
        d_sink   = Delta * np.arctan(rp / (rn + Delta))
        S = d_source / (d_source - d_sink)
        
        curvature_term = (2/rp) + (1/Delta) - (1/rn)
        
        growth_rate = zeta * S * curvature_term
        decay_rate = (Z_theory / rho) * SCALING_NECK
        
        drn_dt = growth_rate - decay_rate
        
        if rn >= rp: 
            rn = rp
            drn_dt = drp_dt 
            
        return np.array([drp_dt, drn_dt])
    return f

# ==============================================================================
# 5. SIMULATION SETUP (GRAIN SIZE STUDY)
# ==============================================================================

T_study = 80.0  # Temperature (Enceladus Plains)

# Grain sizes to test
grain_sizes = [0.5e-6, 5e-6, 20e-6, 50e-6, 100e-6]
colors = ['#00BFFF', '#1E90FF', '#4169E1', '#0000FF', '#000080']
labels = ['0.5 µm', '5 µm', '20 µm', '50 µm', '100 µm']

plt.figure(figsize=(12, 8))
func_system = system_derivatives(T_study)

# Loop over each grain size
for i, rp_val in enumerate(grain_sizes):
    
    # 1. Define SPECIFIC time array for each grain
    t_start = 1e13
    t_max   = 1e23 
    
    Time_array = np.geomspace(t_start, t_max, 5000)
    
    # 2. Initial Conditions (JKR)
    rn_val = JKR_theory(rp_val)
    Y0 = [rp_val, rn_val]
    
    # 3. Run Simulation
    res = runge_kutta_4_vec(func_system, Y0, Time_array)
    
    # 4. Clean Data (Stop plotting if solidified)
    neck_radius = res[:, 1]
    particle_radius = res[:, 0]
    
    valid_mask = neck_radius < (0.999 * particle_radius)
    
    # Plotting only valid part
    plt.loglog(Time_array[valid_mask], neck_radius[valid_mask], 
               color=colors[i], linewidth=2.5, label=f'$r_p$ = {labels[i]}')
    
    # JKR Level line
    plt.axhline(y=rn_val, color=colors[i], linestyle=':', alpha=0.6, linewidth=1)

# ==============================================================================
# 6. FORMATTING
# ==============================================================================

plt.title(f"Influence of Grain Size on Sintering at T = {T_study} K (Enceladus Plains)", fontsize=16)
plt.xlabel("Time [s]", fontsize=14, fontweight='bold')
plt.ylabel("Sinter Neck Radius ($r_n$) [m]", fontsize=14, fontweight='bold')

plt.axvline(x=1.45e17, color='black', linestyle='--', linewidth=2)
plt.text(5e16, 1.2e-6, "Age of the Solar System", rotation=90, fontsize=12, va='bottom')

plt.legend(title="Initial Grain Size", loc='upper left', fontsize=12)
plt.grid(True, which="both", alpha=0.4, color='grey')

ax = plt.gca()
ax.tick_params(which='both', direction='in', top=True, right=True, labelsize=12)
ax.tick_params(which='major', length=8, width=1.5)
ax.tick_params(which='minor', length=4, width=1)
loc_majeurs_y = ticker.LogLocator(base=10.0, numticks=10)
ax.yaxis.set_major_locator(loc_majeurs_y)

plt.tight_layout()
plt.savefig('Fig_Granulometry_80K.png', dpi=300)
plt.show()

# End
