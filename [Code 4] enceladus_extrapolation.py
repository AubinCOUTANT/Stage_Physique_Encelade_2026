# COUTANT Aubin L3 Internship: Modeling physical properties of Enceladus' icy surface
# Extrapolation of the sintering model to astronomical time scales.
# Application to the specific thermal zoning of Enceladus.

# Start

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

# ==============================================================================
# 1. NUMERICAL SOLVER VARIABLE STEP RK4
# ==============================================================================

def runge_kutta_4_dynamic(f, y0, T_array):
    """ 
    Vectorized Runge-Kutta 4 solver adapted for logarithmic time scales.
    Allows covering a range from milliseconds to the age of the Solar System.
    """
    N = len(T_array) - 1
    y0 = np.asarray(y0)
    Y = np.zeros((N + 1, len(y0)))
    Y[0] = y0
    
    for i in range(N):
        t_curr = T_array[i]
        # The time step h is calculated dynamically based on the array spacing
        h = T_array[i+1] - T_array[i]
        y_curr = Y[i]
        
        k1 = h * f(t_curr, y_curr)
        k2 = h * f(t_curr + h/2, y_curr + 0.5 * k1)
        k3 = h * f(t_curr + h/2, y_curr + 0.5 * k2)
        k4 = h * f(t_curr + h, y_curr + k3)
        
        # State update
        Y[i+1] = y_curr + (1/6)*(k1 + 2*k2 + 2*k3 + k4)
        
        # Numerical safety
        if Y[i+1][0] < 1e-10: Y[i+1][0] = 1e-10
        if Y[i+1][1] < 1e-10: Y[i+1][1] = 1e-10
              
    return Y

# ==============================================================================
# 2. PHYSICAL PARAMETERS AND CONSTANTS
# ==============================================================================

# Universal Constants
R       = 8.31447       # Ideal gas constant [J/(mol·K)]

# Water Ice Properties
mu      = 18.015e-3     # Molar mass [kg/mol]
gamma   = 0.2           # Surface energy [J/m^2]
omega   = 2.0e-5        # Molar volume [m^3/mol]
alpha   = np.pi/2       # Packing factor (Porosity ~0.5)

# Scaling Factors (Model Calibration)
SCALING_PARTICULE = 0.9
SCALING_NECK      = 0.01

# Conversion Factor: Seconds to Years
SEC_TO_YEAR = 1 / (365.25 * 24 * 3600)

# ==============================================================================
# 3. THERMODYNAMIC FUNCTIONS 
# ==============================================================================

def rho_ice(T):
    """ Ice density as a function of T [kg/m^3] """
    return 916 - 0.175 * T - 5.0e-4 * T**2

def p_sat(T):
    """ 
    Saturation vapor pressure according to the Andreas (2007) model.
    """
    return np.exp(9.550426 - 5723.265/T + 3.53068*np.log(T) - 0.00728332*T)

def hertz_knudsen_flux(T, P_sat):
    """ 
    Theoretical maximum mass flux Z [kg/(m^2 s)] (Hertz-Knudsen Law).
    """
    return P_sat * np.sqrt(mu / (2 * np.pi * R * T))

# ==============================================================================
# 4. MODEL CORE (DIFFERENTIAL SYSTEM)
# ==============================================================================

def system_derivatives(T_env):
    """
    Generates the derivative function f(t, Y) for a given temperature.
    """
    # Calculation of static thermodynamic properties
    rho = rho_ice(T_env)
    P_sat = p_sat(T_env)
    
    # Theoretical flux and sintering efficiency
    Z_theory = hertz_knudsen_flux(T_env, P_sat)
    zeta = (omega**2 * gamma * P_sat) / (R * T_env * np.sqrt(2 * np.pi * mu * R * T_env))
    
    def f(t, Y):
        rp = Y[0]
        rn = Y[1]
        
        # 1. Particle Evolution 
        drp_dt = - (Z_theory / rho) * SCALING_PARTICULE
        
        # 2. Neck Evolution 
        
        # Geometric factors (Gundlach Eqs. 4-6)
        denom = 2 * (rp - rn)
        if denom <= 1e-20: denom = 1e-20
        
        Delta = rn**2 / denom
        
        # Calculation of angle phi
        ratio = rn/rp
        if ratio > 1.0: ratio = 1.0
        phi = np.arcsin(ratio)
        
        # Source/sink distance terms
        d_source = rp * (alpha/2 + phi - np.pi/2)
        d_sink   = Delta * np.arctan(rp / (rn + Delta))
        S = d_source / (d_source - d_sink)
        
        # Curvature term 
        curvature_term = (2/rp) + (1/Delta) - (1/rn)
        
        # A. Constructive Term (Sintering)
        growth_rate = zeta * S * curvature_term
        
        # B. Destructive Term (Neck Sublimation)
        decay_rate = (Z_theory / rho) * SCALING_NECK
        
        drn_dt = growth_rate - decay_rate
        
        # 3. Physical Constraint
        if rn >= rp: 
            rn = rp
            drn_dt = drp_dt 
            
        return np.array([drp_dt, drn_dt])
    
    return f

# ==============================================================================
# 5. FIGURE INITIALIZATION
# ==============================================================================

plt.figure(figsize=(12, 8)) 

# Initial Conditions (Standard 0.5 µm grains for Fig 12)
rp_init = 0.5e-6        
rn_init = 3.6e-8        # Initial adhesion 
Y0 = [rp_init, rn_init]

y_txt = 1.45e-7         # Vertical position for temperature labels

# ==============================================================================
# 6. SIMULATIONS 
# ==============================================================================

# Note: All "Temps" arrays are in seconds for the solver, 
# but plotted as "Temps * SEC_TO_YEAR"

# 1: PLAINS (80 K)
T = 80.0
Temps = np.geomspace(1e14, 5e19, 5000)
res = runge_kutta_4_dynamic(system_derivatives(T), Y0, Temps)
plt.loglog(Temps * SEC_TO_YEAR, res[:, 1], color='#00BFFF', linewidth=3.5) 
plt.text(2e19 * SEC_TO_YEAR, y_txt, "80 K", fontsize=11, fontweight='bold', color='black', ha='center')

# 2: TIGER STRIPES MARGINS (100 K - 160 K)

# 100 K
T = 100.0
Temps = np.geomspace(1e6, 1.8e13, 5000)
res = runge_kutta_4_dynamic(system_derivatives(T), Y0, Temps)
plt.loglog(Temps * SEC_TO_YEAR, res[:, 1], color='#1E90FF', linewidth=3)
plt.text(5e12 * SEC_TO_YEAR, y_txt, "100 K", fontsize=11, fontweight='bold', color='black', ha='center')

# 120 K
T = 120.0
Temps = np.geomspace(1e3, 9e8, 5000)  
res = runge_kutta_4_dynamic(system_derivatives(T), Y0, Temps)
plt.loglog(Temps * SEC_TO_YEAR, res[:, 1], color='#4169E1', linewidth=3)
plt.text(4e8 * SEC_TO_YEAR, y_txt, "120 K", fontsize=11, fontweight='bold', color='black', ha='center')

# 140 K 
T = 140.0
Temps = np.geomspace(1, 7e5, 5000)
res = runge_kutta_4_dynamic(system_derivatives(T), Y0, Temps)
plt.loglog(Temps * SEC_TO_YEAR, res[:, 1], color='#0000FF', linewidth=3)
plt.text(6e5 * SEC_TO_YEAR, y_txt, "140 K", fontsize=11, fontweight='bold', ha='center', color='black')

# 160 K 
T = 160.0
Temps = np.geomspace(1e-2, 3.3e3, 5000)
res = runge_kutta_4_dynamic(system_derivatives(T), Y0, Temps)
plt.loglog(Temps * SEC_TO_YEAR, res[:, 1], color='#0000CC', linewidth=3)
plt.text(4e3 * SEC_TO_YEAR, y_txt, "160 K", fontsize=11, fontweight='bold', ha='center', color='black')

# 3: HOT SPOT (180 K)
T = 180.0
Temps = np.geomspace(1e-3, 5.1e1, 5000)
res = runge_kutta_4_dynamic(system_derivatives(T), Y0, Temps)
plt.loglog(Temps * SEC_TO_YEAR, res[:, 1], color='#000099', linewidth=3)
plt.text(6e1 * SEC_TO_YEAR, y_txt, "180 K", fontsize=11, fontweight='bold', ha='center', color='black')

# 4: THEORETICAL LIMITS (200 K - 220 K) 

# 200 K 
T = 200.0
Temps = np.geomspace(1e-3, 1.8, 5000)
res = runge_kutta_4_dynamic(system_derivatives(T), Y0, Temps)
plt.loglog(Temps * SEC_TO_YEAR, res[:, 1], color='#000066', linewidth=3)
plt.text(1e0 * SEC_TO_YEAR, y_txt, "200 K", fontsize=11, fontweight='bold', ha='center', color='black')

# 220 K 
T = 220.0
Temps = np.geomspace(1e-3, 1.1e-1, 5000)
res = runge_kutta_4_dynamic(system_derivatives(T), Y0, Temps)
plt.loglog(Temps * SEC_TO_YEAR, res[:, 1], color='#000033', linewidth=3)
plt.text(1.5e-2 * SEC_TO_YEAR, y_txt, "220 K", fontsize=11, fontweight='bold', ha='center', color='black')


# ==============================================================================
# 7. ANNOTATIONS (Positions converted to Years)
# ==============================================================================

y_bar = 1.75e-7
y_petit = 1.6e-7

# Annotation: Enceladus Plains
x_text = 2e19 * SEC_TO_YEAR
plt.text(x_text, 2e-7, "Plains in \n South Pole", fontsize=14, fontweight='bold', color='black', ha='center')
plt.plot([x_text, x_text], [1.6e-7, 1.9e-7], color='black', linewidth=1.5)

# Annotation: Tiger Stripes Margins (100-160 K)
x_start_TS = 4e3 * SEC_TO_YEAR    # 160K Position
x_end_TS = 2e13 * SEC_TO_YEAR     # 100K Position
plt.plot([x_start_TS, x_end_TS], [y_bar, y_bar], color='black', linewidth=1.5)       
plt.plot([x_start_TS, x_start_TS], [y_bar, y_petit], color='black', linewidth=1.5)    
plt.plot([x_end_TS, x_end_TS], [y_bar, y_petit], color='black', linewidth=1.5)        
plt.text(3e8 * SEC_TO_YEAR, 2e-7, "Tiger Stripes\nMargins", fontsize=14, fontweight='bold', color='black', ha='center')
plt.plot([3e8 * SEC_TO_YEAR, 3e8 * SEC_TO_YEAR], [1.9e-7, y_bar], color='black', linewidth=1.5) 

# Annotation: Hot Spot (180 K)
x_text_hot = 6e1 * SEC_TO_YEAR
plt.plot([x_text_hot, x_text_hot], [1.6e-7, 1.9e-7], color='black', linewidth=1.5)
plt.text(x_text_hot, 2e-7, "Hottest\nPoint", fontsize=14, fontweight='bold', color='black', ha='center')

# Annotation: Theoretical Limits (200-220 K)
x_start_Lim = 1.5e-2 * SEC_TO_YEAR # 220K Position
x_end_Lim = 1e0 * SEC_TO_YEAR      # 200K Position
plt.plot([x_start_Lim, x_end_Lim], [y_bar, y_bar], color='black', linewidth=1.5)
plt.plot([x_start_Lim, x_start_Lim], [y_bar, y_petit], color='black', linewidth=1.5)
plt.plot([x_end_Lim, x_end_Lim], [y_bar, y_petit], color='black', linewidth=1.5)
plt.text(1e-1 * SEC_TO_YEAR, 2e-7, "Theoretical\nLimits", fontsize=14, fontweight='bold', color='black', ha='center')
plt.plot([1e-1 * SEC_TO_YEAR, 1e-1 * SEC_TO_YEAR], [1.9e-7, y_bar], color='black', linewidth=1.5)


# ==============================================================================
# 8. FINAL FORMATTING
# ==============================================================================

# Title
plt.title("Surface Consolidation on Enceladus: Thermal Zoning Analysis ( $r_p = 0.5$ µm )", fontsize=16)
plt.xlabel("Time [Years]", fontsize=14, fontweight='bold') # LABEL CHANGÉ
plt.ylabel("Sinter Neck Radius [m]", fontsize=14, fontweight='bold')

# Limit lines (Solar System Age / Initial Radius)
# 4.5e9 years instead of 1.45e17 seconds
plt.axvline(x=4.5e9, color='black', linestyle='--', linewidth=2)
plt.axhline(y=rn_init, color='gray', linestyle=':', linewidth=2 )

plt.text(1.5e9, 1.3e-7, "Age of the Solar System", rotation=90, fontsize=12, va='bottom')
plt.text(1e-10, 3.3e-8, "Initial Neck Radius (Adhesion)", ha='left', fontsize=11) # Adjusted X position for years

# Axis parameters (Log-Log)
ax = plt.gca()
ax.tick_params(which='both', direction='in', top=True, right=True, labelsize=12)
ax.tick_params(which='major', length=8, width=1.5)
ax.tick_params(which='minor', length=4, width=1)
loc_majeurs_y = ticker.LogLocator(base=10.0, numticks=10)
ax.yaxis.set_major_locator(loc_majeurs_y)

# Limits ajusted for Years (1e-3 sec becomes ~3e-11 years)
plt.xlim(1e-11, 1e13) 
plt.ylim(3e-8, 3e-7)
plt.grid(True, which="both", alpha=0.4, color='grey')

plt.tight_layout()
plt.show()

# End
