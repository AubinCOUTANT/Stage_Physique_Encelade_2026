# COUTANT Aubin L3 Internship: Modeling physical properties of Enceladus' icy surface
# Extrapolation of the sintering model (validated in Fig. 8) to astronomical
# time scales for different Solar System bodies (Comets, Europa, Enceladus).
# This script reproduces Figure 12 of Gundlach et al. (2018).

# Start

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

# ==============================================================================
# 1. NUMERICAL SOLVER (VARIABLE STEP RK4)
# ==============================================================================

def runge_kutta_4_dynamic(f, y0, T_array):
    """ 
    Vectorized Runge-Kutta 4 solver adapted for logarithmic time scales.
    Allows covering a range from milliseconds (cometary dynamics) 
    to the age of the Solar System (icy moons).
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
        k2 = h * f(t_curr + h/2, y_curr + 0.5*k1)
        k3 = h * f(t_curr + h/2, y_curr + 0.5*k2)
        k4 = h * f(t_curr + h, y_curr + k3)
        
        # State update
        Y[i+1] = y_curr + (1/6)*(k1 + 2*k2 + 2*k3 + k4)
        
        # Numerical safety: avoid negative or zero radii which would cause divergence
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
# ------------------------------------------
# These factors adapt the theoretical flux to experimental reality (roughness, geometry).
# - SCALING_PARTICULE (0.9) : Surface accommodation coefficient (Andreas 2007).
# - SCALING_NECK (0.01) : Geometric confinement factor. In the neck cavity,
#   99% of sublimated molecules re-condense immediately.
SCALING_PARTICULE = 0.9
SCALING_NECK      = 0.01

# ==============================================================================
# 3. THERMODYNAMIC FUNCTIONS (ANDREAS LAW)
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
    Assumption: Spatial vacuum (P_gas ~ 0), no return flux.
    """
    return P_sat * np.sqrt(mu / (2 * np.pi * R * T))

# ==============================================================================
# 4. MODEL CORE (DIFFERENTIAL SYSTEM)
# ==============================================================================

def system_derivatives(T_env):
    """
    Generates the derivative function f(t, Y) for a given temperature.
    Y[0] = Particle radius (rp)
    Y[1] = Neck radius (rn)
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
        
        # Angle phi calculation
        ratio = rn/rp
        if ratio > 1.0: ratio = 1.0
        phi = np.arcsin(ratio)
        
        # Source/sink distance terms
        d_source = rp * (alpha/2 + phi - np.pi/2)
        d_sink   = Delta * np.arctan(rp / (rn + Delta))
        
        # Geometric transport factor S 
        S = d_source / (d_source - d_sink)
        
        # Curvature term 
        curvature_term = (2/rp) + (1/Delta) - (1/rn)
        
        # A. Constructive Term (Sintering)
        growth_rate = zeta * S * curvature_term
        
        # B. Destructive Term (Neck Sublimation)
        # Application of the geometric confinement factor (0.01)
        decay_rate = (Z_theory / rho) * SCALING_NECK
        
        drn_dt = growth_rate - decay_rate
        
        # 3. Physical Constraint
        # If the neck joins the particle, they evolve together (Solidification)
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
rn_init = 3.6e-8       # Initial adhesion 
Y0 = [rp_init, rn_init]

y_txt = 1.45e-7        # Vertical position for temperature labels

# ==============================================================================
# 6. SIMULATION: ENCELADUS (Inert Regime)
# ==============================================================================

# Simulation at 80 K over the age of the solar system
T = 80.0
Temps = np.geomspace(1e14, 5e19, 5000)
res = runge_kutta_4_dynamic(system_derivatives(T), Y0, Temps)

# Plot Enceladus
plt.loglog(Temps, res[:, 1], color='#00BFFF', linewidth=3.5) 
plt.text(2e19, y_txt, "80 K", fontsize=11, fontweight='bold', color='black', ha='center')
plt.text(2e19, 2e-7, "Enceladus\n(plains)", fontsize=16, fontweight='bold', color='black', ha='center')
plt.plot([2e19, 2e19], [1.6e-7, 1.9e-7], color='black', linewidth=1.5)

# ==============================================================================
# 7. SIMULATION: EUROPA (Intermediate Regime)
# ==============================================================================

# 100 K
T = 100.0
Temps = np.geomspace(1e6, 1.8e13, 5000)
res = runge_kutta_4_dynamic(system_derivatives(T), Y0, Temps)
plt.loglog(Temps, res[:, 1], color='#1E90FF', linewidth=3)
plt.text(5e12, y_txt, "100 K", fontsize=11, fontweight='bold', color='black', ha='center')

# 120 K
T = 120.0
Temps = np.geomspace(1e3, 9e8, 5000) 
res = runge_kutta_4_dynamic(system_derivatives(T), Y0, Temps)
plt.loglog(Temps, res[:, 1], color='#4169E1', linewidth=3)
plt.text(4e8, y_txt, "120 K", fontsize=11, fontweight='bold', color='black', ha='center')

# Europa Annotation
x_start, x_end = 9e7, 2e13     
y_bar, y_petit = 1.75e-7, 1.6e-7   
plt.plot([x_start, x_end], [y_bar, y_bar], color='black', linewidth=1.5)      
plt.plot([x_start, x_start], [y_bar, y_petit], color='black', linewidth=1.5)    
plt.plot([x_end, x_end], [y_bar, y_petit], color='black', linewidth=1.5)        
plt.plot([4.5e10, 4.5e10], [y_bar, 1.9e-7], color='black', linewidth=1.5)
plt.text(4.5e10, 2e-7, "Europa", fontsize=16, fontweight='bold', color='black', ha='center')

# ==============================================================================
# 8. SIMULATION: COMETS (Active Regime)
# ==============================================================================

# 140 K 
T = 140.0
Temps = np.geomspace(1, 7e5, 5000)
res = runge_kutta_4_dynamic(system_derivatives(T), Y0, Temps)
plt.loglog(Temps, res[:, 1], color='#0000FF', linewidth=3)
plt.text(6e5, y_txt, "140 K", fontsize=11, fontweight='bold', ha='center', color='black')

# 160 K 
T = 160.0
Temps = np.geomspace(1e-2, 3.3e3, 5000)
res = runge_kutta_4_dynamic(system_derivatives(T), Y0, Temps)
plt.loglog(Temps, res[:, 1], color='#0000CC', linewidth=3)
plt.text(4e3, y_txt, "160 K", fontsize=11, fontweight='bold', ha='center', color='black')

# 180 K 
T = 180.0
Temps = np.geomspace(1e-3, 5.1e1, 5000)
res = runge_kutta_4_dynamic(system_derivatives(T), Y0, Temps)
plt.loglog(Temps, res[:, 1], color='#000099', linewidth=3)
plt.text(6e1, y_txt, "180 K", fontsize=11, fontweight='bold', ha='center', color='black')

# 200 K 
T = 200.0
Temps = np.geomspace(1e-3, 1.8, 5000)
res = runge_kutta_4_dynamic(system_derivatives(T), Y0, Temps)
plt.loglog(Temps, res[:, 1], color='#000066', linewidth=3)
plt.text(1e0, y_txt, "200 K", fontsize=11, fontweight='bold', ha='center', color='black')

# 220 K 
T = 220.0
Temps = np.geomspace(1e-3, 1.1e-1, 5000)
res = runge_kutta_4_dynamic(system_derivatives(T), Y0, Temps)
plt.loglog(Temps, res[:, 1], color='#000033', linewidth=3)
plt.text(1.5e-2, y_txt, "220 K", fontsize=11, fontweight='bold', ha='center', color='black')

# Comets Annotation
x_start, x_end = 4e-3, 3e4      
y_bar, y_petit = 1.75e-7, 1.6e-7
plt.plot([x_start, x_end], [y_bar, y_bar], color='black', linewidth=1.5)      
plt.plot([x_start, x_start], [y_bar, y_petit], color='black', linewidth=1.5)  
plt.plot([x_end, x_end], [y_bar, y_petit], color='black', linewidth=1.5)      
plt.plot([0.9e1, 0.9e1], [y_bar, 1.9e-7], color='black', linewidth=1.5)
plt.text(0.9e1, 2e-7, "Comets @ 1-3 AU", fontsize=16, fontweight='bold', color='black', ha='center')

# ==============================================================================
# 9. FINAL GRAPH FORMATTING
# ==============================================================================

plt.title("Reproduction des modèles théoriques (Fig. 12 Gundlach 2018)", fontsize=16)
plt.xlabel("Time [s]", fontsize=14, fontweight='bold')
plt.ylabel("Sinter Neck Radius [m]", fontsize=14, fontweight='bold')

# Limit Lines (Solar System Age / Initial Radius)
plt.axvline(x=1.45e17, color='black', linestyle='--', linewidth=2)
plt.axhline(y=rn_init, color='gray', linestyle=':', linewidth=2 )

plt.text(5e16, 1.3e-7, "Age of the Solar System", rotation=90, fontsize=12, va='bottom')
plt.text(1e-2, 3.3e-8, "Initial Neck Radius (Adhesion)", ha='left', fontsize=11)
plt.text(2e-2, 2.5e-7, "$r_p = 0.5$ µm", fontsize=22, color='black')

# Axis Parameters (Log-Log)
ax = plt.gca()
ax.tick_params(which='both', direction='in', top=True, right=True, labelsize=12)
ax.tick_params(which='major', length=8, width=1.5)
ax.tick_params(which='minor', length=4, width=1)
loc_majeurs_y = ticker.LogLocator(base=10.0, numticks=10)
ax.yaxis.set_major_locator(loc_majeurs_y)

plt.xlim(1e-3, 1e21)
plt.ylim(3e-8, 3e-7)
plt.grid(True, which="both", alpha=0.4, color='grey')

plt.tight_layout()
plt.savefig('Fig_12_reproduction_final.png', dpi=300)
plt.show()

# End
