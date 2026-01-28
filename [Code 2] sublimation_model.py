# COUTANT Aubin L3 Internship: Modeling physical properties of Enceladus' icy surface
# Figure 7 Reproduction: Sublimation Validation
# Model: Gundlach et al. (2018)

# Start

import numpy as np
import matplotlib.pyplot as plt


# ==========================================
# 1. NUMERICAL SOLVER 
# ==========================================

def runge_kutta_4(f, y0, a, b, N):
    """
    Solves the differential equation y'(t) = f(t, y) using the RK4 method.
    
    NOTE: This numerical method was studied during the L3 Numerical Analysis course.
    It provides a 4th-order accuracy, which is superior to the Euler method 
    for long-duration simulations like this one.

    """
    h = (b - a) / N                                      # Time step [s]
    T = np.linspace(a, b, N + 1)
    Y = np.zeros(N + 1)
    Y[0] = y0                                            # Initial condition
    
    for i in range(N):
        # Calculate the 4 intermediate slopes (RK4 algorithm)
        k1 = h * f(T[i], Y[i])
        k2 = h * f(T[i] + h/2, Y[i] + 0.5 * k1)
        k3 = h * f(T[i] + h/2, Y[i] + 0.5 * k2)
        k4 = h * f(T[i] + h, Y[i] + k3)
        
        # Weighted average to advance to the next step
        Y[i+1] = Y[i] + (1/6) * (k1 + 2*k2 + 2*k3 + k4)
        
    return T, Y

# ==========================================
# 2. PARAMETERS AND CONSTANTS 
# ==========================================

R = 8.31447       # Ideal gas constant [J/(mol·K)]
mu = 18.015e-3    # Molar mass of water [kg/mol]
gamma = 0.2       # Specific surface energy [J/m^2]
T_gas = 300       # Chamber wall temperature [K]

# Residual pressure in the vacuum chamber (adjusted to experimental data)
p_gas = 7.65e-5   

# Empirical coefficient to match Gundlach's experimental data (Fitting parameter)
Scaling_factor = 0.17 

def p_sat(T):
    """
    Computes the saturation vapor pressure of water ice (Andreas 2007).
    """
    return np.exp(9.550426 - 5723.265/T + 3.53068*np.log(T) - 0.00728332*T)

def rho_ice(T):
    """
    Computes ice density as a function of temperature.
    """
    return 916 - 0.175 * (T-273.15) - 5.0e-4 * (T-273.15)**2

def hertz_knudsen(T_surface, p_gas_val, T_gas_val):
    """ 
    Calculates the net MASS flux Z (Eq. 8 in Gundlach et al., 2018).
    
    The term sqrt(mu / (2*pi*R*T)) converts the pressure into a mass flux [kg m^-2 s^-1].
    """
    # Outgoing flux (Sublimation) 
    flux_out = p_sat(T_surface) * np.sqrt(mu / (2 * np.pi * R * T_surface))
    
    # Incoming flux (Condensation) -> Depends on Gas Back-pressure
    flux_in = p_gas_val * np.sqrt(mu / (2 * np.pi * R * T_gas_val))
    
    return (flux_out - flux_in)

# ==========================================
# 3. PHYSICAL MODEL 
# ==========================================

def deriv_fonc(T_pa):
    """
    Prepares the derivative function dr/dt for a given temperature (T_pa).
    Implements Equation 7 from Gundlach et al. (2018).
    """
    rho = rho_ice(T_pa)
    
    # Theoretical mass flux Z (Eq. 8) scaled by the empirical factor
    Z_particule = Scaling_factor * hertz_knudsen(T_pa, p_gas, T_gas)
    
    # Surface recession velocity [m/s] (dr/dt without curvature)
    v_dim = Z_particule / rho
    
    def f(t, r):
        # Critical radius where curvature doubles the sublimation rate
        rc_val = 2 * mu * gamma / (rho * R * T_pa)
        
        # Curvature correction factor (Eq. 7 term in brackets)
        facteur_courbure = Y0 / (Y0 - rc_val)
        
        return - v_dim * facteur_courbure
    return f

# ==========================================
# 4. SIMULATION
# ==========================================

# Initial Conditions (Based on Gundlach Fig. 7)
Y0 = 1.5e-6        # Initial particle radius [m] 
t_max = 1e6        # Simulation duration [s] 
N_pas = 200000     # Number of time steps

# Simulation A (159.6 K)
f_A = deriv_fonc(159.6)
Temps_A, Rayon_A = runge_kutta_4(f_A, Y0, 0, t_max, N_pas)
Masse_Norm_A = (Rayon_A / Y0)**3

# Simulation B (163.4 K)
f_B = deriv_fonc(163.4)
Temps_B, Rayon_B = runge_kutta_4(f_B, Y0, 0, t_max, N_pas)
Masse_Norm_B = (Rayon_B / Y0)**3


# ==========================================
# 5. PLOTTING
# ==========================================

plt.figure(figsize=(8, 5))

# Plotting the results
plt.semilogx(Temps_A, Masse_Norm_A, label='Model T = 159.6 K', linestyle='--', linewidth=2)
plt.semilogx(Temps_B, Masse_Norm_B, label='Model T = 163.4 K', linestyle='--', linewidth=2)

# Graph styling
plt.xlabel('Time [s]', fontsize=12, fontweight='bold')
plt.ylabel('Normalized Mass (m(t)/m(0))', fontsize=12, fontweight='bold')
plt.title("Validation: Sublimation Model (Gundlach 2018, Fig. 7)", fontsize=13)
plt.text(4e2, 0.81, "$r_p$ = 1.5 µm", ha='center', fontsize=14, color='black', bbox=dict(facecolor='white', alpha=0.8))

plt.legend(fontsize=11)
plt.grid(True, which="both", linestyle=':', alpha=0.7)
plt.xlim(1e2, 1e5)
plt.ylim(0.2, 1.1)

plt.tight_layout()
plt.show()


#End
