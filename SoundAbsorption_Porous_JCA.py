## Sound absorption by porous materials (Johnson-Champoux-Allard model)
## Dr. Iván Herrero Durá
## July 24, 2021


import numpy as np
import matplotlib.pyplot as plt

# Definition of the frequency range

freq_min = 1 # Lower cutoff frequency of the duct [Hz]
freq_max = 10000 # Higher cutoff frequency of the duct [Hz]
freq = np.linspace(freq_min,freq_max,num=(freq_max-freq_min)+1) # Frequency [Hz]
w = 2*np.pi*freq # Angular frequency [rad^-1]

T = 293.15 # Temperature [K]
p = 101325 # Absolute pressure [atm]


# Fluid Properties

rho0 = 1.21 # Density [kg/m^3]
c0 = 343.2 # Speed of sound [m/s]
Z0 = rho0*c0 # Acoustic impedance of air
Cp = 1000 # Heat capacity at constant pressure [J/(kg·K)]
gamma = 1.401 # Ratio of specific heats [1]
kappa = 0.0257 # Thermal conductivity [W/(m·K)]
mu = 1.67581e-5 # Dynamic viscosity [Pa·s]


# Porous Matrix Properties

phi = 0.99 # Porosity [1]
Rf = 14000 # Flow resistivity [Pa·s/m^2]
Lv = 100e-6 # Viscous characteristic length [m]
Lth = 300e-6 # Thermal characteristic length [m]
tau = 1.02 # Tortuosity factor [1]
d = 0.051 # Length of the sample [m] 

Pr = (mu*Cp)/kappa # Prandtl number
rho_eff = ((tau*rho0)/phi)*(1+((Rf*phi)/(1j*w*rho0*tau))*np.sqrt(1+((4*1j*w*tau**2*mu*rho0)/(Rf**2*Lv**2*phi**2)))) # Effective density [kg/m^3]
k_eff = ((gamma*p)/phi)*(gamma-(gamma-1)*(1+((8*mu)/(1j*w*Lth**2*Pr*rho0))*np.sqrt(1+((1j*w*Lth**2*Pr*rho0)/(16*mu))))**(-1))**(-1) # Effective bulk modulus
c_eff = np.sqrt(k_eff/rho_eff) # Effective speed of sound [m/s]

k = w/c_eff

Zi = -1j*rho_eff*c_eff/np.tan(k*d)
R = (Zi-Z0)/(Zi+Z0)
alpha = 1-abs(R)**2


# Plots

plt.plot(freq,alpha,'b')
plt.axis([freq_min, freq_max, 0, 1])
plt.grid()
plt.xlabel("f [Hz]")
plt.ylabel(r"$\alpha$")