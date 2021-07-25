## Sound absorption by viscothermal losses in square ducts
## Dr. Iván Herrero Durá
## July 25, 2021

import numpy as np
import matplotlib.pyplot as plt

# Definition of the frequency range

freq_min = 1 # Lower cutoff frequency of the duct [Hz]
freq_max = 5000 # Higher cutoff frequency of the duct [Hz]
freq = np.linspace(freq_min,freq_max,num=(freq_max-freq_min)+1) # Frequency [Hz]
w = 2*np.pi*freq # Angular frequency [rad^-1]


# Definition of the ambient parameters

p0 = 101320 # Atmospheric pressure [Pa]
mu = 1.568e-5 # Viscosity coefficient [kg/(m·s)]
gamma = 1.4 # Heat capacity ratio [1]
rho0 = 1.204 # Density [kg/m^3]
kappa0 = gamma*p0 # Bulk modulus of air [N/m^2]
c0 = np.sqrt(kappa0/rho0) # Wave speed [m/s]
Pr = 0.702 # Prandtl number at atmospheric pressure [1]
delta = np.sqrt((2*mu)/(rho0*w)) # Viscous boundary layer thickness [m]


# Definition of the geometrical parameters of the duct

d = 0.165 # Diameter of the duct [m]
L = 0.6 # Length of the duct [m]
A = np.pi*(d/2)**2 # Cross-section area [m]
Pt = np.pi*d # Wetted perimeter [m]
Rh = A/Pt # Hydraulic radius [m]
chi = np.sqrt(Pr)
kappa = (1-1j)/np.sqrt(2)
si = Rh/delta


# Definition of the viscothermal losses

k0 = (w/c0)*(1+(kappa/si)*(1+(gamma-1)/chi)) # Wavenumber [m^-1]
Z0 = ((rho0*c0)/A)*(1+(kappa/si)*(1-(gamma-1)/chi)) # Acoustic impedance of air [Pa*s/m^3]


# Definition of the Transfer Matrix

Zp = -1j*Z0/np.tan(k0*L) # Acoustic impedance of the rigid wall [Pa*s/m^3]
R = (Zp-Z0)/(Zp+Z0) # Reflection coefficient [1]
alpha = 1-abs(R**2) # Absorption coefficient [1]


# Plots

plt.plot(freq,alpha, 'b', label = r"$\alpha$")
plt.plot(freq,abs(R**2), 'r--', label = "$|R|^2$")
plt.axis([freq_min, freq_max, 0, 1])
plt.grid()
plt.legend()
plt.xlabel("f [Hz]")
plt.ylabel("$|R|^2$, "r"$\alpha$")