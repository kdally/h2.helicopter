import numpy as np
import matplotlib.pyplot as plt


# Data
MTOW    = 4000 # kg
V_cr    = 350  # km/h
FM      = 0.75 # -
N_prop  = 2    # -
D       = 6    # m
rho     = 1.05 # kg/m3 at 5000 ft
eta_tr  = 0.96 # -
eta_ins = 0.98 # -
g       = 9.81 # kgm/s2

C_d     = 0.1
# S       =

def ideal_power(T, D, rho):
    return T*np.sqrt((2*T)/(rho*np.pi*D**2))

p_idle = ideal_power(MTOW*g, D, rho)
print(p_idle/1000/FM)


# P_ind   =
# P_prof  =
# P_para  =

