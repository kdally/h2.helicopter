from math import *
import numpy as np
import matplotlib.pyplot as plt


def density(height):
    temp = 288.15-0.0065*height
    pressure = 101325*pow(temp/288.15,(-9.80665/(-0.0065*287.05)))
    return pressure/(287.05*temp)

MTOM = 4000     # kg
g    = 9.80665  # m/s2
MTOW = MTOM * g # N
C_l  = 0.6      # -
C_d0  = 0.02    # -
C_dr = 0.3      # -
R    = 2        # m
S    = 12       # m2
wingchord = 1.5 # m
wingspan  = S/wingchord # m
A    = wingspan/wingchord # -
e    = 0.8      # -
FM   = 0.65     # -
h    = 5000     # m
v_cr = 111.1111 # m/s
RoC  = 7        # m/s

def lift(MTOW, rho, V, S, C_lmax):
    L = 0.5 * rho * V**2 * S * C_lmax
    return  L if L<MTOW else MTOW

def drag(rho, V, S, C_d0, C_l):
    C_d = C_d0 + (C_l**2)/(pi*A*e)
    return  0.5 * rho * V**2 * S * C_d

def p_parasite(rho, V, S, C_d, C_l):
    D = drag(rho, V, S, C_d, C_l)
    p_parasite = D * V
    return p_parasite/1000

def p_induced(MTOW,rho,V,S,C_l,wingchord,R,C_dr,FM):
    L_rotor = MTOW - 0.5 * rho * V**2 * S * C_l
    A_p = 2 * wingchord * R
    A = 2 * pi * R**2
    Thrust = (1 + (C_dr * A_p)/A) * MTOW
    V_induced = np.sqrt(-((V**2)/2) + np.sqrt((V**4)/4 + (Thrust/(2 * rho * A))**2))
    p_induced = L_rotor * V_induced / FM
    return np.where(p_induced>0, p_induced, 0)/1000

def p_climb(MTOW, RoC, P_r):
    return MTOW*RoC/1000

V = np.arange(0,112)



parasite = p_parasite(density(h), V, S, C_d0, C_l)
induced = p_induced(MTOW, density(h), V, S, C_l, wingchord, R, C_dr, FM)
total_level = parasite + induced
climb = p_climb(MTOW, RoC, total_level)

total = total_level+climb





plt.plot(V, parasite, label= 'parasite')
plt.plot(V, induced, label= 'induced')
plt.plot(total, label= 'total')
plt.legend()
plt.show()
