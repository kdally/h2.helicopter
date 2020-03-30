import math
import numpy as np
import matplotlib.pyplot as plt

MTOM = 4000
MTOW = MTOM*9.81
dl = 1.15
M = 1
velocityTAS = 0
N = 2
B = 3


DL = np.arange(20, 100, 1)

r = np.sqrt(MTOM/(DL*2*math.pi))
c = 0.2358
omega = 85

def density(height):
    temp = 288.15-0.0065*height
    pressure = 101325*math.pow(temp/288.15,(-9.80665/(-0.0065*287.05)))
    return pressure/(287.05*temp)

def vHover(thrust,rho):
    #velocity induced in hover
    return np.sqrt((dl*thrust)/(2*N*rho*np.pi*r**2))

def c_l_blade(L_rotor, rho):
    print(6*L_rotor / (N*B*rho*omega**2*c*r**3))
    return 6*L_rotor / (N*B*rho*omega**2*c*r**3)

def c_d_blade_1(x):
    return 20441 * x**6 -1505.9 * x**5 + 442.64 *x**4 - 63.716 * x**3 + 4.7524 * x**2 - 0.0961 * x + 0.0058

def p_profile(c_d_blade, rho):
    return np.where(c_d_blade != 0.0058, N * B * c_d_blade * rho * omega**3 * r**4  * c /8, 0)

def inducedPower(thrust,rho):
    #find induced power
    eta_trans = 1#0.96
    eta_inst = 0.98
    induced = (thrust*dl/M)*np.sqrt((-0.5*velocityTAS**2)+np.sqrt(((0.5*velocityTAS**2)**2)+vHover(thrust,rho)**4))*(1/(eta_trans*eta_inst))
    C_l_blade = c_l_blade(thrust,rho)
    aoaob = C_l_blade/(2*math.pi)
    C_d_blade = c_d_blade_1(aoaob)
    P_profile = p_profile(C_d_blade,rho)
    return induced/1000, P_profile/1000.


cruise_altitude = 2000 #m
hover_altitude  = 5000 * 0.3048


P_induced, P_profile = inducedPower(MTOW, density(0))
P_tot_SL = P_induced + P_profile

P_induced, P_profile = inducedPower(MTOW, density(hover_altitude))
P_tot_H = P_induced + P_profile

# plt.plot(DL, MTOW/P_tot_H, label= 'Hover altitude')

# plt.plot(DL, P_tot_SL, label= 'Sea level')
# plt.plot(DL, P_tot_H, label= 'Hover altitude')

plt.plot(DL, P_induced, label= 'induced')
plt.plot(DL, P_profile, label= 'profile')
# plt.plot(DL, P_tot, label= 'total')
# plt.plot(DL, r)
plt.xlabel('DL')
plt.ylabel('MTOW/P_tot')
plt.legend()
plt.show()

# print(inducedPower(MTOW, 1.225))
