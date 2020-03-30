from math import *
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def density(height):
    temp = 288.15-0.0065*height
    pressure = 101325*pow(temp/288.15,(-9.80665/(-0.0065*287.05)))
    return pressure/(287.05*temp)

MTOM    = 4000      # kg
g       = 9.80665   # m/s2
MTOW    = MTOM * g  # N

C_lcr   = 0.4       # -     cruise c_l
C_lmax  = 1.2       # -     max c_l
C_d0    = 0.02      # -     zero lift drag coefficient
C_dr    = 0.0       # -     drag coefficient downwards air

S_t     = 33        # m2    closed wing area
w_f     = 1.6       # m     width fuselage
A       = 6         # -     aspect ratio
ff      = 0.95      # -     fill factor
c       = 0.2       # m     blade chord
e       = 0.5       # -     oswald factor
FM      = 0.65      # -     figure of merit
FM      = 1
c_w     = (w_f+np.sqrt(w_f**2+4*A*S_t))/(2*A) # m     wing chord
b       = A * c_w   # m     wing span
b_1     = (b-w_f)/2 # m     length one wing
R       = min([b_1*ff, c_w*ff])/2           # m     rotor radius
S_r     = 2*pi*R**2 # m2    rotor area
S_o     = S_t-S_r     # m2    open wing area

t_ho    = 60
t_cl    = 500
t_cr    = 2500
t_de    = 500

h_cr    = 4000      # m     cruise altitude
h  = h_cr
h_h     = 1524      # m     hovering altitude
V_cr    = 111.111        # m/s   cruise speed
RoC     = 7         # m/s   rate of climb
V_climb = 220/3.6   # m/s   climb speed
t_a_climb = 80      # s     time to accelerate to climb speed
t_a_cruise = 80     # s     time to accelerate to cruise speed
a_climb = V_climb/t_a_climb # m/s2 acceleration to climb speed
a_cruise = (V_cr - V_climb) / t_a_cruise # m/s2 acceleration to cruise speed
xrange  = 300000    # m    cruise range

omega   = 150       # rad/s angular speed rotor
N       = 2         # -     number of rotors
B       = 30         # -     blades per rotor
theta_t = 10*pi/180 # rad   tip pitch
theta_r = 18*pi/180 # rad   root pitch
theta_1 = theta_t - theta_r

eta_tr  = 1         # -     efficiency transmission
eta_inst= 0.98      # -     efficiency installation
eta_duct= 0.38       # -     efficiency duct
eta_prop= 0.8       # -     efficiency back propeller


print('R =', R)
print('b/2 = ', (b-w_f)/2)
print('b = ', b)
print('c = ', c_w)
print('A = ', A)
print('tip speed', R*omega/343)
print('dL', MTOW/S_r)


def lift(MTOW, rho, V, S, S_o, C_lmax):
    if V == 0:
        return 0, 0, S_o
    L_mc = 0.5 * rho * V**2 * S * C_lmax             # max lift produced by closed wing
    L_mo = 0.5 * rho * V**2 * S_o * C_lmax           # max lift produced by open wing

    L_w = MTOW if L_mc > MTOW else L_mo
    C_l = 2 * MTOW / (rho * V**2 * S) if L_mc > MTOW else C_lmax
    S_w = S if L_mc > MTOW else S_o

    # L_w = np.where(L_mc > MTOW, MTOW, L_mo)
    # C_l = np.where(L_mc > MTOW, 2 * MTOW / (rho * V**2 * S), C_lmax)
    # S_w   = np.where(L_mc > MTOW, S, S_o)

    return L_w, C_l, S_w
    # return L_mo, C_lmax, S_o

def c_d(C_l, C_d0, A, e):
    return C_d0 + C_l**2 / (pi * A * e)

def drag(rho, V, S, C_d, ):
    return  0.5 * rho * V**2 * S * C_d


def c_l_blade(L_rotor, N, B, rho, omega, c, R, C_lmax):
    c_l_b = 6*L_rotor / (N*B*rho*omega**2*c*R**3)
    # if c_l_b > C_lmax:
    #     print('Simulation not correct, value of C_l is too big')
    if c_l_b > C_lmax:
        print('Error: Simulation not valid')
        print('C_l blade is too big: ', max(c_l_b))
    return c_l_b

def collective(L_rotor, N, B, rho, omega, c, R, theta_t, v_i):
    num_1 = L_rotor / (N * B * rho * omega**2 * pi * c)
    num_2 = R**2 * (theta_t * R - v_i / omega) / 2
    col   = (num_1 - num_2) * R**3 / 3
    return np.where(v_i != 0, col, 0)

def v_induced(L_rotor, rho, V, c, R, C_dr, N):
    A_p = N * c * R
    A = N * pi * R**2
    Thrust = (1 + (C_dr * A_p)/A) * L_rotor
    # v_induced = np.sqrt(-((V**2)/2) + np.sqrt((V**4)/4 + (Thrust/(2 * rho * A))**2))
    v_induced = np.sqrt(Thrust/(2*rho*A))
    return np.where(Thrust != 0, v_induced, 0), Thrust

def aoa_cuad(theta_t, R, v_i, omega, alfa_c):
    return theta_t * R - v_i / omega + alfa_c

def aoa_lin(rho, R_o, R_f, omega, c, v_i, N, B, L_rotor, theta_t, theta_r):
    return (L_rotor / (N * B * rho * omega**2 * c * pi) - theta_r * (R_f**3/3 - R_o**3/3) - theta_t * (R_f**3/4 -R_o**3/4) + v_i*(R_f**2/(2*omega)- R_o**2/(2*omega))) / (R_f**3/3 - R_o**3/3)

def omegas_cuad(rho, R, c, v_i, N, B, L_rotor, theta_t):
    A = 0.5 * rho * R**3 * c * pi * theta_t * N * B
    B = -0.5 * rho * R**2 * c * pi * v_i * N * B
    C = -L_rotor
    omega_1 = (-B - np.sqrt(B**2 - 4 * A * C)) / (2 * A)
    omega_2 = (-B + np.sqrt(B**2 - 4 * A * C)) / (2 * A)
    return omega_1, omega_2

def omegas_lin(rho, R, c, v_i, N, B, L_rotor, theta_t, theta_r):
    A = N * B * rho * c * pi * (R**3 * theta_r / 3 + R**3 * theta_t / 4)
    B = N * B * rho * c * pi * (R**2 * v_i/ 2)
    C = -L_rotor
    omega_1 = (-B - np.sqrt(B**2 - 4 * A * C)) / (2 * A)
    omega_2 = (-B + np.sqrt(B**2 - 4 * A * C)) / (2 * A)
    return omega_1, omega_2

def c_d_blade(alfa_c, omega, R, v_i, theta_t):
    x = theta_t * R - v_i / omega + alfa_c
    c_d = 20441 * x**6 -1505.9 * x**5 + 442.64 *x**4 - 63.716 * x**3 + 4.7524 * x**2 - 0.0961 * x + 0.0058
    return np.where(v_i != 0, c_d, 0)

def c_d_blade_1(x):
    c_d = 20441 * x**6 -1505.9 * x**5 + 442.64 *x**4 - 63.716 * x**3 + 4.7524 * x**2 - 0.0961 * x + 0.0058
    return c_d

def p_parasite(D, V, eta_prop):
    p_p = D * V/eta_prop
    return p_p

def p_induced(Thrust, v_induced, FM, eta_tr, eta_inst, eta_duct):
    return Thrust*v_induced/FM * (1/(eta_tr*eta_inst))*eta_duct*1.15**1.5

def p_profile(N, B, c_d_blade, rho, omega, R, c):
    return np.where(c_d_blade != 0.0058, N * B * c_d_blade * rho * omega**3 * R**4  * c /8, 0)

def p_acceleration(F_acc, V, eta_prop):
    return F_acc*V/eta_prop

def torque(N, B, c_d_blade, rho, omega, R, c):
    return np.where(c_d_blade != 0.0058, N * B * c_d_blade * rho * omega**2 * R**4  * c /8, 0)


def p_climb(RoC, L_w, L_rotor, eta_duct, rho, S_r):
    # v_i = -RoC/2 + np.sqrt((RoC/2)**2 + L_rotor/(2*rho*S_r))
    # p_r = L_rotor * v_i
    # p_w = L_w * RoC
    # return p_r + p_w

    return RoC*(L_w + L_rotor/2*eta_duct)

def l_check(rho, omega, R, c, N, B, theta_t, v_i):
    return 0.5 * rho * omega**2 * R**2 * c * pi * N * B * (theta_t * R - v_i/omega)


# V = np.linspace(0.01, 120, 1000)
#
# rho         = density(h)
# L_w, C_l, S = lift(MTOW, rho, V, S, S_o, C_lmax)
# C_d         = c_d(C_l, C_d0, A, e)
# D           = drag(rho, V, S, C_d)
# L_rotor     = MTOW-L_w
# C_l_blade   = c_l_blade(L_rotor, N, B, rho, omega, c, R, C_lmax)
# aoaob       = C_l_blade/(2*pi)
# C_d_blade   = c_d_blade_1(aoaob)
# v_i, Thrust = v_induced(L_rotor, rho, V, c, R, C_dr, N)
#
# P_profile   = p_profile(N, B, C_d_blade, rho, omega, R, c)/1000
# P_induced   = p_induced(Thrust, v_i, FM, eta_tr, eta_inst, eta_duct)/1000
# P_climb     = p_climb(0, L_w, L_rotor, eta_duct)/1000
# P_parasite  = p_parasite(D, V, eta_prop)/1000
# P_level     = P_induced + P_parasite + P_profile



# Hover
rho = density(h_h)
V   = 0
RoC = 0

L_w, C_l, S = lift(MTOW, rho, V, S_t, S_o, C_lmax)
C_d         = c_d(C_l, C_d0, A, e)
D           = drag(rho, V, S, C_d)
L_rotor     = MTOW-L_w
C_l_blade   = c_l_blade(L_rotor, N, B, rho, omega, c, R, C_lmax)
aoaob       = C_l_blade/(2*pi)
C_d_blade   = c_d_blade_1(aoaob)
v_i, Thrust = v_induced(L_rotor, rho, V, c, R, C_dr, N)

P_profile   = p_profile(N, B, C_d_blade, rho, omega, R, c)/1000
P_induced   = p_induced(Thrust, v_i, FM, eta_tr, eta_inst, eta_duct)/1000
P_parasite  = p_parasite(D, V, eta_prop)/1000
P_climb     = p_climb(RoC, L_w, L_rotor, eta_duct, rho, S_r)/1000
F_acc       = 0
P_acc       = 0
P_total_h   = P_induced + P_parasite + P_profile + P_climb + P_acc
print(P_induced, P_parasite, P_profile, P_climb, P_acc, P_total_h)

# plt.bar(1, P_induced)
# plt.bar(2, P_profile)
# plt.bar(3, P_parasite)
# plt.bar(4, P_climb)
# plt.bar(5, P_acc)
# # plt.bar(3, P_induced+P_profile)
# plt.show()

# # rho         = density(h)
# # L_w, C_l, S = lift(MTOW, rho, V, S, S_o, C_lmax)
# # C_d         = c_d(C_l, C_d0, A, e)
# # D           = drag(rho, V, S, C_d)
# # L_rotor     = MTOW-L_w
# # C_l_blade   = c_l_blade(L_rotor, N, B, rho, omega, c, R, C_lmax)
# # aoaob       = C_l_blade/(2*pi)
# # C_d_blade   = c_d_blade_1(aoaob)
# # v_i, Thrust = v_induced(L_rotor, rho, V, c, R, C_dr, N)
# #
# # P_profile   = p_profile(N, B, C_d_blade, rho, omega, R, c)/1000
# # P_induced   = p_induced(Thrust, v_i, FM, eta_tr, eta_inst, eta_duct)/1000
# # P_climb     = p_climb(0, L_w, L_rotor, eta_duct)/1000
# # P_parasite  = p_parasite(D, V, eta_prop)/1000
# # P_level     = P_induced + P_parasite + P_profile
#
#
# Climb at speed
rho = density(0)
V   = V_climb
a = a_climb
RoC = 8
F_acc = MTOM * a

L_w, C_l, S = lift(MTOW, rho, V, S_t, S_o, C_lmax)
C_d         = c_d(C_l, C_d0, A, e)
D           = drag(rho, V, S, C_d)
L_rotor     = MTOW-L_w
C_l_blade   = c_l_blade(L_rotor, N, B, rho, omega, c, R, C_lmax)
aoaob       = C_l_blade/(2*pi)
C_d_blade   = c_d_blade_1(aoaob)
v_i, Thrust = v_induced(L_rotor, rho, V, c, R, C_dr, N)
P_profile   = p_profile(N, B, C_d_blade, rho, omega, R, c)/1000
P_induced   = p_induced(Thrust, v_i, FM, eta_tr, eta_inst, eta_duct)/1000
P_parasite  = p_parasite(D, V, eta_prop)/1000
P_acc       = p_acceleration(F_acc, V, eta_prop)/1000
P_climb     = p_climb(RoC, L_w, L_rotor, eta_duct, rho, S_r)/1000
P_total_cf  = P_induced + P_parasite + P_profile + P_climb + P_acc

# print('MTOW', MTOW)
# print('L_w, C_l, S', L_w, C_l, S)
# print('C_d', C_d)
# print('D', D)
# print('L-r', L_rotor)
# print("C-l-b", C_l_blade)
# print('aoaob', aoaob)
# print(v_i, Thrust)

# print(P_induced, P_parasite, P_profile, P_climb, P_acc)

# print(P_induced, P_parasite, P_profile, P_climb, P_acc)

# Climb vertically
rho = density(h_h)
V   = 0
a = a_climb
a = 0
RoC = 8
F_acc = MTOM * a

L_w, C_l, S = lift(MTOW, rho, V, S, S_o, C_lcr)
C_d         = c_d(C_l, C_d0, A, e)
D           = drag(rho, V, S, C_d)
L_rotor     = MTOW-L_w
C_l_blade   = c_l_blade(L_rotor, N, B, rho, omega, c, R, C_lmax)
aoaob       = C_l_blade/(2*pi)
C_d_blade   = c_d_blade_1(aoaob)
v_i, Thrust = v_induced(L_rotor, rho, V, c, R, C_dr, N)
P_profile   = p_profile(N, B, C_d_blade, rho, omega, R, c)/1000
P_induced   = p_induced(Thrust, v_i, FM, eta_tr, eta_inst, eta_duct)/1000
P_parasite  = p_parasite(D, V, eta_prop)/1000
P_acc       = p_acceleration(F_acc, V, eta_prop)/1000
# print('L_w', L_w)
P_climb     = p_climb(RoC, L_w, L_rotor, eta_duct, rho, S_r)/1000
P_total_cv  = P_induced + P_parasite + P_profile + P_climb + P_acc

# print(P_induced, P_parasite, P_profile, P_climb, P_acc, P_total_cv)
# P_total_cv = 0


# Cruise
rho = density(h_cr)
V   = V_cr
F_acc = 0
RoC = 0

L_w, C_l, S = lift(MTOW, rho, V, S_t, S_o, C_lmax)
C_d         = c_d(C_l, C_d0, A, e)
D           = drag(rho, V, S, C_d)
L_rotor     = MTOW-L_w
C_l_blade   = c_l_blade(L_rotor, N, B, rho, omega, c, R, C_lmax)
aoaob       = C_l_blade/(2*pi)
C_d_blade   = c_d_blade_1(aoaob)
v_i, Thrust = v_induced(L_rotor, rho, V, c, R, C_dr, N)
P_profile   = p_profile(N, B, C_d_blade, rho, omega, R, c)/1000
P_induced   = p_induced(Thrust, v_i, FM, eta_tr, eta_inst, eta_duct)/1000
P_parasite  = p_parasite(D, V, eta_prop)/1000
P_acc       = p_acceleration(F_acc, V, eta_prop)/1000
P_climb     = p_climb(0, L_w, L_rotor, eta_duct, rho, S_r)/1000
P_total_cr  = P_induced + P_parasite + P_profile + P_climb + P_acc




# Descent at speed
rho = density(h_cr)
V   = V_cr
a = 0
RoC = -8
F_acc = MTOM * a

L_w, C_l, S = lift(MTOW, rho, V, S_t, S_o, C_lmax)
C_d         = c_d(C_l, C_d0, A, e)
D           = drag(rho, V, S, C_d)
L_rotor     = MTOW-L_w
C_l_blade   = c_l_blade(L_rotor, N, B, rho, omega, c, R, C_lmax)
aoaob       = C_l_blade/(2*pi)
C_d_blade   = c_d_blade_1(aoaob)
v_i, Thrust = v_induced(L_rotor, rho, V, c, R, C_dr, N)
P_profile   = p_profile(N, B, C_d_blade, rho, omega, R, c)/1000
P_induced   = p_induced(Thrust, v_i, FM, eta_tr, eta_inst, eta_duct)/1000
P_parasite  = p_parasite(D, V, eta_prop)/1000
P_acc       = p_acceleration(F_acc, V, eta_prop)/1000
P_climb     = p_climb(RoC, L_w, L_rotor, eta_duct, rho, S_r)/1000
P_total_de  = P_induced + P_parasite + P_profile + P_climb + P_acc

# print(P_induced, P_parasite, P_profile, P_climb, P_acc)



print('Energy total = ', (P_total_cr * t_cr + P_total_cf * t_cl + P_total_h * t_ho + P_total_de * t_de)/1000)
print('Max power = ', np.max([P_total_cr, P_total_cf, P_total_h, P_total_de, P_total_cv]))


# plt.bar(1, P_total_h, label= 'hover')
# plt.bar(2, P_total_cl, label= 'climb')
# plt.bar(3, P_total_cr, label= 'cruise')

# size the figure frame: the combination [10,7.5] and the fontsizes with the scaling '0.85\linewidth' in latex works best
fig=plt.figure(figsize=(13,8))

# initialize single subplot
ax = plt.subplot(111)

# set the values to plt in the bar chart

power_phases=pd.Series({'Hover': P_total_h,'Climb forward ':P_total_cf, 'Cruise':P_total_cr,'Descent':P_total_de})

power_phases.plot(kind='bar')

# set the y tick format such that thousands are separated by commas for easier reading
ax.set_yticklabels(['{:,}'.format(int(x)) for x in ax.get_yticks().tolist()])

# set the axis labels; use latex format for units
plt.xlabel(r'Flight Phases', fontsize=18)
plt.ylabel(r'Power ($kW$)',fontsize=18)

# set the font size, style and insert grid
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.grid(True)
# plt.savefig('power_phases.eps')


# plt.legend()
plt.show()
