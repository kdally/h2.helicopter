from math import *
import numpy as np
import matplotlib.pyplot as plt

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

S       = 33        # m2    closed wing area
w_f     = 1.7       # m     width fuselage
A       = 5.2       # -     aspect ratio
ff      = 0.85      # -     fill factor
c       = 0.2       # m     blade chord
e       = 0.8       # -     oswald factor
FM      = 0.65      # -     figure of merit
FM      = 1
c_w     = (w_f+np.sqrt(w_f**2+4*A*S))/(2*A) # m     wing chord
b       = A * c_w   # m     wing span
b_1     = (b-w_f)/2 # m     length one wing
R       = min([b_1*ff, c_w*ff])/2           # m     rotor radius
S_r     = 2*pi*R**2 # m2    rotor area
S_o     = S-S_r     # m2    open wing area


h_cr    = 4000      # m     cruise altitude
h_ho    = 1524      # m     hover altitude
# h       = 0         # m     sea level altitude

V_cr    = 97.22222  # m/s   cruise speed
V_max   = 120       # m/s   max speed
RoC     = 7         # m/s   rate of climb
t_hover = 60        # s     hover time
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


# eta_tr = 1
# eta_inst = 1
# R = 2.72
# eta_duct = 1

print('R =', R)
print('b/2 = ', (b-w_f)/2)
print('b = ', b)
print('c = ', c_w)
print('A = ', A)
print('tip speed', R*omega/343)
print('dL', MTOW/S_r)


def lift(MTOW, rho, V, S, S_o, C_lmax):
    L_mc = 0.5 * rho * V**2 * S * C_lmax             # max lift produced by closed wing
    L_mo = 0.5 * rho * V**2 * S_o * C_lmax           # max lift produced by open wing
    L_w = np.where(L_mc > MTOW, MTOW, L_mo)
    C_l = np.where(L_mc > MTOW, 2 * MTOW / (rho * V**2 * S), C_lmax)
    S_w   = np.where(L_mc > MTOW, S, S_o)

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
    if max(c_l_b) > C_lmax:
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

def torque(N, B, c_d_blade, rho, omega, R, c):
    return np.where(c_d_blade != 0.0058, N * B * c_d_blade * rho * omega**2 * R**4  * c /8, 0)


def p_climb(RoC, L_w, L_rotor, eta_duct):
    return RoC*(L_w + L_rotor/2*eta_duct)

def l_check(rho, omega, R, c, N, B, theta_t, v_i):
    return 0.5 * rho * omega**2 * R**2 * c * pi * N * B * (theta_t * R - v_i/omega)


def vals(w_f, A, S, ff):
    c_w     = (w_f+np.sqrt(w_f**2+4*A*S))/(2*A) # m     wing chord
    b       = A * c_w   # m     wing span
    b_1     = (b-w_f)/2 # m     length one wing
    R       = min([b_1*ff, c_w*ff])/2           # m     rotor radius
    S_r     = 2*pi*R**2 # m2    rotor area
    S_o     = S-S_r     # m2    open wing area
    return c_w, b, b_1, R, S_r, S_o




V = np.linspace(0.01, V_cr, 1000)
h           = h_ho

rho         = density(h)
L_w, C_l, S = lift(MTOW, rho, V, S, S_o, C_lmax)
C_d         = c_d(C_l, C_d0, A, e)
D           = drag(rho, V, S, C_d)
L_rotor     = MTOW-L_w
C_l_blade   = c_l_blade(L_rotor, N, B, rho, omega, c, R, C_lmax)
aoaob       = C_l_blade/(2*pi)
C_d_blade   = c_d_blade_1(aoaob)
v_i, Thrust = v_induced(L_rotor, rho, V, c, R, C_dr, N)

P_profile   = p_profile(N, B, C_d_blade, rho, omega, R, c)/1000
P_climb     = p_climb(0, L_w, L_rotor, eta_duct)/1000
P_induced   = p_induced(Thrust, v_i, FM, eta_tr, eta_inst, eta_duct)/1000
P_parasite  = p_parasite(D, V, eta_prop)/1000
P_level     = P_induced + P_parasite + P_profile
# print(L_w[-1], D[-1], L_rotor[-1], C_l_blade[-1], aoaob[-1], C_d_blade[-1], v_i[-1], Thrust[-1], P_profile[-1], P_induced[-1], P_parasite[-1], P_level[-1])

# print(C_l_blade)
# print(C_d)
# print(D)
# print(density(h))
# print(C_l_blade)

print(P_induced[0], P_parasite[0], P_profile[0], P_climb[0])#, P_acc, P_total_h)


print('P_level', P_level[0])
print('P_level', max(P_level))
print('D', D[-1])

""""Powers Curve Template"""

# size the figure frame: the combination [10,7.5] and the fontsizes with the scaling '0.85\linewidth' in latex works best
fig=plt.figure(figsize=(8,7))

# initialize single subplot
ax = plt.subplot(111)


# plot all the graphs in one figure
plt.plot(V, P_induced, 'r-')
plt.plot(V, P_parasite, 'r-.')
plt.plot(V, P_profile, 'r:')
plt.plot(V, P_level, 'b--')
plt.legend(['Induced Power', 'Parasite Power', 'Profile Power', 'Total Power'], fontsize=30)


# set the y tick format such that thousands are separated by commas for easier reading
ax.set_yticklabels(['{:,}'.format(int(x)) for x in ax.get_yticks().tolist()])

# set the axis labels; use latex format for units
plt.xlabel(r'Forward Speed [$\frac{m}{s}$]', fontsize=35)
plt.ylabel(r'Power [$kW$]',fontsize=35)

# set the font size, style and
# insert grid
plt.xticks(fontsize=30)
plt.yticks(fontsize=30)
plt.rc('font', family='serif')
plt.grid(True)
plt.show()

print('S_r', S_r)
print('S_r/S_o', S_r/S_o)

print('P_ind_hover', P_induced[0])
print('P_prof_hover', P_profile[0])
# save picture as eps
# plt.savefig('power_curve_embedded.eps')



# plt.plot(V, P_parasite, label= 'parasite')
# plt.plot(V, P_induced, label= 'induced')
# plt.plot(V, P_profile, label= 'profile')
# plt.plot(P_level, label= 'total')
# # plt.plot(V, torque_1)
#
# plt.legend()
# plt.show()
