import matplotlib.pyplot as plt
import math
from scipy.optimize import fmin, fsolve
# import IterationTools.control_stability.EmpennageSizing as es
import control_stability.structural_loads as st
import numpy as np
# materials
AL7075 = {'E':72.7e9, 'nu':0.33, 'sigmay':424e6, 'price':3.98, 'density': 2810, 'G':27.4e9}
AL2024 = {'E':73.1e9, 'nu':0.33, 'sigmay':360e6, 'price':2.16, 'density': 2780, 'G':28.5e9}



def class_II_cessna(MTOW, S, n_ult, A):
    return 0.04674*(MTOW)**0.397*(S)**0.360*(n_ult)**0.397*(A)**1.712


def class_II_USAF(MTOW, n_ult, A, S, t_c, V_h):
    return 96.948*((MTOW*n_ult/(10**5))**0.65*(A)**0.57*(S/100)**0.61*(2/(2*t_c))**0.36*(1+V_h/500)**0.5)**0.993


def class_II_torenbeek(MTOW, b, n_ult, S, t_r):
    return 0.00125*MTOW*b**0.75*(1+6.3/b)**0.5*n_ult**0.55*(b*S/(t_r*MTOW))**0.3


def bendingMomentStressWing(t,b,M,h,l,tspar,Istif,Astif): # h = height, # l = skin l in c direction (per section)
                                                          # t and b are skin thickness and stringer pitch (variables)
    n = 2 * (math.ceil(l/b) - 1)
    Istring = n*2*(math.pow(h/2,2)*Astif+Istif)
    Iskin = 2*(2*l*t*math.pow(h/2,2)+2*l*math.pow(t,3)/12)
    Ispar = tspar*math.pow(h,3)/12 * 4
    I = Istring + Iskin + Ispar
    return M*(h/2)/I


def equations(x, M, C, h, l, tspar, tstr, Istif, material):
    sigma_cr_else, sigma_cc, two_welse, Astif = st.Stress(C, material['nu'], material['E'], tstr, material['sigmay'])
    bendingstress = bendingMomentStressWing(x[0], x[1], M, h, l, tspar, Istif, Astif)
    cripplingstress = st.Sigma_panel(x[0], x[1], sigma_cr_else, sigma_cc, two_welse, Astif)
    f1 = bendingstress - material['sigmay']
    f2 = cripplingstress - material['sigmay']
    # print('simga else',sigma_cr_else)
    # print('sigma cc', sigma_cc*1e-6)
    # print('two_welse', two_welse)
    if x[0] <= 0 or x[0] >= 0.01:
        f1 = 1000000000000000
        f2 = 1000000000000000
    if x[1] <= 0 or x[1] >= 0.5:
        f1 = 1000000000000000
        f2 = 1000000000000000
    print('bendingstress', bendingstress*1e-6)
    print('crippling stress', cripplingstress*1e-6)
    print('f1', f1*1e-6)
    print('f2', f2*1e-6)
    return [f1, f2]


def area(t,b,h,l,tspar,Astif):
    Askin = l*t*2
    Aspar = h*tspar*2
    n = 2 * (math.ceil(l / b) - 1)
    Astif_tot = Astif*2*n
    return Askin+Aspar+Astif_tot


def mass_per_length(area,material):
    return area*material['density']


def shearFlows(Sy, T, t, tspar, b, l, h, Astif, G):
    # assume positive counter-clockwise direction, positive upward y value
    # calculate moment of inertia
    n = (math.ceil(l / b) - 1)
    Istring = n * 2 * (math.pow(h / 2, 2) * Astif)
    Iskin = 2 * (l * t * math.pow(h / 2, 2) + l * math.pow(t, 3) / 12)
    Ispar = tspar * math.pow(h, 3) / 12 * 2
    Ixx = Istring + Iskin + Ispar

    Bnormal = t * (l / n) + Astif  # stringer booms
    Bcorner = Astif + (t * (l / n) / 2) + (tspar * h / 6)  # spar cap booms
    # Bnormal = 0.0012
    # Bcorner = 0.0009
    B = [Bcorner] + [Bnormal] * n + [Bcorner] * 2 + [
        Bnormal] * n  # list of booms starting from bottom left corner to last before top left
    y = [-h / 2] * (2 + n) + [h / 2] * (1 + n)  # distances in y direction from centroid to booms
    q_base = []  # initialise list for base shear flows
    temp = 0  # temporary holder for shear flow from previous step

    # loop calculating shear flow in each gap between booms
    for i in range(len(B)):
        qb = (-Sy / Ixx) * B[i] * y[i]
        q_base.append(qb + temp)
        temp += qb
    # print('q_base', q_base)

    # set up for finding qs0
    s = [l / (n + 1)] * (n + 1) + [h] + [l / (n + 1)] * (n + 1)  # length of each interval between booms
    l_list = [h / 2] * (n + 1) + [l / 2] + [h / 2] * (n + 1)  # moment arms of different segments between booms
    sum_l = 0  # sum counter for calculating torque shear flow around area

    # sum up shear flows
    for i in range(len(l_list)):
        sum_l += q_base[i] * l_list[i] * s[i]

    qs0 = (T - sum_l) / (2 * l * h)  # shear flow at cut
    q = q_base.copy()  # create copy of base shear flows to calculate actual shear flows in segments
    q.append(0)  # base shear flow on cut segment
    s.append(h)  # length of cut segment

    # add shear flow at cut to each base shear flow
    for i in range(len(q)):
        q[i] += qs0

    # intialise rate of twist
    rate_of_twist = 0
    ts = [t] * (n + 1) + [tspar] + [t] * (n + 1) + [tspar]  # thickness of different segments

    # sum up components of loop for rate of twist
    for i in range(len(l_list)):
        rate_of_twist += q[i] * s[i] / ts[i]

    # divide by 2AG to find rate of twist
    rate_of_twist = rate_of_twist / (2 * l * h * G)

    # find stress in each section
    stress = []
    for i in range(len(q)):
        stress.append(q[i] / ts[i])
    return q, rate_of_twist, stress


def vonMisesStress(bendingstress,stress):
    return math.sqrt(bendingstress**2+6*max(stress)**2)/math.sqrt(2)

def wingloads(phase,n_hover,span,f_w,mtow,nacelle_weight,rad_weight, n_cruise, lt, a, rho_cruise, V_crit, S, c, cmac):
    # initial wing weight
    wing_weight = 300
    M = 0
    Sy = 0
    T = 0
    if phase == 'hover':
        M_thrust = n_hover * (span/2 - f_w/2) * mtow * 9.81/2
        M_nacelle = - (span/2 - f_w/2) * nacelle_weight * 9.81/2
        M_radiators = - (span/4 - f_w/2) * rad_weight * 9.81/2
        M_wings = - (span/4 - f_w/2) * wing_weight * 9.81/2
        Sy = nacelle_weight * 9.81/2 + rad_weight * 9.81/2 + wing_weight * 9.81/2 + mtow * 9.81/2*n_hover
        M = M_thrust + M_nacelle + M_radiators + M_wings
        T = 0
    # initialize loads cruise:
    if phase == 'cruise':
        Lw = (n_cruise*lt/(lt+a)*mtow*9.81+cmac*(0.5*rho_cruise*V_crit**2)*S*c/(lt+a))/2
        q_Lw = Lw/(0.75+0.25/2)
        M1 = 0.75*q_Lw * 0.75/2*span/2
        M2 = 0.25/2*q_Lw * (0.75+0.25/3)*span/2
        M_nacelle = - (span/2 - f_w/2) * nacelle_weight * 9.81/2
        M_radiators = - (span/4 - f_w/2) * rad_weight * 9.81/2
        M_wings = - (span/4 - f_w/2) * wing_weight * 9.81/2
        M = M1 + M2 + M_nacelle + M_radiators + M_wings
        Sy = Lw - nacelle_weight * 9.81/2 - rad_weight * 9.81/2 - wing_weight * 9.81/2
        T = Lw*0.2
    return M, Sy, T


def wingstructure(M, Sy, T, span):

    C = 4
    tstr = 0.003
    h = 0.3
    l = 0.5
    tspar = 0.0025
    Istif = 0
    material = AL2024


    # initialize loads hover:


    # initialize variables
    xlow, xhigh = 0.0008, 0.003
    yhigh, ylow = 0.4, 0.05
    t = np.linspace(xhigh, xlow, 200)
    b = np.linspace(yhigh, ylow, 200)

    # initialize lists
    tlist = []
    blist = []
    bs = []
    cs = []
    sol = []
    aarr = []
    mass = []
    vmises = []
    sigma_cr_else, sigma_cc, two_welse, Astif = st.Stress(C, material['nu'], material['E'], tstr, material['sigmay'])

    for i in t:
        for j in b:
            bendingstress = bendingMomentStressWing(i, j, M, h, l, tspar, Istif, Astif)
            cripplingstress = st.Sigma_panel(i, j, sigma_cr_else, sigma_cc, two_welse, Astif)
            q, rateoftwist, stress = shearFlows(Sy, T, i, tspar, j, l, h, Astif, material['G'])
            vm = vonMisesStress(bendingstress,stress)
            a = area(i, j, h, l, tspar, Astif)
            m = a * material['density'] * span

            tlist.append(i)
            blist.append(j)
            bs.append(bendingstress / (1e6))
            cs.append(cripplingstress / (1e6))
            sol.append(cripplingstress / (1e6) - bendingstress / (1e6))
            aarr.append(a)
            mass.append(m)
            vmises.append(vm / (1e6))
            # score.append(s)

    # reshape all

    np_tlist = np.array(tlist).reshape(200, 200)
    np_blist = np.array(blist).reshape(200, 200)
    np_vmises = np.array(vmises).reshape(200, 200)
    np_bs    = np.array(bs).reshape(200,200)
    np_cs    = np.array(cs).reshape(200,200)
    np_mass = np.array(mass).reshape(200, 200)

    # evaluate:

    final = np.where(np_cs - np_bs < 0.0, 100000, np_mass)
    final2 = np.where(np_vmises > material['sigmay'] / (1e6), 100000, final)
    ind = np.unravel_index(np.argmin(final2), final2.shape)

    # print('final mass', final2[ind])

    fmass = final2[ind]
    finalt = np_tlist[ind]
    print('thickness',finalt)
    finalb = np_blist[ind]
    print('pitch',finalb)
    print('von Mises',np_vmises[ind])
    print('bending stress',np_bs[ind])
    print('crippling',np_cs[ind])
    area_skin = l*2*finalt
    area_spar = h*tspar*4
    n = (math.ceil(l / finalb) - 1)
    area_stringers = n * 2 * Astif
    rib_weight = 10*2*1.5*0.25*material['density']*0.002
    open_section_weight = span*0.45*0.0015*material['density']*2
    open_section_area   = 0.45*finalt*2
    final_mass = fmass + rib_weight + open_section_weight





    # plt.subplot(2, 2, 1)
    # plt.title('Von Mises Stress [MPa]')
    # plt.xlabel('Thickness [mm]')
    # plt.ylabel('Stiffner Pitch [m]')
    # plt.scatter([tlist[i]*1000 for i in range(len(tlist))], blist, c=vmises, cmap='Oranges')
    # plt.colorbar()
    # plt.xlim(xlow*1000, xhigh*1000)
    # plt.ylim(ylow, yhigh)
    # # plt.show()
    # plt.subplot(2, 2, 2)
    # plt.title('Buckling Stress [MPa]')
    # plt.xlabel('Thickness [mm]')
    # plt.ylabel('Stiffner Pitch [m]')
    # plt.scatter([tlist[i]*1000 for i in range(len(tlist))], blist, c=cs, cmap='Oranges')
    # plt.colorbar()
    # plt.xlim(xlow*1000, xhigh*1000)
    # plt.ylim(ylow, yhigh)
    # # plt.show()
    # plt.subplot(2, 2, 3)
    # plt.title('Buckling - Bending Stresss Difference [MPa]')
    # plt.xlabel('Thickness [mm]')
    # plt.ylabel('Stiffner Pitch [m]')
    # plt.scatter([tlist[i]*1000 for i in range(len(tlist))], blist, c=sol, cmap='RdBu', vmin=-300, vmax=300)
    # plt.colorbar()
    # plt.xlim(xlow*1000, xhigh*1000)
    # plt.ylim(ylow, yhigh)
    # # plt.show()
    # plt.subplot(2, 2, 4)
    # plt.title('Wingbox mass [MPa]')
    # plt.xlabel('Thickness [mm]')
    # plt.ylabel('Stiffner Pitch [m]')
    # plt.scatter([tlist[i]*1000 for i in range(len(tlist))], blist, c=mass, cmap='Blues')
    # plt.colorbar()
    # plt.scatter(np_tlist[ind]*1000, np_blist[ind], c='r')
    # plt.xlim(xlow*1000, xhigh*1000)
    # plt.ylim(ylow, yhigh)
    # plt.show()

    return final_mass, [fmass, rib_weight, open_section_weight], [area_skin, area_spar, area_stringers,
                                                                  open_section_area], np_tlist[ind], np_blist[ind]

# if __name__ == '__main__':
#
#     M = 5 * 3.5 * 4000 * 9.81 / 2
    # T = 4000 * 9.81 / 2 * 3.8 * 0.2
    # Sy = 4000 * 9.81 / 2 * 3.8
    # span = 11.61
    # M_c = 145213.103206156
    # Sy_c = 60253.4071766799
    # T_c = 13132.2408415760
    # b = 11.165728974548776
    # wing_mass_c, mass_fractions_c, areas_c = wingstructure(M_c, Sy_c, T_c, b)
    # print(wing_mass_c)
    # print(mass_fractions_c)
    # print(areas_c)
    # wingloads(phase, n_hover, span, f_w, mtow, nacelle_weight, rad_weight, n_cruise, lt, a, rho_cruise, V_crit, S, c)
    # print(wingstructure(M, Sy, T, span))
    # # initialize values
    # MTOW = 4000*2.20462 # kg --> lbs
    # S    = 22*10.7639 # m2 --> ft2
    # n_ult= 3.8 # -
    # A    = 5.37 # -
    # t_c  = 0.18 # -
    # t_r  = 0.36*3.28084 # m --> ft
    # b    = 10.74*3.28084 # m --> ft
    # V_h  = 50*1.94384 # m/s --> kts
    # cessna = class_II_cessna(MTOW, S, n_ult, A)/2.20462
    # USAF   = class_II_USAF(MTOW, n_ult, A, S, t_c, V_h)/2.20462
    # toren  = class_II_torenbeek(MTOW, b, n_ult, S, t_r)/2.20462
    # print((cessna+toren+USAF)/3)
    # plt.bar([1,2,3,4],[cessna,USAF,toren,243.55])
    # plt.ylabel('Structural Weight [kg]')
    # plt.tick_params(axis='x', which='both', labelbottom=False)
    # plt.show()


