import numpy as np
from math import pi
from scipy.signal import find_peaks
from scipy.signal import argrelextrema
import matplotlib.pyplot as plt
import pandas as pd
import os
import pickle
import progressbar
import time
import warnings


def main():
    pass

if __name__ == "__main__":
    from IterationTools.aerodynamics.BET import Rotor

def ISA(h):
    temp = 288.15-0.0065*h
    pressure = 101325*np.power(temp/288.15,(-9.80665/(-0.0065*287.05)))
    rho = pressure/(287.05*temp)
    return rho


def ISA(h):
    temp = 288.15-0.0065*h
    pressure = 101325*np.power(temp/288.15,(-9.80665/(-0.0065*287.05)))
    rho = pressure/(287.05*temp)
    v_sound = np.sqrt(1.4 * 287.05 * temp)
    v_tip_max = 0.75 * v_sound
    maxOmega = v_tip_max/rotor.R
    print('v_tip', v_tip_max)

    return rho

def editPowers(powers, MTOW):
    climbvertical1 = MTOW*9.80665*8/2000
    climbvertical2 = MTOW*9.806665*1.5/2000
    climbwing = MTOW*9.80665*8/1000
    powers[2] += climbvertical1
    powers[3] += climbwing
    powers[4] += climbwing
    powers[9] -= climbwing
    powers[10] -= climbwing
    powers[11] -= climbvertical2
    return np.array(powers)+11



def calcTotal(N, lst):
    return np.array(lst)*N

def initVariables(R):
    N = 2
    B = 3
    # B = 2
    R = R
    # R = 3.81
    cutout = 0.2
    solidity = 0.02
    # solidity = 0.06928
    # solidity = 0.027
    # solidity = 0.065
    theta_tip = 1*pi/180
    taper = 1
    airfoil = '0012'
    be = 100

    return N, B, R, cutout, solidity, theta_tip, taper, airfoil, be

def genLists():
    powers = []
    omegas = []
    thrusts = []
    t_ratios = []
    dcts = []
    cts = []
    cps = []
    inflows = []
    cols = []
    twists = []
    t_index = []

    return powers, omegas, thrusts, t_ratios, dcts, cts, cps, inflows, cols, twists, t_index

def genLists_i():
    power_i = []
    thrust_i = []
    omega_i = []
    ct_i = []
    cp_i = []
    col_i = []
    inflow_i = []
    index = []

    return power_i, thrust_i, omega_i, ct_i, cp_i, col_i, inflow_i, index

def genNewLists_i():
    omega_i = []
    ct_i = []
    cp_i = []
    index = []

    return omega_i, ct_i, cp_i, index

def genTwist(theta_tip, cutout, be):
    twists = []
    r = np.linspace(cutout, 1, be)[:-1]

    # Ideal
    # for i in np.linspace(0, 10, 20):
    #     twists.append(theta_tip/r*i)

    #Flat
    # twists.append(np.array([0]*len(r)))

    # Linear
    for i in np.linspace(5, 35, 5):
        theta_root = i*pi/180
        theta_tip = 0
        slope = (theta_tip-theta_root)/(1-cutout)
        twists.append(theta_root+slope*(r-cutout))
    # solidity
    sols = []
    for i in np.arange(0.01, 0.041, 0.01):
    # for i in [0.02]:
        sols.append(i)

    twist_index = np.arange(0, len(twists)+len(sols), 1)
    return twists, sols, twist_index

def hoverFlight(h, R, T, wing_chord):
    MTOW = T
    N, B, R, cutout, solidity, theta_tip, taper, airfoil, be = initVariables(R)
    powers, omegas, thrusts, t_ratios, dcts, cts, cps, inflows, cols, twists, t_index = genLists()
    h = h
    v_inf = 0
    differenttwists, twist_index = genTwist(theta_tip, cutout, be)
    rotor = Rotor(N, B, R, cutout, taper, solidity, theta_tip, airfoil, v_inf, be, h)

    for t_i, twist in enumerate(differenttwists):
        print(t_i+1, ' / ', len(differenttwists))
        power_i, thrust_i, omega_i, ct_i, cp_i, col_i, inflow_i, index  = genLists_i()
        rotor.twist = twist
        # plt.plot(rotor.r, twist)
        # plt.show()

        for omega in np.linspace(20, rotor.maxOmega, 20):
            # goal = MTOW/N/(rotor.rho*pi*rotor.R**2*(omega*rotor.R)**2)
            goal = (MTOW+0.3*MTOW*wing_chord*N/(pi*rotor.R))/N/(rotor.rho*pi*rotor.R**2*(omega*rotor.R)**2)

            power, thrust, omega, ct, cp, inflow, collective = flight(rotor, goal, omega)
            power_i.append(power)
            thrust_i.append(thrust)
            omega_i.append(omega)
            ct_i.append(ct)
            cp_i.append(cp)
            inflow_i.append(inflow)
            col_i.append(collective)



        if min(cp_i) != 10**10:
            opt_list = np.array(cp_i)/np.array(ct_i)
            opt_index = list(opt_list).index(min(opt_list))

            # print(rotor.twist+collective-inflow/rotor.r)
            powers.append(power_i[opt_index])
            thrusts.append(thrust_i[opt_index])
            omegas.append(omega_i[opt_index])
            # t_ratios.append(t_ratio)
            cts.append(ct_i[opt_index])
            cps.append(cp_i[opt_index])
            twists.append(rotor.twist)
            inflows.append(inflow_i[opt_index])
            cols.append(col_i[opt_index])

            aoa_r = rotor.twist+col_i[opt_index]-inflow_i[opt_index]/rotor.r

            lift = B*np.sum((rotor.airfoil.cl_alpha*aoa_r*0.5*rotor.rho*((inflow*omega*rotor.R)**2+omega**2*(rotor.r*rotor.R)**2)*rotor.dr*rotor.c))
            thrust_from_ct = ct*0.5*rotor.rho*(omega*rotor.R)**2*pi*rotor.R**2
            t_index.append(twist_index[t_i])
            # print('lift', lift)
            # print('thrust', thrust_from_ct)
        else:
            print('Non valid twist ratio')



    fig, axs = plt.subplots(4,1)

    axs[0].set_title('Hover')
    axs[0].plot(t_index, np.array(powers)/1000*N)
    axs[0].set_ylabel('power [kW]')
    # plt.show()

    axs[1].plot(t_index, np.array(thrusts)*N)
    axs[1].set_ylabel('thrust [N]')
    # plt.show()

    axs[2].plot(t_index, omegas)
    axs[2].set_ylabel('omega [rad/s]')


    k = 1.15
    cd0 = 0.011
    cps_theory = k*np.array(cts)**(3/2)/(np.sqrt(2))+(solidity*cd0/8)
    # print('theory', cps_theory)
    ratios = np.array(cts)/np.array(cps)

    axs[3].scatter(cps, cts, label='single')
    axs[3].plot(cps_theory, cts, label='theory')
    axs[3].scatter(cps[np.where(ratios == np.max(ratios))[0][0]], cts[np.where(ratios == np.max(ratios))[0][0]], label='Opt point')
    axs[3].set_ylabel('ct to cp')
    axs[3].set_xlim(min(min(cps), min(cps_theory))-0.00005, max(max(cps), max(cps_theory))+0.00005)
    axs[3].set_ylim(min(cts)-0.001, max(cts)+0.001)
    axs[3].legend()
    # plt.show()

    opt_index = powers.index(min(powers))
    collective = cols[opt_index]
    inflow = inflows[opt_index]
    omega = omegas[opt_index]
    twist = twists[opt_index]
    # dct = dcts[opt_index]

    aoa_r = twist+collective-np.arctan(inflow/rotor.r)
    lift = (rotor.airfoil.cl_alpha*aoa_r*0.5*rotor.rho*((inflow*omega*rotor.R)**2+omega**2*(rotor.r*rotor.R)**2)*rotor.dr*rotor.c)

    plt.figure(1)
    fig2, ax = plt.subplots(4,1)

    ax[0].set_title('Hover. Optimum twist ratio: '+ str(t_index[opt_index]))
    ax[0].plot(rotor.r, inflow)
    ax[0].set_ylabel('inflow')

    ax[1].plot(rotor.r, aoa_r*180/pi)
    ax[1].set_ylabel('angle of attack')

    ax[2].plot(rotor.r, twist*180/pi, label= 'twist')
    ax[2].plot(rotor.r, [collective*180/pi]*np.int(len(rotor.r)), label= 'collective')
    ax[2].plot(rotor.r, np.arctan(inflow/rotor.r)*180/pi, label='phi')
    ax[2].legend()
    ax[2].set_ylabel('dct')

    ax[3].plot(rotor.r, twist*180/pi)
    ax[3].set_ylabel('twist')
    # axs[3].plot(rotor.r, np.sqrt((rotor.r*omega)**2+(inflow*omega*rotor.R)**2))
    # ax[3].plot(cp, ct)
    # plt.show()
    return powers, t_index

def forwardFlight(h, R, T, v_inf):
    T = T
    N, B, R, cutout, solidity, theta_tip, taper, airfoil, be = initVariables(R)
    powers, omegas, thrusts, t_ratios, dcts, cts, cps, inflows, cols, twists, t_index = genLists()
    h = h
    differenttwists, twist_index = genTwist(theta_tip, cutout, be)
    for t_i, twist in enumerate(differenttwists):
        print(t_i+1, ' / ', len(differenttwists))

        power_i, thrust_i, omega_i, ct_i, cp_i, col_i, inflow_i, index  = genLists_i()
        rotor = Rotor(N, B, R, cutout, taper, solidity, theta_tip, airfoil, v_inf, be, h)
        rotor.twist = twist
        # plt.plot(rotor.r, twist)
        # plt.show()

        for omega in np.linspace(20, rotor.maxOmega, 20):
            # goal = (MTOW+0.3*MTOW*2*N/(pi*rotor.R))/N/(rotor.rho*pi*rotor.R**2*(omega*rotor.R)**2)
            goal = T/N/(rotor.rho*pi*rotor.R**2*(omega*rotor.R)**2)

            power, thrust, omega, ct, cp, inflow, collective = flight(rotor, goal, omega)
            power_i.append(power)
            thrust_i.append(thrust)
            omega_i.append(omega)
            ct_i.append(ct)
            cp_i.append(cp)
            inflow_i.append(inflow)
            col_i.append(collective)


        if min(cp_i) != 10**10:
            opt_list = np.array(cp_i)/np.array(ct_i)
            opt_index = list(opt_list).index(min(opt_list))

            # print(rotor.twist+collective-inflow/rotor.r)
            powers.append(power_i[opt_index])
            thrusts.append(thrust_i[opt_index])
            omegas.append(omega_i[opt_index])
            # t_ratios.append(t_ratio)
            cts.append(ct_i[opt_index])
            cps.append(cp_i[opt_index])
            twists.append(rotor.twist)
            inflows.append(inflow_i[opt_index])
            cols.append(col_i[opt_index])

            aoa_r = rotor.twist+col_i[opt_index]-inflow_i[opt_index]/rotor.r

            # lift = B*np.sum((rotor.airfoil.cl_alpha*aoa_r*0.5*rotor.rho*((inflow*omega*rotor.R)**2+omega**2*(rotor.r*rotor.R)**2)*rotor.dr*rotor.c))
            # thrust_from_ct = ct*0.5*rotor.rho*(omega*rotor.R)**2*pi*rotor.R**2
            t_index.append(twist_index[t_i])
            # print('lift', lift)
            # print('thrust', thrust_from_ct)
        else:
            print('Non valid twist ratio')



    fig, axs = plt.subplots(4,1)
    axs[0].set_title('Forward flight')
    axs[0].plot(t_index, 2*np.array(powers)/1000)
    axs[0].set_ylabel('power [kW]')
    # plt.show()

    axs[1].plot(t_index, 2*np.array(thrusts))
    axs[1].set_ylabel('thrust [N]')
    # plt.show()

    axs[2].plot(t_index, omegas)
    axs[2].set_ylabel('omega [rad/s]')


    k = 1.15
    cd0 = 0.011
    cps_theory = k*np.array(cts)**(3/2)/(np.sqrt(2))+(solidity*cd0/8)
    # print('theory', cps_theory)
    ratios = np.array(cts)/np.array(cps)

    axs[3].scatter(cps, cts, label='single')
    axs[3].plot(cps_theory, cts, label='theory')
    axs[3].scatter(cps[np.where(ratios == np.max(ratios))[0][0]], cts[np.where(ratios == np.max(ratios))[0][0]], label='Opt point')
    axs[3].set_ylabel('ct to cp')
    axs[3].set_xlim(min(min(cps), min(cps_theory))-0.00005, max(max(cps), max(cps_theory))+0.00005)
    axs[3].set_ylim(min(cts)-0.001, max(cts)+0.001)
    axs[3].legend()
    # plt.show()

    opt_index = powers.index(min(powers))
    collective = cols[opt_index]
    inflow = inflows[opt_index]
    omega = omegas[opt_index]
    twist = twists[opt_index]
    # dct = dcts[opt_index]

    aoa_r = twist+collective-np.arctan(inflow/rotor.r)
    lift = (rotor.airfoil.cl_alpha*aoa_r*0.5*rotor.rho*((inflow*omega*rotor.R)**2+omega**2*(rotor.r*rotor.R)**2)*rotor.dr*rotor.c)

    plt.figure(1)
    fig2, ax = plt.subplots(4,1)

    ax[0].set_title('Forward flight. Optimum twist ratio: '+ str(t_index[opt_index]))
    ax[0].plot(rotor.r, inflow)
    ax[0].set_ylabel('inflow')

    ax[1].plot(rotor.r, aoa_r*180/pi)
    ax[1].set_ylabel('angle of attack')

    ax[2].plot(rotor.r, twist*180/pi, label= 'twist')
    ax[2].plot(rotor.r, [collective*180/pi]*np.int(len(rotor.r)), label= 'collective')
    ax[2].plot(rotor.r, np.arctan(inflow/rotor.r)*180/pi, label='phi')
    ax[2].legend()
    ax[2].set_ylabel('dct')

    ax[3].plot(rotor.r, twist*180/pi)
    ax[3].set_ylabel('twist')
    return powers
    # plt.show()

def combinedFlight(h_h, h_f, h_c, h_cl, R, MTOM, L_D_cr, L_D_cl, wing_chord):
    MTOW = MTOM*9.81
    T_f = MTOW/L_D_cr
    T_cl = MTOW/L_D_cl
    N, B, R, cutout, solidity, theta_tip, taper, airfoil, be = initVariables(R)
    powers_h, omegas_h, thrusts_h, t_ratios_h, dcts_h, cts_h, cps_h, inflows_h, cols_h, twists_h, t_index_h = genLists()
    powers_f, omegas_f, thrusts_f, t_ratios_f, dcts_f, cts_f, cps_f, inflows_f, cols_f, twists_f, t_index_f = genLists()
    powers_c, omegas_c, thrusts_c, t_ratios_c, dcts_c, cts_c, cps_c, inflows_c, cols_c, twists_c, t_index_c = genLists()
    powers_cl, omegas_cl, thrusts_cl, t_ratios_cl, dcts_cl, cts_cl, cps_cl, inflows_cl, cols_cl, twists_cl, t_index_cl = genLists()

    h_h = h_h
    h_f = h_f
    h_c = h_c
    h_cl = h_cl

    v_inf_h = 0
    v_inf_f = 100
    v_inf_c = 8
    v_inf_cl = 60


    rotor_h = Rotor(N, B, R, cutout, taper, solidity, theta_tip, airfoil, v_inf_h, be, h_h)
    rotor_f = Rotor(N, B, R, cutout, taper, solidity, theta_tip, airfoil, v_inf_f, be, h_f)
    rotor_c = Rotor(N, B, R, cutout, taper, solidity, theta_tip, airfoil, v_inf_c, be, h_c)
    rotor_cl = Rotor(N, B, R, cutout, taper, solidity, theta_tip, airfoil, v_inf_cl, be, h_cl)

    differenttwists, twist_index = genTwist(theta_tip, cutout, be)
    # theta_root = 20*pi/180
    # theta_tip  = 1*pi/180
    # slope = (theta_tip-theta_root)/(self.R-self.R0-self.dr*self.R)
    # slope = (theta_tip-theta_root)/(1-cutout)
    # one_twist = theta_root + slope*(rotor_h.r-rotor_h.cutout)
    #
    # one_twist = [one_twist]
    for t_i, twist in enumerate(differenttwists):
    # for t_i, twist in enumerate(one_twist):
        print(t_i+1, ' / ', len(differenttwists))

        power_i_h, thrust_i_h, omega_i_h, ct_i_h, cp_i_h, col_i_h, inflow_i_h, index_h  = genLists_i()
        power_i_f, thrust_i_f, omega_i_f, ct_i_f, cp_i_f, col_i_f, inflow_i_f, index_f  = genLists_i()
        power_i_c, thrust_i_c, omega_i_c, ct_i_c, cp_i_c, col_i_c, inflow_i_c, index_c  = genLists_i()
        power_i_cl, thrust_i_cl, omega_i_cl, ct_i_cl, cp_i_cl, col_i_cl, inflow_i_cl, index_cl  = genLists_i()

        rotor_h.twist = twist
        rotor_f.twist = twist
        rotor_c.twist = twist
        rotor_cl.twist = twist
        # plt.plot(rotor.r, twist)
        # plt.show()

        # for omega in np.linspace(20, min(rotor_h.maxOmega, rotor_f.maxOmega), 10):
        for omega in [rotor_f.maxOmega]:
            print(omega)
           # Objective ct
            goal_h = (MTOW+0.3*MTOW*wing_chord*N/(pi*rotor_h.R))/N/(rotor_h.rho*pi*rotor_h.R**2*(omega*rotor_h.R)**2)
            goal_f = T_f/N/(rotor_f.rho*pi*rotor_f.R**2*(omega*rotor_f.R)**2)
            goal_c = (MTOW+0.3*MTOW*wing_chord*N/(pi*rotor_c.R))/N/(rotor_h.rho*pi*rotor_c.R**2*(omega*rotor_c.R)**2)
            goal_cl = T_cl/N/(rotor_cl.rho*pi*rotor_cl.R**2*(omega*rotor_cl.R)**2)
            # print('goal_f')
            # print(goal_f)
            # print('goal_cl')
            # print(goal_cl)
            # Simulation
            # print('hover')
            power_h, thrust_h, omega_h, ct_h, cp_h, inflow_h, collective_h = flight(rotor_h, goal_h, omega)
            # print('forward flight')
            power_f, thrust_f, omega_f, ct_f, cp_f, inflow_f, collective_f = flight(rotor_f, goal_f, omega)
            # print('vertical climb')
            power_c, thrust_c, omega_c, ct_c, cp_c, inflow_c, collective_c = flight(rotor_c, goal_c, omega)
            # print('horizontal climb')
            power_cl, thrust_cl, omega_cl, ct_cl, cp_cl, inflow_cl, collective_cl = flight(rotor_cl, goal_cl, omega)

            # Append values for hover
            power_i_h.append(power_h)
            thrust_i_h.append(thrust_h)
            omega_i_h.append(omega_h)
            ct_i_h.append(ct_h)
            cp_i_h.append(cp_h)
            inflow_i_h.append(inflow_h)
            col_i_h.append(collective_h)
            # Append values for forward flight
            power_i_f.append(power_f)
            thrust_i_f.append(thrust_f)
            omega_i_f.append(omega_f)
            ct_i_f.append(ct_f)
            cp_i_f.append(cp_f)
            inflow_i_f.append(inflow_f)
            col_i_f.append(collective_f)
            # Append values for vertical climb
            power_i_c.append(power_c)
            thrust_i_c.append(thrust_c)
            omega_i_c.append(omega_c)
            ct_i_c.append(ct_c)
            cp_i_c.append(cp_c)
            inflow_i_c.append(inflow_c)
            col_i_c.append(collective_c)
            # Append values for horizontal climb
            power_i_cl.append(power_cl)
            thrust_i_cl.append(thrust_cl)
            omega_i_cl.append(omega_cl)
            ct_i_cl.append(ct_cl)
            cp_i_cl.append(cp_cl)
            inflow_i_cl.append(inflow_cl)
            col_i_cl.append(collective_cl)
        # print('cp_i_h')
        # print(cp_i_h)



        opt_list_h =np.array(np.array(cp_i_h)/np.array(ct_i_h))
        opt_list_f = np.array(np.array(cp_i_f)/np.array(ct_i_f))
        opt_list_c = np.array(np.array(cp_i_c)/np.array(ct_i_c))
        opt_list_cl = np.array(np.array(cp_i_cl)/np.array(ct_i_cl))
        opt_val_h = np.min(opt_list_h)
        opt_val_f = np.min(opt_list_f)
        opt_val_c = np.min(opt_list_c)
        opt_val_cl = np.min(opt_list_cl)


        if opt_val_h != 10**10 and not np.isnan(opt_val_h) and opt_val_f != 10**10 and not np.isnan(opt_val_f) and opt_val_c != 10**10 and not np.isnan(opt_val_c) and opt_val_cl != 10**10 and not np.isnan(opt_val_cl):
            # Find opt index
            opt_index_h = list(opt_list_h).index(opt_val_h)
            opt_index_f = list(opt_list_f).index(opt_val_f)
            opt_index_c = list(opt_list_c).index(opt_val_c)
            opt_index_cl = list(opt_list_cl).index(opt_val_cl)


            powers_h.append(power_i_h[opt_index_h])
            thrusts_h.append(thrust_i_h[opt_index_h])
            omegas_h.append(omega_i_h[opt_index_h])
            # t_ratios.append(t_ratio)
            cts_h.append(ct_i_h[opt_index_h])
            cps_h.append(cp_i_h[opt_index_h])
            twists_h.append(twist)
            inflows_h.append(inflow_i_h[opt_index_h])
            cols_h.append(col_i_h[opt_index_h])

            powers_f.append(power_i_f[opt_index_f])
            thrusts_f.append(thrust_i_f[opt_index_f])
            omegas_f.append(omega_i_f[opt_index_f])
            # t_ratios.append(t_ratio)
            cts_f.append(ct_i_f[opt_index_f])
            cps_f.append(cp_i_f[opt_index_f])
            twists_f.append(twist)
            inflows_f.append(inflow_i_f[opt_index_f])
            cols_f.append(col_i_f[opt_index_f])

            powers_c.append(power_i_c[opt_index_c])
            thrusts_c.append(thrust_i_c[opt_index_c])
            omegas_c.append(omega_i_c[opt_index_c])
            # t_ratios.append(t_ratio)
            cts_c.append(ct_i_c[opt_index_c])
            cps_c.append(cp_i_c[opt_index_c])
            twists_c.append(twist)
            inflows_c.append(inflow_i_c[opt_index_c])
            cols_c.append(col_i_c[opt_index_c])

            powers_cl.append(power_i_cl[opt_index_cl])
            thrusts_cl.append(thrust_i_cl[opt_index_cl])
            omegas_cl.append(omega_i_cl[opt_index_cl])
            # t_ratios.append(t_ratio)
            cts_cl.append(ct_i_cl[opt_index_cl])
            cps_cl.append(cp_i_cl[opt_index_cl])
            twists_cl.append(twist)
            inflows_cl.append(inflow_i_cl[opt_index_cl])
            cols_cl.append(col_i_cl[opt_index_cl])

            aoa_r_h = twist+col_i_h[opt_index_h]-inflow_h/rotor_h.r
            aoa_r_f = twist+col_i_f[opt_index_f]-inflow_f/rotor_f.r
            aoa_r_c = twist+col_i_c[opt_index_c]-inflow_c/rotor_c.r
            aoa_r_cl = twist+col_i_cl[opt_index_cl]-inflow_cl/rotor_cl.r

            # lift = B*np.sum((rotor.airfoil.cl_alpha*aoa_r*0.5*rotor.rho*((inflow*omega*rotor.R)**2+omega**2*(rotor.r*rotor.R)**2)*rotor.dr*rotor.c))
            # thrust_from_ct = ct*0.5*rotor.rho*(omega*rotor.R)**2*pi*rotor.R**2
            t_index_h.append(twist_index[t_i])
            t_index_f.append(twist_index[t_i])
            t_index_c.append(twist_index[t_i])
            t_index_cl.append(twist_index[t_i])
            # print('lift', lift)
            # print('thrust', thrust_from_ct)
        else:
            print('Non valid twist ratio')

    # Thrust corrections for number of rotors
    thrust_h = list(np.array(thrusts_h)*rotor_h.N)
    thrust_f = list(np.array(thrusts_f)*rotor_f.N)
    thrust_c = list(np.array(thrusts_c)*rotor_c.N)
    thrust_cl = list(np.array(thrusts_cl)*rotor_cl.N)


    # Power corrections for number of rotors and climb
    powers_h = list(np.array(powers_h)*rotor_h.N)
    powers_f = list(np.array(powers_f)*rotor_f.N)
    powers_c = list(np.array(powers_c)*rotor_c.N + MTOW*8/2)
    powers_cl =list(np.array(powers_cl)*rotor_cl.N + MTOW*8)


    plot_hover = True
    plot_cruise = True
    plot_v_climb = True
    plot_h_climb = True
    write = True

    if plot_hover:
        '''
        Plots for hover
        '''
        fig, axsh = plt.subplots(4,1)

        axsh[0].set_title('Hover')
        axsh[0].plot(t_index_h, powers_h)
        axsh[0].set_ylabel('power [W]')
        # plt.show()

        axsh[1].plot(t_index_h, thrusts_h)
        axsh[1].set_ylabel('thrust [N]')
        # plt.show()

        axsh[2].plot(t_index_h, omegas_h)
        axsh[2].set_ylabel('omega [rad/s]')


        k = 1.15
        cd0 = 0.011
        cps_theory_h = k*np.array(cts_h)**(3/2)/(np.sqrt(2))+(solidity*cd0/8)
        # print('theory', cps_theory)
        ratios_h = np.array(cts_h)/np.array(cps_h)

        axsh[3].scatter(cps_h, cts_h, label='single')
        axsh[3].plot(cps_theory_h, cts_h, label='theory')
        axsh[3].scatter(cps_h[np.where(ratios_h == np.max(ratios_h))[0][0]], cts_h[np.where(ratios_h == np.max(ratios_h))[0][0]], label='Opt point')
        axsh[3].set_ylabel('ct to cp')
        axsh[3].set_xlim(min(min(cps_h), min(cps_theory_h))-0.00005, max(max(cps_h), max(cps_theory_h))+0.00005)
        axsh[3].set_ylim(min(cts_h)-0.001, max(cts_h)+0.001)
        axsh[3].legend()
        # plt.show()

        opt_index_h = powers_h.index(min(powers_h))
        collective_h = cols_h[opt_index_h]
        inflow_h = inflows_h[opt_index_h]
        omega_h = omegas_h[opt_index_h]
        twist_h = twists_h[opt_index_h]
        # dct = dcts[opt_index]

        aoa_r_h = twist_h+collective_h-inflow_h/rotor_h.r
        # lift = (rotor_h.airfoil.cl_alpha*aoa_r*0.5*rotor_h.rho*((inflow_h*omega_h*rotor_h.R)**2+omega_h**2*(rotor_h.r*rotor_h.R)**2)*rotor_h.dr*rotor_h.c)

        # plt.figure(1)
        fig2, axh = plt.subplots(4,1)

        axh[0].set_title('Hover.Optimum twist ratio: '+ str(t_index_h[opt_index_h]))
        axh[0].plot(rotor_h.r, inflow_h)
        axh[0].set_ylabel('inflow')

        axh[1].plot(rotor_h.r, aoa_r_h*pi/180)
        axh[1].set_ylabel('angle of attack')

        axh[2].plot(rotor_h.r, twist_h*180/pi, label= 'twist')
        axh[2].plot(rotor_h.r, [collective_h*180/pi]*np.int(len(rotor_h.r)), label= 'collective')
        axh[2].plot(rotor_h.r, inflow_h/rotor_h.r*180/pi, label='phi')
        axh[2].legend()
        axh[2].set_ylabel('dct')

        axh[3].plot(rotor_h.r, twist*180/pi)
        axh[3].set_ylabel('twist')
        # plt.show()

    if plot_cruise:
        '''
        Plots for forward flight
        '''
        fig, axsf = plt.subplots(4,1)

        axsf[0].set_title('Forward flight')
        axsf[0].plot(t_index_f, powers_f)
        axsf[0].set_ylabel('power [W]')
        # plt.show()

        axsf[1].plot(t_index_f, thrusts_f)
        axsf[1].set_ylabel('thrust [N]')
        # plt.show()

        axsf[2].plot(t_index_f, omegas_f)
        axsf[2].set_ylabel('omega [rad/s]')


        k = 1.15
        cd0 = 0.011
        cps_theory_f = k*np.array(cts_f)**(3/2)/(np.sqrt(2))+(solidity*cd0/8)
        # print('theory', cps_theory)
        ratios_f = np.array(cts_f)/np.array(cps_f)

        axsf[3].scatter(cps_f, cts_f, label='single')
        axsf[3].plot(cps_theory_f, cts_f, label='theory')
        axsf[3].scatter(cps_f[np.where(ratios_f == np.max(ratios_f))[0][0]], cts_f[np.where(ratios_f == np.max(ratios_f))[0][0]], label='Opt point')
        axsf[3].set_ylabel('ct to cp')
        axsf[3].set_xlim(min(min(cps_f), min(cps_theory_f))-0.00005, max(max(cps_f), max(cps_theory_f))+0.00005)
        axsf[3].set_ylim(min(cts_f)-0.001, max(cts_f)+0.001)
        axsf[3].legend()
        # plt.show()

        opt_index_f = powers_f.index(min(powers_f))
        collective_f = cols_f[opt_index_f]
        inflow_f = inflows_f[opt_index_f]
        omega_f = omegas_f[opt_index_f]
        twist_f = twists_f[opt_index_f]
        # dct = dcts[opt_index]

        aoa_r_f = twist_f+collective_f-inflow_f/rotor_f.r
        # lift = (rotor_f.airfoil.cl_alpha*aoa_r_f*0.5*rotor_f.rho*((inflow_f*omega*rotor_f.R)**2+omega_f**2*(rotor_f.r*rotor_f.R)**2)*rotor_f.dr*rotor_f.c)

        # plt.figure(1)
        fig2, axf = plt.subplots(4,1)

        axf[0].set_title('Forward flight. Optimum twist ratio: '+ str(t_index_f[opt_index_f]))
        axf[0].plot(rotor_f.r, inflow_f)
        axf[0].set_ylabel('inflow')

        axf[1].plot(rotor_f.r, aoa_r_f*pi/180)
        axf[1].set_ylabel('angle of attack')

        axf[2].plot(rotor_f.r, twist_f*180/pi, label= 'twist')
        axf[2].plot(rotor_f.r, [collective_f*180/pi]*np.int(len(rotor_f.r)), label= 'collective')
        axf[2].plot(rotor_f.r, inflow_f/rotor_f.r*180/pi, label='phi')
        axf[2].legend()
        axf[2].set_ylabel('dct')

        axf[3].plot(rotor_f.r, twist_f*180/pi)
        axf[3].set_ylabel('twist')
        # plt.show()

    if plot_v_climb:
        '''
        Plots for vertical climb
        '''
        fig, axsc = plt.subplots(4,1)

        axsc[0].set_title('Vertical climb')
        axsc[0].plot(t_index_c, powers_c)
        axsc[0].set_ylabel('power [W]')
        # plt.show()

        axsc[1].plot(t_index_c, thrusts_c)
        axsc[1].set_ylabel('thrust [N]')
        # plt.show()

        axsc[2].plot(t_index_c, omegas_c)
        axsc[2].set_ylabel('omega [rad/s]')


        k = 1.15
        cd0 = 0.011
        cps_theory_c = k*np.array(cts_c)**(3/2)/(np.sqrt(2))+(solidity*cd0/8)
        # print('theory', cps_theory)
        ratios_c = np.array(cts_c)/np.array(cps_c)

        axsc[3].scatter(cps_c, cts_c, label='single')
        axsc[3].plot(cps_theory_c, cts_c, label='theory')
        axsc[3].scatter(cps_c[np.where(ratios_c == np.max(ratios_c))[0][0]], cts_c[np.where(ratios_c == np.max(ratios_c))[0][0]], label='Opt point')
        axsc[3].set_ylabel('ct to cp')
        axsc[3].set_xlim(min(min(cps_c), min(cps_theory_c))-0.00005, max(max(cps_c), max(cps_theory_c))+0.00005)
        axsc[3].set_ylim(min(cts_c)-0.001, max(cts_c)+0.001)
        axsc[3].legend()
        # plt.show()

        opt_index_c = powers_c.index(min(powers_c))
        collective_c = cols_c[opt_index_c]
        inflow_c = inflows_c[opt_index_c]
        omega_c = omegas_c[opt_index_c]
        twist_c = twists_c[opt_index_c]
        # dct = dcts[opt_index]

        aoa_r_c = twist_c+collective_c-inflow_c/rotor_c.r
        # lift = (rotor_c.airfoil.cl_alpha*aoa_r_c*0.5*rotor_c.rho*((inflow_c*omega_c*rotor_c.R)**2+omega_c**2*(rotor_c.r*rotor_c.R)**2)*rotor_c.dr*rotor_c.c)

        # plt.figure(1)
        fig2, axc = plt.subplots(4,1)

        axc[0].set_title('Vertical climb. Optimum twist ratio: '+ str(t_index_c[opt_index_c]))
        axc[0].plot(rotor_c.r, inflow_c)
        axc[0].set_ylabel('inflow')

        axc[1].plot(rotor_c.r, aoa_r_c*pi/180)
        axc[1].set_ylabel('angle of attack')

        axc[2].plot(rotor_c.r, twist_c*180/pi, label= 'twist')
        axc[2].plot(rotor_c.r, [collective_c*180/pi]*np.int(len(rotor_c.r)), label= 'collective')
        axc[2].plot(rotor_c.r, inflow_c/rotor_c.r*180/pi, label='phi')
        axc[2].legend()
        axc[2].set_ylabel('dct')

        axc[3].plot(rotor_c.r, twist_c*180/pi)
        axc[3].set_ylabel('twist')
        # plt.show()

    if plot_h_climb:
        '''
        Plots for horizontal climb
        '''
        fig, axscl = plt.subplots(4,1)

        axscl[0].set_title('Horizontal climb')
        axscl[0].plot(t_index_cl, powers_cl)
        axscl[0].set_ylabel('power [W]')
        # plt.show()

        axscl[1].plot(t_index_cl, thrusts_cl)
        axscl[1].set_ylabel('thrust [N]')
        # plt.show()

        axscl[2].plot(t_index_cl, omegas_cl)
        axscl[2].set_ylabel('omega [rad/s]')


        k = 1.15
        cd0 = 0.011
        cps_theory_cl = k*np.array(cts_cl)**(3/2)/(np.sqrt(2))+(solidity*cd0/8)
        # print('theory', cps_theory)
        ratios_cl = np.array(cts_cl)/np.array(cps_cl)

        axscl[3].scatter(cps_cl, cts_cl, label='single')
        axscl[3].plot(cps_theory_cl, cts_cl, label='theory')
        axscl[3].scatter(cps_cl[np.where(ratios_cl == np.max(ratios_cl))[0][0]], cts_cl[np.where(ratios_cl == np.max(ratios_cl))[0][0]], label='Opt point')
        axscl[3].set_ylabel('ct to cp')
        axscl[3].set_xlim(min(min(cps_cl), min(cps_theory_cl))-0.00005, max(max(cps_cl), max(cps_theory_cl))+0.00005)
        axscl[3].set_ylim(min(cts_cl)-0.001, max(cts_cl)+0.001)
        axscl[3].legend()
        # plt.show()

        opt_index_cl = powers_cl.index(min(powers_cl))
        collective_cl = cols_cl[opt_index_cl]
        inflow_cl = inflows_cl[opt_index_cl]
        omega_cl = omegas_cl[opt_index_cl]
        twist_cl = twists_cl[opt_index_cl]
        # dct = dcts[opt_index]

        aoa_r_cl = twist_cl+collective_cl-inflow_cl/rotor_cl.r
        # lift = (rotor_c.airfoil.cl_alpha*aoa_r_c*0.5*rotor_c.rho*((inflow_c*omega_c*rotor_c.R)**2+omega_c**2*(rotor_c.r*rotor_c.R)**2)*rotor_c.dr*rotor_c.c)

        # plt.figure(1)
        fig2, axcl = plt.subplots(4,1)

        axcl[0].set_title('Horizontal climb. Optimum twist ratio: '+ str(t_index_cl[opt_index_cl]))
        axcl[0].plot(rotor_cl.r, inflow_cl)
        axcl[0].set_ylabel('inflow')

        axcl[1].plot(rotor_cl.r, aoa_r_cl*pi/180)
        axcl[1].set_ylabel('angle of attack')

        axcl[2].plot(rotor_cl.r, twist_cl*180/pi, label= 'twist')
        axcl[2].plot(rotor_cl.r, [collective_cl*180/pi]*np.int(len(rotor_cl.r)), label= 'collective')
        axcl[2].plot(rotor_cl.r, inflow_cl/rotor_cl.r*180/pi, label='phi')
        axcl[2].legend()
        axcl[2].set_ylabel('dct')

        axcl[3].plot(rotor_cl.r, twist_cl*180/pi)
        axcl[3].set_ylabel('twist')
        # plt.show()

    '''
    Energy and power plots
    '''
    print('powers_h')
    print(powers_h)
    print('powers_f')
    print(powers_f)
    print('powers_vc')
    print(powers_c)
    print('powers_hcl')
    print(powers_cl)

    mission_power, time, TOe, Tendcruise = calcEnergy(powers_h, powers_f, powers_c, powers_cl)

    # fig, img = plt.subplots(2,1)
    # img[0].plot(t_index_h, energy)
    # img[0].set_ylabel('energy')
    # img[1].plot(t_index_h, powers_c)
    # img[1].set_ylabel('power')
    print(mission_power)
    print(time)
    print(TOe)
    print(Tendcruise)

    # powers_h = [powers_h]
    # powers_f = [powers_f]
    # powers_c = [powers_c]
    # powers_cl = [powers_cl]
    #
    # powers = np.append(powers_h, powers_f, axis=0)
    # powers = np.append(powers, powers_c, axis=0)
    # powers = np.append(powers, powers_cl, axis=0)
    if write:
        powers = pd.DataFrame(mission_power, columns=t_index_h)


        pre = os.path.dirname(os.path.realpath(__file__))
        fname = 'mission_profiles.xlsx'
        path = os.path.join(pre, fname)
        writer = pd.ExcelWriter(path, engine='xlsxwriter')
        # for i, df in enumerate(dfs):
        # for pow in range(mission_power):
        powers.to_excel(writer, sheet_name='power', index=False)

        writer.save()
    return mission_power, time, TOe, Tendcruise

def calcEnergy(p_hover, p_cruise, p_climb_v, p_climb_h):
    t_ho    = 60
    t_cl_v  = 30
    t_cl_h  = 470
    t_cr    = 2500
    t_de    = 500

    TOe = t_ho+t_cl_v+t_cl_h
    Tendcruise = TOe+t_cr+t_de

    time = np.arange(0, t_ho+t_cl_v+t_cl_h+t_cr+t_de+t_ho+0.1, 0.1)

    p_hover_lst = [p_hover]
    p_cruise_lst = [p_cruise]
    p_climb_v_lst = [p_climb_v]
    p_climb_h_lst = [p_climb_h]
    p_descent_lst = [list(2*np.array(p_cruise) -np.array(p_climb_h))]

    p_hover_lst = np.array(p_hover_lst)
    p_cruise_lst = np.array(p_cruise_lst)
    p_climb_v_lst = np.array(p_climb_v_lst)
    p_climb_h_lst = np.array(p_climb_h_lst)
    p_descent_lst = np.array(p_descent_lst)

    mission_power = np.append(np.repeat(p_hover_lst, t_ho*10, axis=0), np.repeat(p_climb_v_lst, t_cl_v*10, axis=0), axis=0) # intial hover and vertical climb
    mission_power = np.append(mission_power, np.repeat(p_climb_h_lst, t_cl_h*10, axis=0), axis=0)                           # horizontal climb
    mission_power = np.append(mission_power, np.repeat(p_cruise_lst, t_cr*10, axis=0), axis=0)                              # cruise
    mission_power = np.append(mission_power, np.repeat(p_descent_lst, t_de*10, axis=0), axis=0)                             # descent
    mission_power = np.append(mission_power, np.repeat(p_hover_lst, t_ho*10, axis=0), axis=0)                               # final hover



    # energy = 2*np.array(p_hover)*t_ho + np.array(p_climb)*t_cl + np.array(p_cruise)*t_cr + (np.array(p_cruise)-(np.array(p_climb)-np.array(p_cruise)))*t_de
    return mission_power, time, TOe, Tendcruise

# def missionFlight(h_h, h_c, h_f, R, MTOM, L_D, wing_chord, v_inf):
#     MTOW = MTOM*9.81
#     N, B, R, cutout, solidity, theta_tip, taper, airfoil, be = initVariables(R)
#     powers_h, omegas_h, thrusts_h, t_ratios_h, dcts_h, cts_h, cps_h, inflows_h, cols_h, twists_h, t_index_h = genLists()
#     powers_h, omegas_h, thrusts_h, t_ratios_h, dcts_h, cts_h, cps_h, inflows_h, cols_h, twists_h, t_index_h = genLists()
#     powers_c, omegas_c, thrusts_c, t_ratios_c, dcts_c, cts_c, cps_c, inflows_c, cols_c, twists_c, t_index_c = genLists()
#
#     differenttwists, twist_index = genTwist(theta_tip, cutout, be)
#     return

def flight(rotor, goal, omega):
            collective = 0
            thrustmatch = False
            iterations = 0
            while not thrustmatch:
                iterations += 1
                thrust, power, ct, cp, inflow, dct, dcp = rotor.simulation(omega, collective)
                if ct == 1000:
                    # print('sim not valid')
                    power_i = 10**10
                    thrust_i = 10**10
                    omega_i = 10**10
                    ct_i = 10**10
                    cp_i = 10**10
                    inflow_i = 10**10
                    col_i = 10**10
                    thrustmatch = True
                    dct_i = 10**10
                    dcp_i = 10**10

                else:
                    diff = ct - goal
                    # if iterations > 100:
                    #     print(iterations)
                    # print(diff)
                    # if diff > 0.0001:
                    #     collective -= 3 * diff
                    if diff > 0.00001:
                        # collective -= 40 /iterations * diff
                        collective -= 10 * diff
                    # elif diff < -0.001:
                    #     collective -= 3 * diff
                    elif diff < -0.00001:
                        # collective -= 40 / iterations * diff
                        collective -= 10 * diff
                    else:
                        thrustmatch = True
                        power_i = power
                        thrust_i = thrust
                        omega_i = omega
                        ct_i = ct
                        cp_i = cp
                        inflow_i = inflow
                        col_i = collective
                        dct_i = dct
                        dcp_i = dcp
            return power_i, thrust_i, omega_i, ct_i, cp_i, inflow_i, col_i, dct_i, dcp_i

def diskloading(h, R, T):
    DL = np.arange(30, 100, 1)
    MTOW = T
    N, B, R, cutout, solidity, theta_tip, taper, airfoil, be = initVariables(R)

def finalsim(initial_list, rotor, wing_chord):
    hs = []
    v_infs = []
    angles = []
    dcts = []
    dcps = []
    Machs = []
    p_outs = []
    powers, omegas, thrusts, t_ratios, dcts, cts, cps, inflows, cols, twists, t_index = genLists()

    bar = progressbar.ProgressBar(maxval=len(initial_list), widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
    bar.start()

    for status_idx, cond in enumerate(initial_list):
        MTOW = cond[0]
        rotor.v_inf = cond[1]
        rotor.h = cond[2]
        angle = cond[3]
        time = cond[4]
        # N, B, R, cutout, solidity, theta_tip, taper, airfoil, be = initVariables(R)
        # differenttwists, twist_index = genTwist(theta_tip, cutout, be)
        # rotor = Rotor(N, B, R, cutout, taper, solidity, theta_tip, airfoil, v_inf, be, h)
        dct_i = []
        dcp_i = []
        Mach_i = []
        power_i, thrust_i, omega_i, ct_i, cp_i, col_i, inflow_i, index  = genLists_i()
        # print('MTOW',MTOW)
        # print('conds:')
        # print('v_inf', rotor.v_inf)
        # print('h', rotor.h)
        # print('angle', angle)
        # print('time', time)

        for omega in np.linspace(25, rotor.maxOmega, 7):
            # print('omega', omega)
            goal = (MTOW+np.sin(angle)*0.3*MTOW*wing_chord*rotor.N/(pi*rotor.R))/rotor.N/(rotor.rho*pi*rotor.R**2*(omega*rotor.R)**2)
            power, thrust, omega, ct, cp, inflow, collective, dct, dcp = flight(rotor, goal, omega)

            power_i.append(power)
            thrust_i.append(thrust)
            omega_i.append(omega)
            ct_i.append(ct)
            cp_i.append(cp)
            inflow_i.append(inflow)
            col_i.append(collective)
            dct_i.append(dct)
            dcp_i.append(dcp)
            Mach_i.append(np.sqrt((omega*rotor.R)**2+rotor.v_inf**2)/rotor.v_sound)

        if min(cp_i) != 10**10:
            opt_list = np.array(cp_i)/np.array(ct_i)
            opt_index = list(opt_list).index(min(opt_list))

            # print(rotor.twist+collective-inflow/rotor.r)
            # if status_idx in [2, 3, 4, 5, 6]:
            #     power_i[opt_index] = power_i[opt_index]+MTOW*8/4
            # elif status_idx in [8, 9, 10]:
            #     power_i[opt_index] = power_i[opt_index]-MTOW*8/4
            powers.append(power_i[opt_index])
            thrusts.append(MTOW)
            omegas.append(omega_i[opt_index])
            t_index.append(time)
            cts.append(ct_i[opt_index])
            cps.append(cp_i[opt_index])
            inflows.append(inflow_i[opt_index])
            cols.append(col_i[opt_index])
            dcts.append(dct_i[opt_index])
            dcps.append(dcp_i[opt_index])
            Machs.append(Mach_i[opt_index])
            hs.append(rotor.h)
            v_infs.append(rotor.v_inf)
            angles.append(angle)
            # plt.plot(rotor.r, rotor.twist+col_i[opt_index]-inflow_i[opt_index]/rotor.r, label=str(status_idx))
            # plt.plot(rotor.r, len(rotor.r)*[col_i[opt_index]], label=str(status_idx))
            # plt.plot(rotor.r, rotor.twist+col_i[opt_index], label=str(status_idx))

        else:
            print('Non valid BET simulation.')
        bar.update(status_idx+1)

        # print('Thrust:', MTOW, '[N]. Power:', powers)
    # plt.plot(rotor.r, rotor.twist, label='twist')
    # plt.legend()
    # plt.show()
    bar.finish()
    powers = np.array(powers)*2/1000

    powers = editPowers(powers, initial_list[0][0]/9.81)


    # print('p_outs', p_outs)
    # print('Propeller efficiency')
    # print(p_eff)

    mission_profile, mission_time = calcMissionProfile(powers, t_index)
    # print(powers, t_index)
    # print(powers)

    # plt.plot(mission_time, mission_profile)
    # plt.show()

    # print('time')
    # print(t_index)

    # return powers
    return powers, t_index, hs, v_infs, angles, mission_profile, mission_time, dcts, omegas, cols, inflows, dcps, Machs, cts, cps

def calcBladeStructure(rotor, dcts, omegas, collectives, inflows):
    for idx, dct in enumerate(dcts):
        # Calculate moments along span due to thrust
        dT = np.array(dct)*pi*rotor.R**4*omegas[idx]**2
        M_x = [0]
        for dT_i in dT[::-1]:
            M_x.append(M_x[-1]+dT_i*rotor.dr)

        M_x = np.abs(M_x)

        # Calculate drag force
        theta = rotor.twist + collectives[idx]
        aoa_r = theta-inflows[idx]/rotor.r
        closest_aoa_r = np.round(aoa_r*180/pi *4)/4
        index = ((closest_aoa_r - rotor.airfoil.info[0, 0])*4-1).astype(int)
        cd = rotor.airfoil.info[index, 2]
        dD = 0.5*cd*1.225*omegas[idx]**2*rotor.r**3*rotor.R**3*rotor.c

        # Calculate cetripetal acceleration due to spin
        a_z = (1-rotor.r)*omegas[idx]**2

        t = rotor.c*rotor.airfoil.t_c/100
        y_max = t/2

        # ct = np.trapz(dct, rotor.r)+dct[-1]*rotor.dr

        # plt.plot(rotor.r, dT, label=str(id)+', ct:'+str(ct))

        # plt.plot(rotor.r, M_x[::-1][1:])
    # plt.plot(rotor.r, dT)
    # plt.legend()
    # plt.show()

def calcRotorWeight(rotor):
    W_b_i = 0.026*(rotor.B)**0.66*(rotor.croot+rotor.ctip)/2*3.28084*(rotor.R*3.28084)**1.3*(rotor.v_tip_max*3.28084)**0.67
    W_b = W_b_i*rotor.B*rotor.N

    W_h_i = 0.0037*(rotor.B)**0.28*(rotor.R*3.28085)**1.5*(rotor.v_tip_max*3.28084)**0.43*(0.67*W_b_i+(422.53*9.81*3.28084)/(rotor.R*3.28084)**2)*0.55
    W_h = W_h_i*rotor.N

    W_b_kg = W_b*0.453592*0.85
    W_h_kg = W_h*0.453592*0.85

    return W_b_kg, W_h_kg, W_b_kg+W_h_kg

def calcMissionProfile(powers, times):
    mission_profile = np.array([])
    time_new = np.array([])
    for i, time in enumerate(times[:-1]):
        x = [times[i], times[i+1]]
        y = [powers[i], powers[i+1]]

        coefficients = np.polyfit(x, y, 1)

        polynomial = np.poly1d(coefficients)
        x_new = np.arange(times[i], times[i+1], 0.1)
        mission_power_i = polynomial(x_new)
        time_new = np.append(time_new, x_new)
        mission_profile = np.append(mission_profile, mission_power_i)
    return mission_profile, time_new

def nonfinalsim(rotor, h, T, wing_chord, v_inf, angle):
    MTOW = T

    rotor.h = h
    rotor.v_inf = v_inf

    power_i, thrust_i, omega_i, ct_i, cp_i, col_i, inflow_i, index  = genLists_i()
    rotor.calcTwist('linear', 20*pi/180, 1*pi/180)

    for omega in np.linspace(50, rotor.maxOmega, 10):
        print(omega)
        goal = (MTOW+0.3*MTOW*wing_chord*rotor.N/(pi*rotor.R)*np.sin(angle))/rotor.N/(rotor.rho*pi*rotor.R**2*(omega*rotor.R)**2)

        power, thrust, omega, ct, cp, inflow, collective = flight(rotor, goal, omega)
        if power != 10**10:
            power_i.append(power)
            thrust_i.append(thrust)
            omega_i.append(omega)
            ct_i.append(ct)
            cp_i.append(cp)
            inflow_i.append(inflow)
            col_i.append(collective)


    # if min(cp_i) != 10**10:
    #     opt_list = np.array(cp_i)/np.array(ct_i)
    #     opt_index = list(opt_list).index(min(opt_list))
    #
    #     # print(rotor.twist+collective-inflow/rotor.r)
    #     powers.append(power_i[opt_index])
    #     thrusts.append(thrust_i[opt_index])
    #     omegas.append(omega_i[opt_index])
    #     # t_ratios.append(t_ratio)
    #     cts.append(ct_i[opt_index])
    #     cps.append(cp_i[opt_index])
    #     twists.append(rotor.twist)
    #     inflows.append(inflow_i[opt_index])
    #     cols.append(col_i[opt_index])
    #
    #     aoa_r = rotor.twist+col_i[opt_index]-inflow_i[opt_index]/rotor.r
    #
    #     lift = B*np.sum((rotor.airfoil.cl_alpha*aoa_r*0.5*rotor.rho*((inflow*omega*rotor.R)**2+omega**2*(rotor.r*rotor.R)**2)*rotor.dr*rotor.c))
    #     thrust_from_ct = ct*0.5*rotor.rho*(omega*rotor.R)**2*pi*rotor.R**2
    #     t_index.append(twist_index[t_i])
    #     # print('lift', lift)
    #     # print('thrust', thrust_from_ct)
    # else:
    #     print('Non valid twist ratio')

    fig, axs = plt.subplots(3,1)

    axs[0].set_title('Hover')
    axs[0].plot(omega_i, np.array(power_i)/1000*rotor.N)
    axs[0].set_ylabel('power [kW]')
    # plt.show()

    axs[1].plot(omega_i, np.array(thrust_i)*rotor.N)
    axs[1].set_ylabel('thrust [N]')
    # plt.show()

    k = 1.15
    cd0 = 0.011
    cps_theory = k*np.array(ct_i)**(3/2)/(np.sqrt(2))+(rotor.sigma*cd0/8)
    # print('theory', cps_theory)
    ratios = np.array(ct_i)/np.array(cp_i)

    axs[2].scatter(cp_i, ct_i, label='single')
    axs[2].plot(cps_theory, ct_i, label='theory')
    axs[2].scatter(cp_i[np.where(ratios == max(ratios))[0][0]], ct_i[np.where(ratios == max(ratios))[0][0]], label='Opt point')
    axs[2].set_ylabel('ct to cp')
    axs[2].set_xlim(min(min(cp_i), min(cps_theory))-0.00005, max(max(cp_i), max(cps_theory))+0.00005)
    axs[2].set_ylim(min(ct_i)-0.001, max(ct_i)+0.001)
    axs[2].legend()
    # plt.show()

    opt_index = power_i.index(min(power_i))
    collective = col_i[opt_index]
    inflow = inflow_i[opt_index]
    omega = omega_i[opt_index]
    # dct = dcts[opt_index]

    aoa_r = rotor.twist+collective-np.arctan(inflow/rotor.r)
    lift = (rotor.airfoil.cl_alpha*aoa_r*0.5*rotor.rho*((inflow*omega*rotor.R)**2+omega**2*(rotor.r*rotor.R)**2)*rotor.dr*rotor.c)

    plt.figure(1)
    fig2, ax = plt.subplots(4,1)

    # ax[0].set_title('Hover. Optimum twist ratio: '+ str(t_index[opt_index]))
    ax[0].plot(rotor.r, inflow)
    ax[0].set_ylabel('inflow')

    ax[1].plot(rotor.r, aoa_r*180/pi)
    ax[1].set_ylabel('angle of attack')

    ax[2].plot(rotor.r, rotor.twist*180/pi, label= 'twist')
    ax[2].plot(rotor.r, [collective*180/pi]*np.int(len(rotor.r)), label= 'collective')
    ax[2].plot(rotor.r, np.arctan(inflow/rotor.r)*180/pi, label='phi')
    ax[2].legend()
    ax[2].set_ylabel('dct')

    ax[3].plot(rotor.r, rotor.twist*180/pi)
    ax[3].set_ylabel('twist')
    # ax[3].plot(rotor.r, np.sqrt((rotor.r*omega)**2+(inflow*omega*rotor.R)**2))
    # ax[3].plot(cp, ct)
    # plt.show()
    min_power = min(power_i)
    # print(power_i)
    # print(min(power_i))
    # print(np.where(power_i == min_power))
    # print(omega_i[np.where(power_i == min_power)[0][0]])
    return 2*min_power, omega_i[np.where(power_i == min_power)[0][0]]

def simplesim(R, cond):
    R = 4.44
    MTOW = cond[0]
    v_inf = cond[1]
    h = cond[2]
    angle = cond[3]
    powers, omegas, thrusts, t_ratios, dcts, cts, cps, inflows, cols, twists, t_index = genLists()
    N, B, R, cutout, solidity, theta_tip, taper, airfoil, be = initVariables(R)
    # v_inf = 7.98770778188738
    # v_inf = 8
    rotor = Rotor(N, B, R, cutout, taper, solidity, theta_tip, airfoil, v_inf, be, h)
    # angle = 0.0018937530751135466
    # angle = pi/2
    rotor.calcTwist('linear',20*pi/180,0*pi/180)
    # thrust, power, ct, cp, inflow, idx = rotor.simulation2(500/12.5, 8*np.pi/180)
    wing_chord = 2

    dct_i = []
    solidities = []
    # print('MTOW',MTOW)
    # print('conds:')
    # print('v_inf', rotor.v_inf)
    # print('h', rotor.h)
    # print('angle', angle)
    # print('time', time)
    for rotsol in np.linspace(0.03, 0.1, 7):#[0.03, 0.04, 0.05, 0.06, 0.07, 0.08]:
    # for rotsol in np.linspace(0.06, 0.1, 15):
    # for rotsol in [0.065, 0.075]:
    # for rotsol in [0.04, 0.08]:
        rotor.sigma = rotsol
        power_i, thrust_i, omega_i, ct_i, cp_i, col_i, inflow_i, index  = genLists_i()
        for omega in np.linspace(15, rotor.maxOmega, 10):
        # for omega in [40,45,50,55,60]:
        # for omega in [rotor.maxOmega]:
            # print('omega', omega)
            goal = (MTOW+np.sin(angle)*0.3*MTOW*wing_chord*rotor.N/(pi*rotor.R))/rotor.N/(rotor.rho*pi*rotor.R**2*(omega*rotor.R)**2)
            power, thrust, omega, ct, cp, inflow, collective, dct = flight(rotor, goal, omega)

            power_i.append(power)
            thrust_i.append(thrust)
            omega_i.append(omega)
            ct_i.append(ct)
            cp_i.append(cp)
            inflow_i.append(inflow)
            col_i.append(collective)
            dct_i.append(dct)
        # print(omega_i)
        # print(power_i)
        # plt.plot(omega_i, power_i)
        # plt.show()
        if min(cp_i) != 10**10:
            opt_list = np.array(cp_i)#/np.array(ct_i)
            opt_index = list(opt_list).index(min(opt_list))
            powers.append(power_i[opt_index])
            thrusts.append(MTOW)
            omegas.append(omega_i[opt_index])
            t_index.append(time)
            cts.append(ct_i[opt_index])
            cps.append(cp_i[opt_index])
            inflows.append(inflow_i[opt_index])
            cols.append(col_i[opt_index])
            dcts.append(dct_i[opt_index])
            solidities.append(rotsol)
            # hs.append(rotor.h)
            # v_infs.append(rotor.v_inf)
            # angles.append(angle)
        else:
            print('Non valid BET')
    fig, ar = plt.subplots(2,1)
    # print(solidities)

    ar[0].plot(solidities, powers)
    ar[1].plot(solidities, omegas)
    idxs = argrelextrema(np.array(powers), np.less)[0]
    p_new = []
    sol_new = []
    for i in idxs:
        sol_new.append(solidities[i])
        p_new.append(powers[i])
    powers = p_new
    solidities = sol_new
    # print(powers)
    # print(powers)

    # fig, ar = plt.subplots(2,1)
    # ar[0].plot(solidities, omegas)
    # ar[1].plot(solidities, 2*np.array(powers)/1000+MTOW*8/2000)
    plt.show()


    opt_list = np.array(powers)
    opt_index = list(opt_list).index(min(opt_list))
    power = powers[opt_index]
    thrust = thrusts[opt_index]
    omega = omegas[opt_index]
    ct = cts[opt_index]
    cp = cps[opt_index]
    solidity = solidities[opt_index]
            # inflow_i[opt_index]
            # col_i[opt_index]
            # dct_i[opt_index]


    #
    # print('ct', ct)
    # print('cp', cp)
    # print('theoretical thrust', ct*rotor.rho*pi*rotor.R**2*(omega*rotor.R)**2/9.81*2)
    # print('theoretical power', cp*rotor.rho*pi*rotor.R**2*(omega*rotor.R)**3/1000*2)

    # goal = (MTOW+0.3*MTOW*wing_chord*rotor.N/(pi*rotor.R)*np.sin(angle))/rotor.N/(rotor.rho*pi*rotor.R**2*(omega*rotor.R)**2)
    # print('omega', omega)
    # power, thrust, omega, ct, cp, inflow, collective, dct = flight(rotor, goal, 10)
    # print('power raw', power/500)
    power = 2*power/1000+MTOW*8/2000
    # thrust, power, ct, cp, inflow, F = rotor.simulation(500/12.5, 8*np.pi/180)
    print('power', power)
    print('solidity', solidity)
    # plt.plot(rotor.r, F)
    # plt.figure(1)
    # plt.plot(rotor.r, inflow)
    # plt.plot()
    # plt.show()




def simplesim2():
    R = 4.44
    MTOW = 4000*9.81
    v_inf = 8
    h =0
    angle = pi/2
    powers, omegas, thrusts, t_ratios, dcts, cts, cps, inflows, cols, twists, t_index = genLists()
    N, B, R, cutout, solidity, theta_tip, taper, airfoil, be = initVariables(R)
    # v_inf = 7.98770778188738
    # v_inf = 8
    rotor = Rotor(N, B, R, cutout, taper, solidity, theta_tip, airfoil, v_inf, be, h)
    # angle = 0.0018937530751135466
    # angle = pi/2
    rotor.calcTwist('linear',20*pi/180,0*pi/180)
    # thrust, power, ct, cp, inflow, idx = rotor.simulation2(500/12.5, 8*np.pi/180)
    wing_chord = 2

    rotor.sigma = 0.02
    dct_i = []
    power_i, thrust_i, omega_i, ct_i, cp_i, col_i, inflow_i, index  = genLists_i()
    for omega in np.linspace(15, rotor.maxOmega, 10):
    # for omega in [40,45,50,55,60]:
    # for omega in [rotor.maxOmega]:
        # print('omega', omega)
        goal = (MTOW+np.sin(angle)*0.3*MTOW*wing_chord*rotor.N/(pi*rotor.R))/rotor.N/(rotor.rho*pi*rotor.R**2*(omega*rotor.R)**2)
        power, thrust, omega, ct, cp, inflow, collective, dct = flight(rotor, goal, omega)

        power_i.append(power)
        thrust_i.append(thrust)
        omega_i.append(omega)
        ct_i.append(ct)
        cp_i.append(cp)
        inflow_i.append(inflow)
        col_i.append(collective)
        dct_i.append(dct)
    # print(omega_i)
    # print(power_i)
    # plt.plot(omega_i, power_i)
    # plt.show()
    if min(cp_i) != 10**10:
        opt_list = np.array(cp_i)#/np.array(ct_i)
        opt_index = list(opt_list).index(min(opt_list))
        powers.append(power_i[opt_index])
        thrusts.append(MTOW)
        omegas.append(omega_i[opt_index])
        t_index.append(time)
        cts.append(ct_i[opt_index])
        cps.append(cp_i[opt_index])
        inflows.append(inflow_i[opt_index])
        cols.append(col_i[opt_index])
        dcts.append(dct_i[opt_index])
        # hs.append(rotor.h)
        # v_infs.append(rotor.v_inf)
        # angles.append(angle)
    else:
        print('Non valid BET')
    fig, ar = plt.subplots(2,1)
    # print(solidities)

    power_i = 2*np.array(power_i)/1000+MTOW*8/2000


    theta = rotor.twist + cols
    aoa_r = theta - inflows/rotor.r
    print(aoa_r[0])

    ar[0].plot(omega_i, col_i)
    ar[1].plot(rotor.r, aoa_r[0])
    # ar[1].plot(solidities, omegas)
    # idxs = argrelextrema(np.array(powers), np.less)[0]
    # p_new = []
    # sol_new = []
    # for i in idxs:
    #     sol_new.append(solidities[i])
    #     p_new.append(powers[i])
    # powers = p_new
    # solidities = sol_new
    # print(powers)
    # print(powers)

    # fig, ar = plt.subplots(2,1)
    # ar[0].plot(solidities, omegas)
    # ar[1].plot(solidities, 2*np.array(powers)/1000+MTOW*8/2000)
    plt.show()


    # opt_list = np.array(powers)
    # opt_index = list(opt_list).index(min(opt_list))
    # power = powers[opt_index]
    # thrust = thrusts[opt_index]
    # omega = omegas[opt_index]
    # ct = cts[opt_index]
    # cp = cps[opt_index]
    # solidity = solidities[opt_index]
            # inflow_i[opt_index]
            # col_i[opt_index]
            # dct_i[opt_index]


    #
    # print('ct', ct)
    # print('cp', cp)
    # print('theoretical thrust', ct*rotor.rho*pi*rotor.R**2*(omega*rotor.R)**2/9.81*2)
    # print('theoretical power', cp*rotor.rho*pi*rotor.R**2*(omega*rotor.R)**3/1000*2)

    # goal = (MTOW+0.3*MTOW*wing_chord*rotor.N/(pi*rotor.R)*np.sin(angle))/rotor.N/(rotor.rho*pi*rotor.R**2*(omega*rotor.R)**2)
    # print('omega', omega)
    # power, thrust, omega, ct, cp, inflow, collective, dct = flight(rotor, goal, 10)
    # print('power raw', power/500)
    power = 2*power/1000+MTOW*8/2000
    # thrust, power, ct, cp, inflow, F = rotor.simulation(500/12.5, 8*np.pi/180)
    print('power', power)
    print('solidity', solidity)
    # plt.plot(rotor.r, F)
    # plt.figure(1)
    # plt.plot(rotor.r, inflow)
    # plt.plot()
    # plt.show()





# def genMission(p_hover, p_climb, p_cruise)

# def tablegen()

def makeMP():
    N, B, R, cutout, solidity, theta_tip, taper, airfoil, be = initVariables(4.45)
    rotor = Rotor(N, B, R, cutout, taper, solidity, theta_tip, airfoil, 0, be, 0)
    twists, sols, t_index = genTwist(1*pi/180, rotor.cutout, rotor.be)
    # twists, t_index =

    # for i in np.linspace(1, 45, 20):
    #     theta_root = i*pi/180
    #     theta_tip = 0
    #     slope = (theta_tip-theta_root)/(1-cutout)
    #     twists.append(theta_root+slope*(r-cutout))

    pre = os.path.dirname(os.path.realpath(__file__))
    fname = 'mission_data.xlsx'
    path = os.path.join(pre, fname)
    data = pd.read_excel(path, sheet_name='mission_data').values

    # data = openPickle('mission_data.xlsx')
    missions = []
    for sol in sols:
        print(sol)
        for twist in twists:
            # rotor.calcTwist('linear',20*math.pi/180,math.pi/180)
            rotor.twist = twist
            rotor.sigma = sol
            powers, times, hs, v_infs, angles, mission_profile, mission_time, dcts, omegas, cols, inflows, dcps, Machs, cts, cps  = finalsim(data,rotor,2)
            missions.append(mission_profile)

    missions = np.array(missions).T
    write = True
    if write:
        powers = pd.DataFrame(missions)
        # time_t = pd.DataFrame(mission_time)


        pre = os.path.dirname(os.path.realpath(__file__))
        fname = 'mission_profiles_fin1.xlsx'
        path = os.path.join(pre, fname)
        # fnamet = 'mission_times2.xlsx'
        writer = pd.ExcelWriter(path, engine='xlsxwriter')
        # writert = pd.ExcelWriter(os.path.join(pre, fname), engine='xlsxwriter')
        # for i, df in enumerate(dfs):
        # for pow in range(mission_power):
        powers.to_excel(writer, sheet_name='power', index=False)
        # mission_time.to_excel

        writer.save()

def makeMP_i():
    N, B, R, cutout, solidity, theta_tip, taper, airfoil, be = initVariables(4.45)
    rotor = Rotor(N, B, R, cutout, taper, solidity, theta_tip, airfoil, 0, be, 0)
    twists, sols, t_index = genTwist(1*pi/180, rotor.cutout, rotor.be)
    # twists, t_index =

    # for i in np.linspace(1, 45, 20):
    #     theta_root = i*pi/180
    #     theta_tip = 0
    #     slope = (theta_tip-theta_root)/(1-cutout)
    #     twists.append(theta_root+slope*(r-cutout))

    pre = os.path.dirname(os.path.realpath(__file__))
    fname = 'mission_data.xlsx'
    path = os.path.join(pre, fname)
    data = pd.read_excel(path, sheet_name='mission_data').values

    # data = openPickle('mission_data.xlsx')

    r = np.linspace(0.2, 1, 100)[:-1]
    theta_root_1 = 1*pi/180
    theta_root_2 = 15*pi/180
    theta_root_3 = 30*pi/180
    theta_tip = 0
    slope1 = (theta_tip-theta_root_1)/(1-0.2)
    slope2 = (theta_tip-theta_root_2)/(1-0.2)
    slope3 = (theta_tip-theta_root_3)/(1-0.2)
    twists.append(theta_root_1+slope1*(r-0.2))



    missions = []
    combinations = [[theta_root_1+slope1*(r-0.2), 0.02],[theta_root_2+slope2*(r-0.2), 0.05], [theta_root_3+slope3*(r-0.2), 0.09]]
    for comb in combinations:
        # rotor.calcTwist('linear',20*math.pi/180,math.pi/180)
        rotor.twist = comb[0]
        rotor.sigma = comb[1]
        powers, times, hs, v_infs, angles, mission_profile, mission_time, dcts, omegas, cols, inflows, dcps, Machs, cts, cps  = finalsim(data,rotor,2)
        missions.append(mission_profile)

    missions = np.array(missions).T

    write = True
    if write:
        powers = pd.DataFrame(missions)
        # time_t = pd.DataFrame(mission_time)


        pre = os.path.dirname(os.path.realpath(__file__))
        fname = 'mission_profiles_i.xlsx'
        path = os.path.join(pre, fname)
        # fnamet = 'mission_times2.xlsx'
        writer = pd.ExcelWriter(path, engine='xlsxwriter')
        # writert = pd.ExcelWriter(os.path.join(pre, fname), engine='xlsxwriter')
        # for i, df in enumerate(dfs):
        # for pow in range(mission_power):
        powers.to_excel(writer, sheet_name='power', index=False)
        # mission_time.to_excel

        writer.save()

def calcPowers(rotor):
    power_table = {}
    bar = progressbar.ProgressBar(maxval=111, widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
    bar.start()

    for v_i in np.arange(0, 111, 1):
    # for v_i in [0]:
        # print(v_i+1, ' / ', 111)
        bar.update(v_i+1)
        power_v_i = {}
        for T in np.arange(0, 30000, 100):
        # for T in [1000]:
        #     print(T/100+1, ' / ', 210)
            if T < 3000:
                powers = []
                for omega in np.linspace(20,rotor.maxOmega, 5):
                    goal = T/rotor.N/(rotor.rho*pi*rotor.R**2*(omega*rotor.R)**2)
                    power, thrust, omega, ct, cp, inflow, collective, dcts = flight(rotor, goal, omega)
                    powers.append(power)
                power_v_i[T] = min(powers)
            else:
                goal = T/(rotor.rho*pi*rotor.R**2*(rotor.maxOmega*rotor.R)**2)
                power, thrust, omega, ct, cp, inflow, collective, dct = flight(rotor, goal, rotor.maxOmega)
                power_v_i[T] = power

        power_table[v_i] = power_v_i
    bar.finish()

    print(power_table[0])
    pre = os.path.dirname(os.path.realpath(__file__))
    fname = 'power_table.pickle'
    path = os.path.join(pre, fname)
    with open(path, 'wb') as handle:
        pickle.dump(power_table, handle, protocol=pickle.HIGHEST_PROTOCOL)

def calcPowersFast(rotor):
    print('Calculating power tables for rotor...')
    time.sleep(0.1)
    power_table = {}

    bar = progressbar.ProgressBar(maxval=111, widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
    bar.start()
    # print(bar)
    t1 = time.time()
    for v_i in np.arange(0, 111, 1):
    # for v_i in [0]:
        # print(v_i+1, ' / ', 111)
        bar.update(v_i+1)
        # print()
        power_v_i = {}
        for T in np.arange(0, 30000, 200):
        # for T in [20000]:
            # print(T/100+1, ' / ', 210)
            goal = T/(rotor.rho*pi*rotor.R**2*(rotor.maxOmega*rotor.R)**2)
            power, thrust, omega, ct, cp, inflow, collective = flight(rotor, goal, rotor.maxOmega)
            power_v_i[T] = power
        power_table[v_i] = power_v_i
    t2 = time.time()
    bar.finish()

    rotor.power_tables = power_table
    # print(power_table)
    # pre = os.path.dirname(os.path.realpath(__file__))
    # fname = 'power_table.pickle'
    # path = os.path.join(pre, fname)
    # with open(path, 'wb') as handle:
    #     pickle.dump(power_table, handle, protocol=pickle.HIGHEST_PROTOCOL)
    print('Power tables completed! Time:', t2-t1, 'sec')
    return

def ctcpCurve(rotor):
    cts = []
    cps = []
    v_is= []
    idxs= []
    for v_inf in np.arange(0, 110, 1):
    # for v_inf in [0, 1, 2]:
        print(v_inf+1, ' / ', 110)
        omega_i, ct_i, cp_i, idx_i = genNewLists_i()
        rotor.v_inf = v_inf
        for idx, col in enumerate(np.linspace(0, 20*pi/180, 20)):
        # for idx, col in enumerate([500/12.5]):
            # print(idx+1, ' / ', 1000)
            omega_j, ct_j, cp_j, _ = genNewLists_i()
            for omega in np.linspace(20, rotor.maxOmega, 10):
                thrust, power, ct, cp, inflow, dct = rotor.simulation(omega, col)
                omega_j.append(omega)
                ct_j.append(ct)
                cp_j.append(cp)

            trim_list = np.where(np.array(ct_j)>0)[0]
            if trim_list != []:
                omega_j = np.array(omega_j)[trim_list]
                ct_j = np.array(ct_j)[trim_list]
                cp_j = np.array(cp_j)[trim_list]

                opt_list = ct_j/cp_j
                opt_idx = np.argmax(opt_list == np.max(opt_list))
                ct_i.append(ct_j[opt_idx])
                cp_i.append(cp_j[opt_idx])
                omega_i.append(omega_j[opt_idx])
                idx_i.append(ct_j[opt_idx]*omega_j[opt_idx]**2)


            # trim_list = np.where(np.array(ct_j)<0)[0]
            # if len(trim_list) != 10:
            #     omega_j = np.array(omega_j)[trim_list]
            #     ct_j = np.array(ct_j)[trim_list]
            #     cp_j = np.array(cp_j)[trim_list]
            #
            #     opt_list = ct_j/cp_j
            #     opt_idx = np.argmax(opt_list == np.max(opt_list))
            #     ct_i.append(ct_j[opt_idx])
            #     cp_i.append(cp_j[opt_idx])
            #     omega_i.append(omega_j[opt_idx])
            #     idx_i.append(ct_j[opt_idx]*omega_j[opt_idx]**2)


        cts.append(ct_i)
        cps.append(cp_i)
        v_is.append(omega_i)
        idxs.append(idx_i)
    # print(idxs[0])
    # print(cts[0])
    # print()
    # print(idxs[1])
    # print()
    # print(idxs[2])

    # print(idxs)

    plot= True


    if plot:
        pre = os.path.dirname(os.path.realpath(__file__))
        fname = 'H1_single_fig2.csv'
        path = os.path.join(pre, fname)
        data = pd.read_csv(path).values
        plot_idx = 0
        pre2 = os.path.dirname(os.path.realpath(__file__))
        fname2 = 'ctcpvalidation.xlsx'
        path2 = os.path.join(pre2, fname2)
        data2 = pd.read_excel(path2).values


        # plt.plot(cps[0], cts[0])
        k = 1.15
        cd0 = 0.011
        cps_theory = k*np.array(cts[0])**(3/2)/(np.sqrt(2))+(rotor.sigma*cd0/8)
        # print('theory', cps_theory)
        # ratios = np.array(cts[0])/np.array(cps[0])


        plt.scatter((cps[plot_idx]), (cts[plot_idx]), label='single')
        plt.plot((cps_theory), (cts[plot_idx]), label='theory')
        plt.plot((data2[:,0])/10, (data2[:,1])*10, label='exp dat')
        plt.scatter((data[:,0]), (data[:,1]), label='experimental data')
        plt.legend()
        plt.show()



        plt.scatter(np.array(cps[plot_idx])/rotor.sigma, np.array(cts[plot_idx])/rotor.sigma, label='single')
        plt.plot(np.array(cps_theory)/rotor.sigma, np.array(cts[plot_idx])/rotor.sigma, label='theory')
        plt.plot(np.array(data2[:,0])/10/rotor.sigma, np.array(data2[:,1])*10/rotor.sigma, label='exp dat')
        plt.scatter(np.array(data[:,0])/rotor.sigma, np.array(data[:,1])/rotor.sigma, label='experimental data')
        # plt.scatter(cps[0][np.where(ratios == np.max(ratios))[0][0]], cts[0][np.where(ratios == np.max(ratios))[0][0]], label='Opt point')
        plt.ylabel('ct to cp')
        plt.xlim(min(min(cps[plot_idx]), min(cps_theory))-0.00005, max(max(cps[plot_idx]), max(cps_theory))+0.00005)
        plt.ylim(min(cts[plot_idx])-0.001, max(cts[plot_idx])+0.001)
        plt.legend()


        # goal = 0.01
        # omega = 80
        # collective = 0
        # thrustmatch = False
        # iterations = 0
        # while not thrustmatch:
        #     iterations += 1
        #     thrust, power, ct, cp, inflow = rotor.simulation(omega, collective)
        #     if ct == 1000:
        #         # print('sim not valid')
        #         power_i = 10**10
        #         thrust_i = 10**10
        #         omega_i = 10**10
        #         ct_i = 10**10
        #         cp_i = 10**10
        #         inflow_i = 10**10
        #         col_i = 10**10
        #         thrustmatch = True
        #
        #     else:
        #         diff = ct - goal
        #         if iterations > 100:
        #             print(iterations)
        #         # print(diff)
        #         # if diff > 0.0001:
        #         #     collective -= 3 * diff
        #         if diff > 0.00001:
        #             collective -= 100 /iterations * diff
        #         # elif diff < -0.001:
        #         #     collective -= 3 * diff
        #         elif diff < -0.00001:
        #             collective -= 100 / iterations * diff
        #         else:
        #             thrustmatch = True
        #             power_i = power
        #             thrust_i = thrust
        #             omega_i = omega
        #             ct_i = ct
        #             cp_i = cp
        #             inflow_i = inflow
        #             col_i = collective
        # return power_i, thrust_i, omega_i, ct_i, cp_i, inflow_i, col_i

def plotMass():
    fs = 20
    data = openPickle('masses.pickle')
    print(data)
    data['Mintots'] = np.array(data['Mintots'])# + np.array([0,0,0,0,0,0,0,0,0,0,0,3,4,5,5,5,5,5,5,15, 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,])
    # plt.plot(np.linspace(5, 35, 10), data['Mintots'])
    plt.plot(np.linspace(5, 35, 10), np.array(data['Mintots'])+11+45.4+301.529, linewidth=4, color= 'g')#+301.529
    plt.title('Mass of power plant for different rotor twists', fontsize=fs)
    # plt.grid()
    # plt.plot(data['Mintots'])
    plt.ylabel('Weight [kg]', fontsize=fs)
    # plt.xlabel('Blade Type', fontsize=fs)
    plt.xlabel('Blade Root Twist [$^{\circ}$]', fontsize=fs)
    plt.yticks(fontsize=fs)
    plt.xticks(fontsize=fs)
    idx = 0
    # for i in np.arange(5, 20, 5):
    #     plt.axvline(i, linestyle='--', c= 'k')
    #     line = r'$\sigma$ = '+str(0.01+idx)
    #     plt.text(i-3.3, 1400, line, fontsize=fs)
    #     idx += 0.01
    # line = r'$\sigma$ = '+str(0.01+idx)
    # plt.text(20-3.3, 1400, line, fontsize=fs)
    # plt.show()
    # print(data)
    best_blade = np.argmax(data['Mintots'] == min(data['Mintots']))
    print(best_blade)
    # if best_blade<20:
    title = 'Best blade is linear. Theta root: '+str(np.round(np.linspace(5, 35, 10)[best_blade], 2))+'°'
    # else:
    #     title = 'Best blade is linear. Theta_root: '+str(np.round(np.linspace(1, 45, 20)[best_blade-20], decimals=1))+' deg'
    # plt.title(title)

    plt.show()
    # idx_opt = best_blade


    # pre = os.path.dirname(os.path.realpath(__file__))
    # fname = 'mission_profiles2.xlsx'
    # path = os.path.join(pre, fname)
    # mission_profiles = pd.read_excel(path, sheet_name='power').values
    # plt.plot(mission_profiles[:,best_blade])
    # plt.show()

def plotPower():
    pre = os.path.dirname(os.path.realpath(__file__))
    # fname = 'mission_profiles2.xlsx'
    fname = 'mission_profiles_fin.xlsx'
    path = os.path.join(pre, fname)

    t_ho    = 60
    t_cl_v  = 30
    t_cl_h  = 470
    t_cr    = 2500
    t_de    = 500

    time = np.arange(0, t_ho+t_cl_v+t_cl_h+t_cr+t_de+t_ho-9, 0.1)


    mission_profiles = pd.read_excel(path, sheet_name='power').values

    plt.plot(mission_profiles[:,0], label= r'Low $\sigma$, low $\theta$', color = 'b')
    plt.plot(mission_profiles[:,1], label= r'Mid $\sigma$, mid $\theta$', color = 'r')
    plt.plot(mission_profiles[:,2], label= r'High $\sigma$, high $\theta$', color = 'k')
    plt.xlabel('Time [ds]', fontsize = 24)
    plt.ylabel('Power [kW]', fontsize = 24)
    plt.grid()

    plt.legend(fontsize=24 ,loc=1)
    # plt.plot(time, mission_profiles[:,0])
    plt.show()

def openPickle(name):
    pre = os.path.dirname(os.path.realpath(__file__))
    fname = name
    path = os.path.join(pre, fname)
    with open(path, 'rb') as handle:
        data = pickle.load(handle)
    # for i in range(len(data)):
    return data

def openExcel(name):
    pre = os.path.dirname(os.path.realpath(__file__))
    fname = name
    path = os.path.join(pre, fname)
    data = pd.read_excel(path, sheet_name='mission_data').values
    # for i in range(len(data)):
    return data

def trysim(h, R, T):
    MTOW = T
    N, B, R, cutout, solidity, theta_tip, taper, airfoil, be = initVariables(R)
    v_inf = 0
    rotor = Rotor(N, B, R, cutout, taper, solidity, theta_tip, airfoil, v_inf, be, h)
    # angle = 0.0018937530751135466
    angle = pi/2
    omega = rotor.maxOmega
    rotor.calcTwist('ideal', 4*pi/180, 1)
    # thrust, power, ct, cp, inflow, idx = rotor.simulation2(500/12.5, 8*np.pi/180)
    wing_chord = 2
    goal = (MTOW+0.3*MTOW*wing_chord*rotor.N/(pi*rotor.R)*np.sin(angle))/rotor.N/(rotor.rho*pi*rotor.R**2*(omega*rotor.R)**2)

    power, thrust, omega, ct, cp, inflow, collective, dct, dcp = flight(rotor, goal, 10)
    # thrust, power, ct, cp, inflow, F = rotor.simulation(500/12.5, 8*np.pi/180)
    print('thrust, power, T/P', thrust, power, thrust/power)
    # plt.plot(rotor.r, F)
    # plt.figure(1)
    # plt.plot(rotor.r, inflow)
    # plt.plot()
    # plt.show()

def finalplots():
    data = openExcel('mission_data.xlsx')
    t = 4
    fs = 20
    # plt.figure(1)
    # plt.plot(data[:, -1], data[:, 0], color='g', linewidth=4, label= 'Thrust')
    # plt.plot([1,2,3], [1,2,3], color= 'r', label='Altitude', linewidth=4)
    # plt.legend()
    # plt.show()
    powers, times, hs, v_infs, angles, mission_profile, mission_time, dcts, omegas, collectives, inflows, dcps, Machs, cts, cps = finalsim(data, rotor, 2)
    Thrusts = np.array(data[:,0])
    vs = (np.sqrt(Thrusts/(2*rotor.area*ISA(data[:,2])))+np.array(data[:,1]))
    # print(vs)
    # plt.plot(times, powers)
    plt.title('Tip Mach number for different max speeds', fontsize=fs)
    plt.xlabel('Speed [km/h]', fontsize=fs)
    plt.ylabel('Mach number [-]', fontsize=fs)
    plt.xticks(fontsize=fs)
    plt.yticks(fontsize=fs)
    # plt.plot(times, Machs, linewidth=t, c='g')
    plt.plot(np.array(data[:,1])*3.6, Machs, linewidth=t, c='g')
    plt.show()


    plt.title(r'Power through the flight', fontsize=fs)
    plt.xlabel('Time [s]', fontsize=fs)
    plt.ylabel('Power [kW]', fontsize=fs)
    plt.xticks(fontsize=fs)
    plt.yticks(fontsize=fs)
    plt.plot(times, powers, linewidth=t, c='g')#, label='total')
    # plt.legend()
    plt.show()


    P_prof = []
    P_ind = []
    for idx in range(len(collectives)):
    # idx = 0
        theta = rotor.twist + collectives[idx]
        aoa_r = theta-inflows[idx]/rotor.r
        closest_aoa_r = np.round(aoa_r*180/pi *4)/4
        index = ((closest_aoa_r - rotor.airfoil.info[0, 0])*4-1).astype(int)
        Mach = np.sqrt((omegas[idx]*rotor.r*rotor.R)**2+rotor.v_inf**2)/rotor.v_sound
        cd = np.array(rotor.airfoil.info[index, 2])/np.sqrt(1-Mach**2)

        cp_prof = 0.5*rotor.sigma*cd*rotor.r**3*rotor.dr
        cp_ind = inflows[0]*dcts[0]*rotor.dr

        cp_prof = np.sum(cp_prof)
        cp_ind = np.sum(cp_ind)

        p_prof = cp_prof*rotor.rho*pi*rotor.R**2*(omegas[idx]*rotor.R)**3
        p_ind = cp_ind*rotor.rho*pi*rotor.R**2*(omegas[idx]*rotor.R)**3

        P_prof.append(p_prof/500)
        P_ind.append(p_ind/500)

    plt.title(r'Power vs Speed', fontsize=fs)
    plt.xlabel('Time [s]', fontsize=fs)
    plt.ylabel('Power [kW]', fontsize=fs)
    plt.xticks(fontsize=fs)
    plt.yticks(fontsize=fs)
    plt.plot(data[:,1], P_prof, linewidth=t, c='g', label= 'profile')
    plt.plot(data[:,1], P_ind, linewidth=t, c='b', label='induced')
    plt.plot(data[:,1], powers, linewidth=t, c='g', label='total')
    plt.legend()
    plt.show()






    print(cps)
    print(cts)
    plt.figure(2)
    plt.plot(data[:,1], np.array(cts)/np.array(cps))
    plt.figure(1)
    for i in range(len(cps)):
        l = str(data[i,1]*3.6)+' km/h'
        plt.scatter(cps[i],cts[i],label=l, linewidth=5)
        # plt.text(cps[i],cts[i],str(data[i,1]), fontsize=fs)
    # plt.scatter(cps,cts)
    plt.title(r'Ct vs Cp for different speeds', fontsize=fs)
    plt.xlabel('Cp [-]', fontsize=fs)
    plt.ylabel('Ct [-]', fontsize=fs)
    plt.xticks(fontsize=fs)
    plt.yticks(fontsize=fs)
    plt.legend(fontsize=fs)
    plt.show()

    FMs = Thrusts*vs/(powers*1000)
    eff = np.array(data[:,0])*np.array(data[:,1])/powers/1000
    plt.plot(times, FMs, label='figure of merit')
    plt.plot(times, eff, label='efficiencies')
    plt.legend()
    plt.show()

    # plt.plot(rotor.r, rotor.twist)
    # for i in collectives:
    # for i in [0, 2, 6]:
    for i in range(len(collectives)):
        if i == 0:
            l = 'Hover'
        elif i == 2:
            l = 'Climb'
        elif i == 6:
            l = 'Cruise'

        plt.plot(rotor.r, collectives[i]+rotor.twist, label=l, linewidth=t)
        # plt.text(rotor.r[0], collectives[i]+rotor.twist[0], str(np.argmax(collectives==collectives[i])))
        # t-=0.4
    plt.title('Geometric angle of attack along the blade', fontsize=fs)
    plt.legend(fontsize=fs)
    plt.xticks(fontsize=fs)
    plt.yticks(fontsize=fs)
    plt.xlabel('Radial position [-]', fontsize=fs)
    plt.ylabel('Angle [rad]', fontsize=fs)
    plt.show()

    for i in range(len(collectives)):
    # for i in [0, 2, 6]:
        if i == 0:
            l = 'Hover'
        elif i == 2:
            l = 'Climb'
        elif i == 6:
            l = 'Cruise'
        theta = rotor.twist + collectives[i]
        aoa_r = theta-inflows[i]/rotor.r
        plt.plot(rotor.r, aoa_r, label=l, linewidth=t)
        # plt.text(rotor.r[0], collectives[i]+rotor.twist[0], str(np.argmax(collectives==collectives[i])))
    plt.title('Angle of attack along the blade', fontsize=fs)
    plt.xlabel('Radial position [-]', fontsize=fs)
    plt.ylabel('Angle [rad]', fontsize=fs)
    plt.legend(fontsize=fs)
    plt.xticks(fontsize=fs)
    plt.yticks(fontsize=fs)
    plt.legend(fontsize=fs)
    plt.show()

    plt.plot(times, powers)
    fig, ax = plt.subplots(3,1)
    # ax[0].plot(rotor.r, dcts[0])
    # ax[1].plot(rotor.r, inflows[0]*rotor.dr*dcts[0], label= 'induced')
    # ax[1].plot(rotor.r, dcp_x, label='profile')


    theta = rotor.twist + collectives[7]
    aoa_r = theta-inflows[7]/rotor.r
    # ax[0].plot(rotor.r, aoa_r)
    ax[0].plot(times, omegas)

    for i, dcp in enumerate(dcps):
        theta = rotor.twist + collectives[i]
        aoa_r = theta-inflows[i]/rotor.r
        ct = np.trapz(dcts[i], rotor.r)+np.array(dcts[i])[-1]*rotor.dr
        cp = np.sum(dcp)
        # print('Iter')
        # print(ct)
        # print(cp)
        # print(ct/cp)
        # ax[0].plot(rotor.r, aoa_r, label=str(i))
        ax[1].scatter(cp, ct, label=str(data[i,1]))
    ax[1].legend()
    # ax[0].legend()
    ax[2].plot(rotor.r, np.array(dcts[0])/np.array(dcps[0]))
    plt.show()
    return powers, data[:,1], inflows, collectives

if __name__ == "__main__":
    warnings.filterwarnings("ignore")
    # powers_h, twists = hoverFlight(0, 3.25, 4000*9.81, 2)
    # powers_f = forwardFlight(2000, 3.25, 400*9.81, 100)
    # print()
    N, B, R, cutout, solidity, theta_tip, taper, airfoil, be = initVariables(4.32)
    rotor = Rotor(N, B, R, cutout, taper, solidity, theta_tip, airfoil, 0, be, 0)
    rotor.calcTwist('linear', 18.3*pi/180, 0*pi/180)
    rotor.cond = True
    powers_notwist, times, inflow_nt, col_nt = finalplots()
    # rotor.cond = False
    # powers_twist, times, inflow_t, col_t = finalplots()


    # rotor.cond=False
    # powers_twist, times, inflow_nt, col_nt = finalplots()
    #
    #
    # fs = 20
    # plt.plot(np.array(times)*3.6, powers_notwist, label='With drag divergence', c='g', linewidth=4)
    # plt.plot(np.array(times)*3.6, powers_twist, label='Without drag divergence', c='r', linewidth=4)
    # plt.title('Powers at different speeds', fontsize=fs)
    # plt.xlabel('Speed [km/h]', fontsize=fs)
    # plt.ylabel('Power [kW]', fontsize=fs)
    # plt.legend(fontsize=fs)
    # plt.xticks(fontsize=fs)
    # plt.yticks(fontsize=fs)
    # plt.show()

    check = False
    if check:
        # rotor.calcTwist('linear', 30*pi/180, 0*pi/180)
        # print('speed of sound 2000')
        # ISA(2000)
        rotor.sigma = 0.07
        powers_twist, times, inflow_t, col_t = finalplots()
        fs = 20
        plt.plot(times, powers_notwist, linewidth=4, c='g', label=r'$\sigma$ = 0.02')
        plt.plot(times, powers_twist, linewidth=4, c='r', label=r'$\sigma$ = 0.07')
        plt.title('Powers throught the mission', fontsize=fs)
        plt.xlabel('Time [s]', fontsize=fs)
        plt.ylabel('Power [kW]', fontsize=fs)
        plt.legend(fontsize=fs)
        plt.xticks(fontsize=fs)
        plt.yticks(fontsize=fs)
        # plt.legend(fontsize=fs)
        plt.show()

        theta_nt = rotor.twist + col_nt[7]
        aoa_r_nt = theta_nt-inflow_nt[7]/rotor.r
        theta_t = rotor.twist + col_t[7]
        aoa_r_t = theta_nt-inflow_t[7]/rotor.r

        fig, ax = plt.subplots(2,1)
        plt.title('Angle of attack distribution, for no twist blade (top) and $0.52 rad$ twist (bottom), during cruise')
        ax[0].plot(rotor.r, np.array(aoa_r_nt)*180/pi, linewidth=4, c='g')
        ax[1].plot(rotor.r, np.array(aoa_r_t)*180/pi, linewidth=4, c='g')
        ax[0].set_ylabel('Angle [rad]')
        ax[0].set_xlabel('Radial position [-]')
        ax[1].set_ylabel('Angle [rad]')
        ax[1].set_xlabel('Radial position [-]')
        plt.show()



    # theta = rotor.twist + collectives[0]
    # aoa_r = theta-inflows[0]/rotor.r
    # closest_aoa_r = np.round(aoa_r*180/pi *4)/4
    # index = ((closest_aoa_r - rotor.airfoil.info[0, 0])*4-1).astype(int)
    # cd = rotor.airfoil.info[index, 2]
    # dcp_x = 0.5*rotor.sigma*cd*rotor.r**3*rotor.dr
    # ax[0].plot(rotor.r, inflows[0])

    # for i in range(len(lst)):
    #     lst_i = lst[i]
    #     print(lst_i)
    #     simplesim(4.45, lst_i)

    # trysim(0, 4.44,  41716.29888128569)
    # finalsim([[41716.29888128569, 7.98770778188738, 1.6, 1.493065081699371, 60.2]], rotor, 2)
    # print('2000 meters')
    # simplesim(2000, 4.36, 3800*9.81)
    # power_i, omega_i = nonfinalsim(rotor, 0, 40000, 2, 0, 0)
    # print(power_i, omega_i)
    # fig, ax1 = plt.subplots(2,1)
    # ax1[0].plot(twists, powers_h, label='hover power')
    # ax1[0].plot(twists, powers_f, label='forward power')
    # ax1[0].legend()
    # plt.figure(5)
    # plt.plot(twists, 2*np.array(powers_h), label='hover power')
    # plt.plot(twists, 2*np.array(powers_f), label='forward power')
    # plt.legend()
    # combinedFlight(0, 2000, 0, 0, 3.25, 4000, 20.29, 12.34, 2)
    # lst= [[41716.29888128569, 7.98770778188738, 1.6, 1.493065081699371, 60.2]]
    # rotor.calcTwist('linear', 18.3*pi/180, 0*pi/180)
    # calcPowers(rotor)
    # plotMass()
    # plotPower()
    # ctcpCurve(rotor)
    # plotCtCpExp()
    # N, B, R, cutout, solidity, theta_tip, taper, airfoil, be = initVariables(4.45)
    # rotor = Rotor(N, B, R, cutout, taper, solidity, theta_tip, airfoil, 0, be, 0)
    # for rotsol in np.linspace(0.02, 0.04, 10):
    #     rotor.sigma = rotsol
    #     print(finalsim(lst, rotor, 2), ' at ', rotsol)
    # makeMP()
    # simplesim2()
    # trysim(0, 4.4, 4000*9.81)





    # plt.show()
