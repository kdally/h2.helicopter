import numpy as np
from math import pi, pow
import matplotlib.pyplot as plt
import time
import pandas as pd
from IterationTools.aerodynamics.BET import Rotor, Airfoil

'''
Code for optimising the rotor performance for both, hover and forward flight.



'''

N = 2
# B = 2
B = 3
R = 3.25
# R = 3.81
# cutout = 0.15
# cutout = 0.00001
cutout = 0.2
solidity = 0.08
solidity = 0.06928
# solidity = 0.027
theta_tip = 10*pi/180
theta_0 = 20*pi/180
# t_ratio = 1
taper = 1
airfoil = '0012'
MTOM = 4000
# MTOM = 400
MTOW = MTOM*9.81
v_inf = 0
be = 100



powers = []
omegas = []
thrusts = []
t_ratios = []
dcts = []
cts = []
cps = []
dcts = []
dcps = []
inflows = []
cols = []
twists = []
h = 0
collective = 0*pi/180
opt_twist = False
one_rpm = False
no_twist = True
real_opt_twist = False
ctcp = False
# omega = 500/12.5
# collective = 0
#
# rotor = Rotor(N, B, R, cutout, solidity, theta_tip, airfoil, h)
# thrust, power, ct, cp, inflow = rotor.simulation(omega, collective)
# print('thrst', thrust)
# print('power', power)
# print('ct', ct)
# print('cp', cp)
# print('inflow', inflow)

# iterating over RPM and twist ratio
# t1 = time.time()
if opt_twist:
    for t_ratio in np.linspace(0., 1, 10):
        # for theta_root in np.arange(theta_tip, 50*pi/180, 1*pi/180):
        print(t_ratio)
        power_i = []
        thrust_i = []
        omega_i = []
        ct_i = []
        cp_i = []
        col_i = []
        inflow_i = []
        index = []

        t1 = time.time()

        rotor = Rotor(N, B, R, cutout, taper, solidity, theta_tip, airfoil, v_inf, be, h)
        # rotor.calcTwist('ideal', theta_tip, t_ratio)
        rotor.calcTwist('linear', theta_0*t_ratio, 0)
        plt.plot(rotor.r, rotor.twist)
        plt.show()
        # rotor.calcTwist('poly', a, b, c)
        print('Maximum rpms', rotor.maxOmega)

        for omega in np.linspace(30, rotor.maxOmega, 15):
            thrustmatch = False
            collective = 0*pi/180
            print(omega)

            # goal = MTOW/N/(0.5*rotor.rho*(omega*rotor.R)**2*rotor.area)
            goal = MTOW/N/(rotor.rho*pi*rotor.R**2*(omega*rotor.R)**2)
            while not thrustmatch:
                print('col', collective)
                thrust, power, ct, cp, inflow = rotor.simulation(omega, collective)
                # print(thrust)
                # print('ct', ct)
                # print(collective)
                if ct == 1000:
                    print('sim not valid')
                    power_i.append(10**10)
                    thrust_i.append(10**10)
                    omega_i.append(10**10)
                    ct_i.append(10**10)
                    cp_i.append(10**10)
                    inflow_i.append(10**10)
                    col_i.append(10**10)
                    thrustmatch = True

                else:
                    diff = ct - goal
                    # diff = thrust-MTOW/N
                    # print('try')
                    # print(diff)
                    # print(diff)
                    if diff > 0.0000001:
                        collective -= 10 * diff
                    elif diff < -0.0000001:
                        collective -= 10 * diff

                    else:
                        thrustmatch = True
                        power_i.append(power)
                        thrust_i.append(thrust)
                        omega_i.append(omega)
                        ct_i.append(ct)
                        cp_i.append(cp)
                        inflow_i.append(inflow)
                        col_i.append(collective)

        # t2 = time.time()
        # print(t2-t1)
        if min(cp_i) != 10**10:
            opt_list = np.array(cp_i)/np.array(ct_i)
            opt_index = list(opt_list).index(min(opt_list))

            # print(rotor.twist+collective-inflow/rotor.r)
            powers.append(power_i[opt_index])
            thrusts.append(thrust_i[opt_index])
            omegas.append(omega_i[opt_index])
            t_ratios.append(t_ratio)
            cts.append(ct_i[opt_index])
            cps.append(cp_i[opt_index])
            twists.append(rotor.twist)
            inflows.append(inflow_i[opt_index])
            cols.append(col_i[opt_index])

            aoa_r = rotor.twist+col_i[opt_index]-inflow/rotor.r

            lift = B*np.sum((rotor.airfoil.cl_alpha*aoa_r*0.5*rotor.rho*((inflow*omega*rotor.R)**2+omega**2*(rotor.r*rotor.R)**2)*rotor.dr*rotor.c))
            thrust_from_ct = ct*0.5*rotor.rho*(omega*rotor.R)**2*pi*rotor.R**2
            # print('lift', lift)
            # print('thrust', thrust_from_ct)
        else:
            print('Non valid twist ratio')

    fig, axs = plt.subplots(4,1)

    axs[0].plot(t_ratios, np.array(powers)/1000)
    axs[0].set_ylabel('power [kW]')
    # plt.show()

    axs[1].plot(t_ratios, thrusts)
    axs[1].set_ylabel('thrust [N]')
    # plt.show()

    axs[2].plot(t_ratios, omegas)
    axs[2].set_ylabel('omega [rad/s]')


    k = 1.15
    cd0 = 0.011
    cps_theory = k*np.array(cts)**(3/2)/(np.sqrt(2))+(solidity*cd0/8)
    print('theory', cps_theory)
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

    aoa_r = rotor.twist+collective-inflow/rotor.r
    lift = (rotor.airfoil.cl_alpha*aoa_r*0.5*rotor.rho*((inflow*omega*rotor.R)**2+omega**2*(rotor.r*rotor.R)**2)*rotor.dr*rotor.c)

    plt.figure(1)
    fig2, ax = plt.subplots(4,1)

    ax[0].set_title('Optimum twist ratio: '+ str(t_ratios[opt_index]))
    ax[0].plot(rotor.r, inflow)
    ax[0].set_ylabel('inflow')

    ax[1].plot(rotor.r, aoa_r)
    ax[1].set_ylabel('angle of attack')

    # ax[2].plot(rotor.r, dct)
    # ax[2].set_ylabel('dct')

    ax[3].plot(rotor.r, twist*180/pi)
    ax[3].set_ylabel('twist')
    # axs[3].plot(rotor.r, np.sqrt((rotor.r*omega)**2+(inflow*omega*rotor.R)**2))
    # ax[3].plot(cp, ct)
    plt.show()


    # plt.plot(cp, ct)
    # plt.show()

# See results for one RPM
elif one_rpm:
    for t_ratio in np.arange(0.1,1.1, 0.1):
        print(t_ratio)
        power_i = []
        thrust_i = []
        omega_i = []
        t1 = time.time()
        # for omega in np.arange(30, 100, 1):
        thrustmatch = False
        collective = 0*pi/180

        rotor = Rotor(N, B, R, cutout,taper, solidity, theta_tip, airfoil, v_inf, be, h)
        # rotor.calcTwist('ideal', theta_tip, t_ratio)
        # rotor.calcTwist('linear', theta_0, theta_tip)
        # rotor.calcTwist('poly', a, b, c)

        omega = min(rotor.maxOmega, 85)
        print(omega)

        # goal = MTOW/N/(0.5*rotor.rho*(omega*rotor.R)**2*rotor.area)
        goal = MTOW/N/(rotor.rho*pi*rotor.R**2*(omega*rotor.R)**2)
        while not thrustmatch:

            thrust, power, ct, cp, inflow = rotor.simulation(omega, collective)
            # print(thrust)
            if ct == 1000:
                print('simulation not valid')
                omega += 1
                collective -= 0.1
                print(omega)
                goal = MTOW/N/(rotor.rho*pi*rotor.R**2*(omega*rotor.R)**2)

            else:
                diff = ct - goal
                # diff = thrust-MTOW/N
                # print('try')
                # print(diff)

                if diff > 0.00001:
                    collective -= 10 * diff
                elif diff < -0.00001:
                    collective -= 10 * diff

                else:
                    thrustmatch = True
            # power_i.append(power)
            # thrust_i.append(thrust)
            # omega_i.append(omega)
        t2 = time.time()
        print(t2-t1)
        # opt_index = power_i.index(min(power_i))
        # power = power_i[opt_index]
        # thrust = thrust_i[opt_index]
        # omega = omega_i[opt_index]

        # print(rotor.twist+collective-inflow/rotor.r)
        powers.append(power)
        cols.append(collective)
        inflows.append(inflow)
        omegas.append(omega)
        t_ratios.append(t_ratio)
        thrusts.append(thrust)
        cps.append(cp)
        cts.append(ct)
        # dcts.append(dct)




        aoa_r = rotor.twist+collective-inflow/rotor.r

        lift = (rotor.airfoil.cl_alpha*aoa_r*0.5*rotor.rho*((inflow*omega*rotor.R)**2+omega**2*(rotor.r*rotor.R)**2)*rotor.dr*rotor.c)

        print('lift', np.sum(lift*rotor.dr))
        print('thrust', (ct*rotor.rho*pi*rotor.R**2*(omega*rotor.R)**2))

    fig, axs = plt.subplots(4,1)
    axs[0].plot(t_ratios, np.array(powers)/1000)
    axs[0].set_ylabel('power [kW]')
    # plt.show()

    axs[1].plot(t_ratios, thrusts)
    axs[1].set_ylabel('thrust [N]')
    # plt.show()

    axs[2].plot(t_ratios, omegas)
    axs[2].set_ylabel('omega [rad/s]')


    k = 1.15
    cd0 = 0.011
    cps_theory = k*np.array(cts)**(3/2)/(np.sqrt(2))+(solidity*cd0/8)
    print('theory', cps_theory)
    ratios = np.array(cts)/np.array(cps)

    axs[3].scatter(cps, cts, label='single')
    axs[3].plot(cps_theory, cts, label='theory')
    axs[3].scatter(cps[np.where(ratios == np.max(ratios))[0][0]], cts[np.where(ratios == np.max(ratios))[0][0]], label='Opt point')
    axs[3].set_ylabel('ct to cp')
    axs[3].set_xlim(min(min(cps), min(cps_theory))-0.00005, max(max(cps), max(cps_theory))+0.00005)
    axs[3].set_ylim(min(cts)-0.001, max(cts)+0.001)
    axs[3].legend()




    opt_index = powers.index(min(powers))
    collective = cols[opt_index]
    inflow = inflows[opt_index]
    omega = omegas[opt_index]
    # dct = dcts[opt_index]

    aoa_r = rotor.twist+collective-inflow/rotor.r
    lift = (rotor.airfoil.cl_alpha*aoa_r*0.5*rotor.rho*((inflow*omega*rotor.R)**2+omega**2*(rotor.r*rotor.R)**2)*rotor.dr*rotor.c)



    fig, axs = plt.subplots(3,1)

    axs[0].set_title(t_ratios[opt_index])
    axs[0].plot(rotor.r, inflow)
    axs[0].set_ylabel('inflow')

    axs[1].plot(rotor.r, aoa_r)
    axs[1].set_ylabel('angle of attack')

    # axs[2].plot(rotor.r, dct)
    # axs[2].set_ylabel('dct')

    # axs[3].plot(rotor.r, np.sqrt((rotor.r*omega)**2+(inflow*omega*rotor.R)**2))
    plt.show()

#Iterating over a range of RPM's
elif no_twist:
    rotor = Rotor(N, B, R, cutout, taper, solidity, theta_tip, airfoil, v_inf, be, h)
    # rotor.calcTwist('ideal', theta_tip)
    power_i = []
    thrust_i = []
    omega_i = []
    dct_i = []
    ct_i = []
    cp_i = []
    col_i = []
    inflow_i = []

    # for omega in np.linspace(30,rotor.maxOmega, 20):
    for omega in [85]:
        print(omega)
        thrustmatch = False
        collective = 0*pi/180

        goal = (1+0.3*2/(pi*3.25))*MTOW/N/(rotor.rho*pi*rotor.R**2*(omega*rotor.R)**2)
        while not thrustmatch:
            thrust, power, ct, cp, inflow = rotor.simulation(omega, collective)
            # print(thrust)
            if ct == 1000:
                print('sim not valid')
                power = 10**10
                thrust = 10**10
                omega = 10**10
                thrustmatch = True

            else:
                diff = ct-goal


                if diff > 0.0000001:
                    collective -= 10 * diff
                elif diff < -0.0000001:
                    collective -= 10 * diff

                else:
                    thrustmatch = True
                    powers.append(power)
                    thrusts.append(thrust)
                    omegas.append(omega)
                    # dcts.append(dct)
                    cts.append(ct)
                    cps.append(cp)
                    inflows.append(inflow)
                    cols.append(collective)
                    # dcts.append(dct)


    fig, axs = plt.subplots(3,1)

    print('power =', powers[0]*2/1000, 'kW')
    axs[0].plot(omegas, np.array(powers)*2)
    axs[0].set_ylabel('power')

    # axs[0].plot(omegas, cps)
    # axs[0].set_ylabel('power coefficient')
    # plt.show()

    axs[1].plot(omegas, np.array(thrusts)*2)
    axs[1].set_ylabel('Thrust')
    # plt.show()

    k = 1.15
    cd0 = 0.011
    cps_theory = k*np.array(cts)**(3/2)/(np.sqrt(2))+(solidity*cd0/8)
    ratios = np.array(cts)/np.array(cps)

    axs[2].scatter(cps, cts, label='single')
    axs[2].plot(cps_theory, cts, label='theory')
    axs[2].scatter(cps[np.where(ratios == np.max(ratios))[0][0]], cts[np.where(ratios == np.max(ratios))[0][0]], label='Opt point')
    axs[2].set_ylabel('ct to cp')
    axs[2].set_xlim(min(min(cps), min(cps_theory))-0.00005, max(max(cps), max(cps_theory))+0.00005)
    axs[2].set_ylim(min(cts)-0.001, max(cts)+0.001)
    axs[2].legend()


    opt_index = powers.index(min(powers))
    collective = cols[opt_index]
    inflow = inflows[opt_index]
    omega = omegas[opt_index]
    power = powers[opt_index]
    # dct = dcts[opt_index]

    aoa_r = rotor.twist+collective-inflow/rotor.r
    lift = (rotor.airfoil.cl_alpha*aoa_r*0.5*rotor.rho*((inflow*omega*rotor.R)**2+omega**2*(rotor.r*rotor.R)**2)*rotor.dr*rotor.c)

    inflow[-1] = 0
    plt.figure(1)
    fig2, ax = plt.subplots(2,1)

    ax[0].set_title('Omega = '+str(np.around(omegas[opt_index], 2))+' [rad/s], Col = '+str(np.around(collective, 2))+' [rad], Power = '+str(np.around(power/1000, 2))+' [kW]')
    ax[0].plot(rotor.r, inflow)
    ax[0].set_ylabel('inflow')

    ax[1].plot(rotor.r, aoa_r)
    ax[1].set_ylabel('angle of attack')

    # ax[2].plot(rotor.r, dct)
    # ax[2].set_ylabel('dct')
    # ax[2].bar(1, power)

    # axs[3].plot(rotor.r, np.sqrt((rotor.r*omega)**2+(inflow*omega*rotor.R)**2))
    # ax[3].plot(cps, cts)
    plt.show()

#Constant angle of attach through blade over a range of RPM's
elif real_opt_twist:
    optimised = False
    rotor = Rotor(N, B, R, cutout, taper, solidity, theta_tip, airfoil, v_inf, be, h)

    while not optimised:

        # rotor.calcTwist('ideal', theta_tip)
        power_i = []
        thrust_i = []
        omega_i = []
        # dct_i = []
        ct_i = []
        cp_i = []
        col_i = []
        inflow_i = []

        for omega in np.arange(30,100, 5):
            print(omega)
            thrustmatch = False
            collective = 0*pi/180

            goal = MTOW/N/(rotor.rho*pi*rotor.R**2*(omega*rotor.R)**2)
            while not thrustmatch:
                thrust, power, ct, cp, inflow = rotor.simulation(omega, collective)
                # print(thrust)

                if ct == 1000:
                    print('sim not valid')
                    power = 10**10
                    thrust = 10**10
                    omega = 10**10
                    thrustmatch = True

                else:
                    diff = ct-goal


                    if diff > 0.0000001:
                        collective -= 10 * diff
                    elif diff < -0.0000001:
                        collective -= 10 * diff

                    else:
                        thrustmatch = True
                        powers.append(power)
                        thrusts.append(thrust)
                        omegas.append(omega)
                        # dcts.append(dct)
                        cts.append(ct)
                        cps.append(cp)
                        inflows.append(inflow)
                        cols.append(collective)

        opt_index = powers.index(min(powers))
        collective = cols[opt_index]
        inflow = inflows[opt_index]

        aoa_r = rotor.twist+collective-inflow/rotor.r
        diff = aoa_r - np.min(aoa_r)
        print(diff)
        print(rotor.twist)


        fig, axs = plt.subplots(3,1)

        axs[0].scatter(omegas, powers)
        axs[0].set_ylabel('power')
        # plt.show()

        axs[1].scatter(omegas, thrusts)
        axs[1].set_ylabel('thrust')
        # plt.show()

        k = 1.09
        cd0 = 0.011
        cps_theory = k*np.array(cts)**(3/2)/(np.sqrt(2))+(solidity*cd0/8)

        axs[2].scatter(cps, cts, label='single')
        axs[2].scatter(cps_theory, cts, label='theory')
        axs[2].set_ylabel('ct to cp')
        axs[2].set_xlim(min(cps)-0.00005, max(cps)+0.00005)
        axs[2].set_ylim(min(cts)-0.001, max(cts)+0.001)
        axs[2].legend()


        opt_index = powers.index(min(powers))
        collective = cols[opt_index]
        inflow = inflows[opt_index]
        omega = omegas[opt_index]
        # dct = dcts[opt_index]

        aoa_r = rotor.twist+collective-inflow/rotor.r
        lift = (rotor.airfoil.cl_alpha*aoa_r*0.5*rotor.rho*((inflow*omega*rotor.R)**2+omega**2*(rotor.r*rotor.R)**2)*rotor.dr*rotor.c)

        inflow[-1] = 0
        fig2, ax = plt.subplots(3,1)

        ax[0].set_title('Omega = '+str(omegas[opt_index]))
        ax[0].plot(rotor.r, inflow)
        ax[0].set_ylabel('inflow')

        ax[1].plot(rotor.r, aoa_r, label='aoa')
        ax[1].plot(rotor.r, rotor.twist, label='twist')
        ax[1].set_ylabel('angle of attack')
        ax[1].legend()

        # ax[2].plot(rotor.r, dct)
        # ax[2].set_ylabel('dct')
        plt.show()

        if np.average(diff)>0.0001:
            rotor.calcTwist('coords', diff)
        else:
            optimised = True

    # plt.figure(0)
    fig, axs = plt.subplots(3,1)

    axs[0].plot(omegas, powers)
    axs[0].set_ylabel('power')
    # plt.show()

    axs[1].plot(omegas, thrusts)
    axs[1].set_ylabel('thrust')
    # plt.show()

    k = 1.09
    cd0 = 0.011
    cps_theory = k*np.array(cts)**(3/2)/(np.sqrt(2))+(solidity*cd0/8)

    axs[2].scatter(cps, cts, label='single')
    axs[2].plot(cps_theory, cts, label='theory')
    axs[2].plot()
    axs[2].set_ylabel('ct to cp')
    axs[2].set_xlim(min(min(cps), min(cps_theory))-0.00005, max(max(cps), max(cps_theory))+0.00005)
    axs[2].set_ylim(min(cts)-0.001, max(cts)+0.001)
    axs[2].legend()


    opt_index = powers.index(min(powers))
    collective = cols[opt_index]
    inflow = inflows[opt_index]
    omega = omegas[opt_index]
    # dct = dcts[opt_index]

    aoa_r = rotor.twist+collective-inflow/rotor.r
    lift = (rotor.airfoil.cl_alpha*aoa_r*0.5*rotor.rho*((inflow*omega*rotor.R)**2+omega**2*(rotor.r*rotor.R)**2)*rotor.dr*rotor.c)

    inflow[-1] = 0
    plt.figure(1)
    fig2, ax = plt.subplots(3,1)

    ax[0].set_title('Omega = '+str(omegas[opt_index]))
    ax[0].plot(rotor.r, inflow)
    ax[0].set_ylabel('inflow')

    ax[1].plot(rotor.r, aoa_r)
    ax[1].set_ylabel('angle of attack')

    # ax[2].plot(rotor.r, dct)
    # ax[2].set_ylabel('dct')

    # axs[3].plot(rotor.r, np.sqrt((rotor.r*omega)**2+(inflow*omega*rotor.R)**2))
    # ax[3].plot(cps, cts)
    plt.show()

#Iterate over range of collectives for flat blades, fixed RPM's
elif ctcp:
    rotor = Rotor(N, B, R, cutout, taper, solidity, theta_tip, airfoil, v_inf, be, h)
    # rotor.calcTwist('ideal', theta_tip)
    power_i = []
    thrust_i = []
    omega_i = []
    ct_i = []
    cp_i = []
    col_i = []
    inflow_i = []

    # for omega in np.arange(30,100, 5):
    omega = 500/12.5
    omega = 40 # rad/s
    collective = 0*pi/180

    goal = MTOW/N/(rotor.rho*pi*rotor.R**2*(omega*rotor.R)**2)
    for t in np.arange(0, 15, 1):
        t = t*pi/180
        rotor.calcTwist('flat', t)
        thrust, power, ct, cp, inflow = rotor.simulation(omega, collective)

        powers.append(power)
        thrusts.append(thrust)
        omegas.append(omega)
        cts.append(ct)
        cps.append(cp)
        inflows.append(inflow)
        cols.append(collective)
        # if 0.0042 < ct <0.044:
        #     cts.append(dct)

    # plt.figure(0)
    fig, axs = plt.subplots(2,1)
    #
    # axs[0].plot(omegas, powers)
    # axs[0].set_ylabel('power')
    # # plt.show()
    #
    # axs[1].plot(omegas, thrusts)
    # axs[1].set_ylabel('thrust')
    # plt.show()

    k = 1.15
    cd0 = 0.011
    cps_theory = k*np.array(cts)**(3/2)/(np.sqrt(2))+(solidity*cd0/8)

    axs[0].scatter(cps, cts, label='single')
    axs[0].plot(cps_theory, cts, label='theory')
    axs[0].set_ylabel('ct to cp')
    # axs[0].set_xlim(min(min(cps), min(cps_theory))-0.00005, max(max(cps), max(cps_theory))+0.00005)
    axs[0].set_xlim(min(min(cps), min(cps_theory))-0.00005, 0.000433)
    axs[0].set_ylim(min(cts)-0.001, max(cts)+0.001)
    axs[0].legend()

    axs[1].plot(rotor.r, inflow)
    axs[1].set_ylim(0, max(inflow)+0.01)
    axs[1].set_ylabel('inflow')

    plt.show()

    # opt_index = powers.index(min(powers))
    # collective = cols[opt_index]
    # inflow = inflows[opt_index]
    # omega = omegas[opt_index]
    #
    # aoa_r = rotor.twist+collective-inflow/rotor.r
    # lift = (rotor.airfoil.cl_alpha*aoa_r*0.5*rotor.rho*((inflow*omega*rotor.R)**2+omega**2*(rotor.r*rotor.R)**2)*rotor.dr*rotor.c)
    #
    # inflow[-1] = 0
    # plt.figure(1)
    # fig2, ax = plt.subplots(3,1)
    #
    # ax[0].set_title('Omega = '+str(omegas[opt_index])+', Col = '+str(collective))
    # ax[0].plot(rotor.r, inflow)
    # ax[0].set_ylabel('inflow')
    #
    # ax[1].plot(rotor.r, aoa_r)
    # ax[1].set_ylabel('angle of attack')
    #
    # ax[2].plot(rotor.r, lift)
    # ax[2].set_ylabel('lift')
    #
    # # axs[3].plot(rotor.r, np.sqrt((rotor.r*omega)**2+(inflow*omega*rotor.R)**2))
    # # ax[3].plot(cps, cts)
    # plt.show()
