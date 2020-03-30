#global imports
import numpy as np
import math
import warnings
import time
import pickle
import os
import matplotlib.pyplot as plt
import pandas as pd
import shelve
# import dill

#tool imports
from Mission import TiltRotor
import fpp.WeightFPPfinal as pp
import control_stability.ControlStabilityMain as con_sta
import load_calculation.loadingDiagrams as ld
import aerodynamics.cabin_sizing as cabin
import load_calculation.landing_gear_sizing as lg
import aerodynamics.Fuselage as body
import aerodynamics.Airfoil as airfoil
import aerodynamics.WingLoading as wl
import aerodynamics.RotorPerformance as rp
from aerodynamics.BET import Rotor
import load_calculation.wing_sizing as ws




"""

Iterative Tool to size different aspects of the tilt-rotor aircraft for Futura, group 19 of the DSE. Includes Class I performance, rotor, powerplant, fuselage aerodynamics, S&C.

Created: 03-06-2019 by Futura

"""

def initialiseClassI():
    MTOW = 4000. #Maximum Take-Off Weight [kg]
    OEW = 2600. #Operational Empty Weight
    diskloading = 60. #Disk Loading [kg/m^2]
    stall = 50 #Stall speed [m/s]
    A = 5.81 #Aspect Ratio [-]
    g = 9.81
    wing_area = 23.24
    return MTOW, OEW, diskloading, stall, A, wing_area, g

def initialiseRotorSizing():
    return

def initialisePowerplant():
    return

def initialiseFuselage():
    fuselage_length = 9.73
    fuselage_mass = 327
    fuselage_width = 1.5
    fuselage_height = 1.4
    return fuselage_length, fuselage_mass, fuselage_width, fuselage_height

def initialiseLoading():
    cl_stall_sealevel = 1.386
    cl_cruise_sealevel = 0.436
    cl_dive_sealevel = 0.436
    cl_stall_cruise = 1.386
    cl_cruise_cruise = 0.436
    cl_dive_cruise = 0.436
    cl_max_flaperons=2
    lateral_drag_area = 2
    front_drag_area = 1.2
    cruise_height = 2000
    Cl_alphaw = 6.607
    lt = 3
    lcg = 0.3
    wing_area = 21.48 #probably change in iteration
    chord = 2
    vc_to = 97.222
    vc_cr = 97.222
    return cl_stall_sealevel, cl_cruise_sealevel, cl_dive_sealevel, cl_stall_cruise, cl_cruise_cruise, cl_dive_cruise,cl_max_flaperons, lateral_drag_area, front_drag_area, cruise_height, Cl_alphaw, lt, lcg, wing_area, chord, vc_to, vc_cr

def initialiseAerodynamics():
    MAC = 2
    Ah=30
    Lambda_ch=20*math.pi/180
    Cl_alphaw=6.0
    Cm_acw=-0.05 # moment coefficient aerodynamic center 
    Cm_acf=-0.15# moment coefficient aerodynamic center
    lh= 6.1
    lf=-0.3
    b=14
    xac_w= 0.25 # position of the ac of the main wing as percentage of the MAC
    Cl_w=0.4
    Cl_f=0.2
    Cl_alphaf=3
    Sw=20
    Sf=32
    e = 0.9
    climb_rate = 15.983
    return MAC,Ah,Lambda_ch,Cl_alphaw,Sw,Cm_acf,Cm_acw,lh,b,Cl_f,Cl_alphaf,xac_w,Cl_w,Sf,lf,e,climb_rate

def initialiseControlStability(fuselage_length):
    fixed = [fuselage_length*0.35, 625.8]
    horizontal = [fuselage_length*0.95,75]
    vertical = [fuselage_length*0.95,43]
    xLEMAC = 0.41*fuselage_length
    payload = 900
    passengermass = 90
    SM = 0.04
    xcg = 0.41*fuselage_length + 1
    lm = 0.28
    ln = 3.4
    return fixed, horizontal, vertical, xLEMAC, payload, passengermass, SM, xcg, lm, ln

def temperature(height):
    return 288.15-0.0065*height

def density(height):
    temp = 288.15-0.0065*height
    pressure = 101325*math.pow(temp/288.15,(-9.80665/(-0.0065*287.05)))
    return pressure/(287.05*temp)

def convertTuple(tup):
    tup = list(tup)
    for i in range(len(tup)):
        tup[i] = float(tup[i])
    return tuple(tup)

#=====================================================================================================================================================================================================
#intialise values
MTOW, OEW, diskloading, stall, A, wing_area, g = initialiseClassI()
fuselage_length, fuselage_mass, fuselage_width, fuselage_height = initialiseFuselage()
fixed, horizontal, vertical, xLEMAC, payload, passengermass, SM, xcg, lm, ln = initialiseControlStability(fuselage_length)
cl_stall_sealevel, cl_cruise_sealevel, cl_dive_sealevel, cl_stall_cruise, cl_cruise_cruise, cl_dive_cruise,cl_max_flaperons, lateral_drag_area, front_drag_area, cruise_height, Cl_alphaw, lt, lcg, wing_area, chord, vc_to, vc_cr = initialiseLoading()
MAC, Ah, Lambda_ch, Cl_alphaw, Sw, Cm_acf, Cm_acw, lh, b, Cl_f, Cl_alphaf, xac_w, Cl_w, Sf, lf, e, climb_rate = initialiseAerodynamics()

MTOW_iterations = [4000]
#=====================================================================================================================================================================================================
percentages = [1, 0.95,0.93,0.91,0.89,0.87]
MTOW_percentages = []
for j in range(len(percentages)):
    for i in range(5):
        fuselage_mass = 327*percentages[j]
        cond = False
        if i == 30: cond = True
        #=====================================================================================================================================================================================================
        ##Wing Loading
        rho = density(cruise_height)
        cl_stall_cruise = wl.CLmaxwing_calc(0.9)
        cl_stall_sealevel=cl_stall_cruise
        Mach = wl.Mcruise_calc(vc_cr, temperature(cruise_height), 1.4, 287.15)
        e_wing = wl.OswaldEfficiencyFActor(A)
        CD0_wing = wl.ZeroLiftDragCoeff()
        ypoint1 = wl.ypoint1(density(cruise_height), stall, cl_stall_cruise, 1.225, 0.75, MAC, A, e_wing, CD0_wing)
        ypoint2 = wl.ypoint2(stall, cl_stall_cruise, 2.5, vc_cr, A, e_wing, density(cruise_height), CD0_wing, 0.75)
        WoverS = wl.dp_calc(density(cruise_height), stall, cl_stall_cruise, 1.225, 0.75, climb_rate, A, 2.5, vc_cr)[0]
        wing_area = wl.OptimumWingArea(MTOW, g, 1.225, 0.75, climb_rate, A, 2.5, vc_cr, density(cruise_height), stall, cl_stall_cruise)

        Cl_alphaw = wl.WingLiftCurveSlope(A, Mach, 0.95)
        CL_list_wing = wl.CL_calc(Cl_alphaw, -19, 20, -1.6)
        CD_list_wing = wl.CD_calc(CD0_wing, A, e_wing, Mach, 0.95, -19, 20, -1.6)
        alpha_trim = wl.TrimAngleCalc(0.436, Cl_alphaw, -1.6)
        # cl_cruise_cruise = wl.CLtrimcalc(Cl_alphaw, -1.6, alpha_trim)
        cl_cruise_sealevel=cl_cruise_cruise
        cl_dive_sealevel=cl_cruise_sealevel
        cl_dive_cruise=cl_cruise_cruise
        stall_angle_wing = wl.StallAngleWing_calc(cl_stall_cruise, Cl_alphaw, -1.6, 2.2)
        delta_CL_max = wl.DeltaCLmax_calc(1.9, cl_stall_cruise)
        cl_max_flaperons=cl_stall_cruise+delta_CL_max
        SwfoverS = wl.SwfoverS_calc(delta_CL_max, 0.9)
        flapchord = wl.cflap_calc(0.25, chord)
        lflap = wl.lflap_calc(SwfoverS, wing_area)
        alpha_stall_flap = wl.alphastall_flap_calc(stall_angle_wing, SwfoverS)

        Cl_alphaw = Cl_alphaw*180/math.pi

        print("Wingloading Complete")
        #=====================================================================================================================================================================================================

        #=====================================================================================================================================================================================================
        ## Airfoil
        rho_cruise = density(cruise_height)
        reynold = airfoil.ReynoldNumber(chord, vc_cr, rho_cruise, 1.789*(10**-5))
        dynamic_press = airfoil.Dynam_press(rho_cruise, vc_cr)
        thick_chord = airfoil.t_over_c(Mach, MTOW, dynamic_press, wing_area, g)
        lift_wing = airfoil.Lwing_calc(MTOW, g)
        desired_lift_of_wing = airfoil.CLdes_calc(dynamic_press, WoverS)
        desired_lift_of_airfoil = airfoil.Cldes_calc(lift_wing, dynamic_press, wing_area)
        thickness_airfoil = airfoil.ThicknessAirfoil(thick_chord, chord)

        b = airfoil.SpanWing(wing_area, chord)


        print("Airfoil Sizing Complete")
        #=====================================================================================================================================================================================================

        #=====================================================================================================================================================================================================
        ##Class I Performance
        pre = os.path.dirname(os.path.realpath(__file__))
        fname = 'power_table.pickle'
        path = os.path.join(pre, fname)
        with open(path, 'rb') as handle:
            propdata = pickle.load(handle)
        aircraft = TiltRotor(propdata,cruise_height,350,270,220,MTOW,g,8,stall,A,wing_area,e,diskloading,1.15,0.75)
        aircraft.performFlight(1.225)
        # aircraft.plotFlight('time','powert')
        aircraft.roundtime()
        alex_data = aircraft.valuesForAlex()

        powers = pd.DataFrame(alex_data)


        # pre = os.path.dirname(os.path.realpath(__file__))
        # fname = 'mission_data.xlsx'
        # path = os.path.join(pre, fname)
        # print(pre)
        writer = pd.ExcelWriter('aerodynamics\mission_data.xlsx', engine='xlsxwriter')
        # for i, df in enumerate(dfs):
        # for pow in range(mission_power):
        powers.to_excel(writer, sheet_name='mission_data', index=False)

        writer.save()

        print("Class I Performance Complete")
        #=====================================================================================================================================================================================================

        #=====================================================================================================================================================================================================
        ## Rotor Sizing
        warnings.filterwarnings("ignore")
        radius = (b-fuselage_width)/2
        print('radius', radius)
        N, B, R, cutout, solidity, theta_tip, taper, airfoilr, be = rp.initVariables(radius)
        rotor = Rotor(N, B, R, cutout, taper, solidity, theta_tip, airfoilr, 0, be, 0)
        # twists, t_index = rp.genTwist(1*math.pi/180)
        # for twist in twists:
        rotor.calcTwist('linear',18.3*math.pi/180,0*math.pi/180)
        # rotor.calcTwist('ideal', 1*math.pi/180, 10)
        powers, times, hs, v_infs, angles, mission_profile, mission_time, dcts, omegas, collectives, inflows = rp.finalsim(alex_data,rotor,chord)
        # rp.calcBladeStructure(rotor, dcts, omegas, collectives, inflows)
        # print('Omegas', omegas)
        # print('Powers', powers)
        # print('Torques', powers/omegas*1000)

        # print('Omega range: ', min(omegas), max(omegas))
        # print('Power range: ', min(powers), max(powers))
        # print('Torque range: ', min(powers/omegas)*1000, max(powers/omegas)*1000)

        # print('rotor weight')
        # print(rp.calcRotorWeight(rotor))
        # fig, ax = plt.subplots(4,1)
        # ax[0].plot(times, hs)
        # ax[0].set_ylabel('hs')
        # ax[1].plot(times,powers)
        # ax[1].set_ylabel('powers')
        # ax[2].plot(times, v_infs)
        # ax[2].set_ylabel('v_infs')
        # ax[3].plot(times, angles)
        # ax[3].set_ylabel('angles')

        # print('profile')
        # print(mission_profile)
        # print(mission_time)

        # plt.figure(2)
        # plt.plot(mission_time, mission_profile)

        # plt.figure(3)
        # plt.plot(np.array(alex_data)[:,4], np.array(alex_data)[:,0])
        # plt.show()


        print("Rotor Sizing Complete")
        #=====================================================================================================================================================================================================

        #=====================================================================================================================================================================================================
        ## Cabin Sizing
        passenger_cabin = cabin.passenger_cabin(MTOW, OEW, fuselage_width, fuselage_height)
        cockpit_nose = cabin.cockpit_nose(MTOW, OEW, fuselage_width, fuselage_height)
        fixed_xcg, mass_cockpit_nose, mass_cabin, mass_fixed_equipment_cabin_cockpit,mass_fixed_equipment_cabin_cockpit_list = cabin.xcg_fixed_equipment(cockpit_nose,passenger_cabin)
        fixed = [fixed_xcg, mass_fixed_equipment_cabin_cockpit+175] #cg location and mass of fixed equipment needed for cg calculation
        mass_main_landing_gear, xcg_main_landing, mass_front_landing_gear, xcg_front_landing = passenger_cabin.landing_gear_calculation(xcg, lm, ln)
        nose = [xcg_front_landing, mass_front_landing_gear]
        landing = [xcg_main_landing, mass_main_landing_gear]

        print('Cabin Sizing Complete')
        #=====================================================================================================================================================================================================

        #=====================================================================================================================================================================================================
        ## Loading

        LiftCoefficients = ld.LiftCoefficients(cl_stall_sealevel, cl_cruise_sealevel, cl_dive_sealevel,cl_stall_cruise, cl_cruise_cruise,cl_dive_cruise,cl_max_flaperons) # assumption that cl dive and cruise are the same
        HorizontalSpeeds = ld.HorizontalSpeeds(vc_to, vc_cr)
        LimitLoadsCS29 = ld.LimitLoadsCS29()   
        Rotorcraft = ld.Rotorcraft(MTOW, lateral_drag_area, front_drag_area, cruise_height, Cl_alphaw, lt, lcg, wing_area, chord, vc_cr)
        ld.maneuver_diagram(HorizontalSpeeds, Rotorcraft, LimitLoadsCS29, LiftCoefficients, plot=False)
        GustSpeedCS25 = ld.GustSpeedCS25(Rotorcraft, vc_to, vc_cr, HorizontalSpeeds.vs_to,HorizontalSpeeds.vs_cr,Rotorcraft.lift_curve,Rotorcraft.wing_loading,Rotorcraft.chord,Rotorcraft.density_sea_level, Rotorcraft.density_cruise, 1.1) # speeds at sea level and cruise are the same
        delta_n_u = ld.gust_diagram(Rotorcraft,GustSpeedCS25, plot=False)
        maximum_load, vb_to_solution, vb_cr_solution = ld.combined_load_diagram(HorizontalSpeeds, Rotorcraft, LimitLoadsCS29, LiftCoefficients, delta_n_u, GustSpeedCS25, plot=False, plot_final=False)

        print('Loading Calculations Complete')
        #=====================================================================================================================================================================================================

        #=====================================================================================================================================================================================================
        ## Powerplant
        mass_powerplant, mass_fuelcell, mass_DCDC, mass_compressor, mass_radiator, mass_electric_motor, mass_motor_controller, mass_tank, mass_fuel, mass_battery, m_tankinnerwall, m_radshield, m_spacer, m_out, m_support = convertTuple(pp.calculateWeightFPP(mission_profile, mission_time, np.around(alex_data[7][4], decimals=1), np.around(alex_data[10][4], decimals=1), cond))
        tank_mass_list=[m_tankinnerwall, m_radshield, m_spacer, m_out, m_support]
        print('propusion system weight', mass_powerplant)
        print("Powerplant Characteristics Complete")
        #=====================================================================================================================================================================================================

        #=====================================================================================================================================================================================================
        ## Wing Structure
        # initialize values
        tilt_rotor_mechanism_mass = 40 #  is the mass of the tilting mechanism,
        gearbox_mass = 252 #  is the mass of the 2 gearboxes
        nacelle_weight = mass_electric_motor+mass_motor_controller+tilt_rotor_mechanism_mass+gearbox_mass
        # hover
        mass_nacelle_list=[mass_electric_motor,mass_motor_controller,gearbox_mass,tilt_rotor_mechanism_mass]
        phase = 'hover'
        M_h, Sy_h, T_h = ws.wingloads(phase,LimitLoadsCS29.upper_bound_m_l_cs_29,b,fuselage_width,MTOW,nacelle_weight,mass_radiator, maximum_load, lt, lcg, rho_cruise, vb_cr_solution, wing_area, chord, Cm_acw)
        wing_mass_h, mass_fractions_h, areas_h, _, _ = ws.wingstructure(M_h, Sy_h, T_h, b)

        # cruise
        # phase = 'cruise'
        # M_c, Sy_c, T_c = ws.wingloads(phase,LimitLoadsCS29.upper_bound_m_l_cs_29,b,fuselage_width,MTOW,nacelle_weight,mass_radiator, maximum_load, lt, lcg, rho_cruise, vb_cr_solution, S, chord, Cm_acw)
        # wing_mass_c, mass_fractions_c, areas_c = ws.wingstructure(M_c, Sy_c, T_c, b)

        # if wing_mass_h > wing_mass_c:
        #     wing_mass = wing_mass_h
        #     wing_mass_fractions = mass_fractions_h
        #     wing_areas = areas_h
        # else:
        #     wing_mass = wing_mass_c
        #     wing_mass_fractions = mass_fractions_c
        #     wing_areas = areas_c

        wing_mass = wing_mass_h

        print('wing weight', wing_mass)
        print('Wing Structure Complete')
        #=====================================================================================================================================================================================================

        #=====================================================================================================================================================================================================
        ## Fuselage


        e_fuselage = body.OswaldEfficiencyFActor(0.1907)
        CD0_fuselage = body.ZeroLiftDragCoeff(0.0045, 4)
        CLa_fuselage = body.FuselageLiftCurveSlope(0.1907, Mach, 0.95)
        CL_list_fuselage = body.CL_calc(CLa_fuselage, -20,20, -1.5)
        CD_list_fuselage = body.CD_calc(CD0_fuselage,0.1907,e_fuselage, vc_cr, temperature(cruise_height), 1.4, 287.15, 0.95, -20, 19, -1.5)
        # CL_optimum_fuselage = body.CLtrimcalc(CLa_fuselage, alpha_trim, -1.5)
        maxthick_fuselage = body.MaxThickness(0.21, fuselage_length)
        stall_angle_fuselage = body.StallAnglefuselage_calc()

        print("Fuselage Characteristics Complete")
        #=====================================================================================================================================================================================================

        #=====================================================================================================================================================================================================
        ## Stability and Control
        fus = [fuselage_length*0.5, fuselage_mass]
        cargo = (fuselage_length*0.2+3)*1.03
        battery = [0.3*fuselage_length,mass_battery]
        tank = [0.95*fuselage_length,mass_tank]
        fuelcell = [0.8*fuselage_length,mass_fuelcell+mass_DCDC+mass_compressor]
        radiator = [0,mass_radiator+15]
        nacelle = [0,nacelle_weight]
        blades_rotor_mass=rp.calcRotorWeight(rotor)[0]
        hubs_rotor_mass = rp.calcRotorWeight(rotor)[1]
        prop_mass_list=[blades_rotor_mass,hubs_rotor_mass]
        prop = [0,rp.calcRotorWeight(rotor)[2]+140]
        wing = [0,wing_mass]
        firstperson = passenger_cabin.xcg_passenger1
        seatpitch = passenger_cabin.xcg_passenger_pitch
        fuelpos = tank[1]
        oew_mass_budget=[fus[1],fixed[1],horizontal[1],vertical[1],nose[1],prop[1],nacelle[1],wing[1],landing[1],battery[1],tank[1],fuelcell[1],radiator[1]]
        cg, Sh_S, xLEMAC, tail_arm, OEW, MTOW = con_sta.runControlStabilitySizing(MAC,fuselage_length,fus,fixed,horizontal,vertical,nose,prop,nacelle,wing,landing,battery,tank,fuelcell,radiator,payload,passengermass,mass_fuel,cargo,firstperson,seatpitch,fuelpos,MAC,Ah,Lambda_ch,Cl_alphaw,wing_area,Cm_acf,Cm_acw,lh,b,Cl_f,Cl_alphaf,xac_w,SM,Cl_w,Sf,lf)
        wing_pos = xLEMAC*fuselage_length+chord*0.25
        prop[0] = xLEMAC*fuselage_length
        wing[0] = wing_pos
        radiator[0] = prop[0] + 0.5*MAC
        nacelle[0] = prop[0] + 0.5*MAC
        print('cg range as fraction of MAC from xLEMAC', cg)
        print('xLEMAC', xLEMAC)
        print('Sh/S', Sh_S)
        print('OEW', OEW)
        print('MTOW', MTOW)
        cgfor = xLEMAC*fuselage_length + cg[0]*MAC
        cgaft = xLEMAC*fuselage_length + cg[1]*MAC
        xcg = (cgfor+cgaft)/2
        landing[0] = xcg+lm
        print ('OEW mass budget', sum(oew_mass_budget))
        print ('OEW mass budget per sub-system', oew_mass_budget)

        print("Stability and Control Sizing Complete")
        #=====================================================================================================================================================================================================

        #=====================================================================================================================================================================================================
        ## Landing Gear Full Sizing
        LandingGear = lg.LandingGear(MTOW,3158, 0.19,0.47,0.8,0.06,0.47,0.0178,2,ln,lm,21,9,6.75,9.1,13.6,6,4.25,6.1)
        print("Landing Gear Sizing Complete")
        #=====================================================================================================================================================================================================
        MTOW_iterations.append(MTOW)
    MTOW_percentages.append(MTOW)

ns = []
for i in range(len(MTOW_iterations)):
    ns.append(i+1)
plt.plot(times, powers)
plt.show()

plt.plot(ns,MTOW_iterations)
plt.show()
# filename = 'optimumworkspace.out'
# my_shelf = shelve.open(filename, 'n')
# for key in dir():
#     try:
#         my_shelf[key] = globals()[key]
#     except TypeError:
#         #
#         # __builtins__, my_shelf, and imported modules can not be shelved.
#         #
#         print('ERROR shelving: {0}'.format(key))
# my_shelf.close()


# filename = 'globalsave.pkl'
# dill.dump_session(filename)

"""
To load shelf 
my_shelf = shelve.open(filename)
for key in my_shelf:
    globals()[key]=my_shelf[key]
my_shelf.close()
"""

session = {}
for key in dir():
    if type(globals()[key]) == type(float()) or type(globals()[key]) == type(int()) or type(globals()[key]) == type(list()) or type(globals()[key]) == type(np.float64()) or type(globals()[key]) == type(np.array) or type(globals()[key]) == type(mission_profile) or type(globals()[key]) == type(dict()):
    # if type(globals()[key]) != type(initialiseAerodynamics) or  type(globals()[key]) != type(math) or type(globals()[key]) != type(Rotor) or type(globals()[key]) != type(ws):
        session[key] = globals()[key]


with open('finaldata_fuselage_test.pickle', 'wb') as handle:
    pickle.dump(session, handle, protocol=pickle.HIGHEST_PROTOCOL)


