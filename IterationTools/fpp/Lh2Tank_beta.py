#Written by Killian Dally
#Cryogenic Pressure Vessel Design
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fmin
import math
import time
import random
import progressbar





#************************************************************************************************************************************************************************************************##

def designMechanicalTank(m_h2_flight,P_fc,list):



    m_h2_flight = m_h2_flight[0,:]
    P_fc = P_fc[0, :]

    ## inputs
    # requirements:
    p_vent = 2.5e5 # maximum allowable pressure in Pa, also operating pressure
    y_max = 0.97 # lh2 liquid volume fraction to be vented suddenly or supplied to the fuel cell
    p_fill = 1.2e5 # min. allowed pressure to avoid explosions according to Joty's paper
    p_atm = 1.2e5  # atm pressure in Pa

    # operating parameters
    empty_frac = 0.15 # minimum volume of fuel in the tank to be used in gaseous form to maintain design pressure according to TNO
    max_frac = 0.85 # max fuel volume to allow for liquid expansion and contingencies w/o rupture according to TNO

    # piping components
    t_pipe = 0.5e-3 # guesses

    # shape parameters
    lamb = 3. # = l_s / r_tank sphere
    vacuum = True #deciding if using vacuum insulation

    # allowances
    c_1 = 0.009 # tank contraction due to cooling and expansion according to Brewer book
    c_2 = 0.006 # equipment inside tank e.g. pumps according to Brewer book
    c_3 = 0.05 # gas to be vented in flight according to Brewer book
    c_4 = 0. # corrosion allowance
    c_5 = 0. # thinning allowance
    weld_eff = 0.8 # welding efficiency according to ASME

    # material characteristics                                                                              # http://www.gb.nrao.edu/electronics/edir/edir306.pdf
    # SF (1.4 - 2) safety factor according to "Review of current state of the art and key design issues with potential solutions for liquid hydrogen cryogenic storage tank structures for aircraft applications."
    alum_2024 = {'poisson': 0.337, 'E_mod': 72.0e9, 'stress_limit': 176e6, 'rho_mat': 2855,'emissivity': 0.7,"SF":2}  # E= 81GPa at cryogenic temp. acc https://books.google.nl/books?hl=en&lr=&id=3U_eQdnmzmwC&oi=fnd&pg=PA1&dq=fatigue+strength+aluminum+2024+low+temperature&ots=1j1Wdi3IQ9&sig=fU2VrLHANov0t7U8x8UhVl0Zs8E#v=onepage&q=fatigue%20strength%20aluminum%202024%20low%20temperature&f=false
    CFRP = {'poisson': 0.27, 'E_mod': 44.2e9, 'stress_limit': 478e6, 'rho_mat': 1540,'emissivity': 0.875,"SF":3}  # QI woven-fabric and https://link.springer.com/article/10.1007/s10765-019-2498-0 and S-N Curve Models for Composite Materials Characterisation: An Evaluative Review
    stain_steel = {'poisson': 0.27, 'emissivity': 0.44,'E_mod': 190e9,"thermal_cond":8.425,'rho_mat':7850,"SF":2,'stress_limit': 550e6,} # Stainless steel austenitic (304) with thermal condctivity average between cold and hot from boom Cryo. Eng. explained and http://www-eng.lbl.gov/~dw/projects/DW4229_LHC_detector_analysis/calculations/emissivity2.pdf
    # spray_on_foam = {'thermal_cond': 0.005,'rho_mat': 38.4443}  # NCFI 24-124 according to "Thermal conductivity of rigid foam insulations for aerospace vehicles" and https://www.nasa.gov/pdf/63758main_TPS_FACT_SHEET.pdf
    MLI = {'emissivity': 0.031,'rho_mat': 140,"N_density":20,"SF":2}  # at 50^C from https://iopscience.iop.org/article/10.1088/1757-899X/101/1/012017/pdf and https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/20110014015.pdf, radiation shield is  = aluminized mylar film = polymer with t=0.006mm
    #and # layers/cm  from https://www.researchgate.net/publication/283072749_Cryogenic_systems_in_Radio_Astronomy and https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/20110014015.pdf
    #compare to aerogel https://www.aerogel.com/products-and-solutions/product-documents/
    G10 = {"thermal_cond":0.288, 'rho_mat': 1800,"SF":2} # https://laminatedplastics.com/g-10.pdf


    ## structural Design

    def designOpeningReinforcement(d, t_s, t_n, wall_mat, pipe_mat):

        # this function design opening reinforcements for any hole and wall size
        # works only on discrete values because of the optimization required

        #guesses
        T_n_option = np.arange(1.2e-3,4e-3,0.5e-4)
        T_s = 2e-3 + t_s
        m_reinf = np.array(d.shape)
        reinf_options = np.ones((d.shape[0],T_n_option.shape[0],5))
        reinf_check = np.ones((d.shape[0], T_n_option.shape[0]))


        for j in range(len(T_n_option)):

            # geometrical parameters
            T_n = T_n_option[j] * np.ones(T_s.shape)
            R_n = d / 2.
            x = np.maximum(d,R_n + t_n + T_n)
            y = np.minimum(2.5 * T_s, 2.5 * T_n)
            A_required = d * t_s
            A_s = np.maximum(d * (T_s - t_s) - 2 * T_n * (T_s - t_s), 2 * (T_s + t_n) * (T_s - t_s) - 2 * t_n * (T_s - t_s))
            A_n = np.minimum(2 * (2.5 * T_s * (T_n - t_n)), 2 * (2.5 * T_n * (T_n - t_n)))
            m_reinf_pipe = pipe_mat["rho_mat"] * np.pi * ((d/2 + T_n)**2 - ((d/2) + t_n)**2) * (y - T_s + t_s)
            m_reinf_tank = wall_mat["rho_mat"] * (T_s - t_s) * (x - d/2 - t_n)
            m_reinf = m_reinf_pipe + m_reinf_tank
            Reinforcement = np.maximum(((A_s + A_n)-A_required),np.zeros(A_s.shape))
            Reinforcement[Reinforcement == 0] = 10000

            reinf_options[:, j, 0] = T_n
            reinf_options[:, j, 1] = T_s
            reinf_options[:, j, 2] = x
            reinf_options[:, j, 3] = y
            reinf_options[:, j, 4] = m_reinf
            reinf_check[:, j] = Reinforcement

        good_design = reinf_check.argmin(1)

        for i in range(d.shape[0]):

            T_n[i] = reinf_options[i,good_design[i],0]
            T_s[i] = reinf_options[i, good_design[i], 1]
            x[i]   = reinf_options[i, good_design[i], 2]
            y[i] = reinf_options[i, good_design[i], 3]
            m_reinf[i] = reinf_options[i, good_design[i], 4]

        return (T_n, T_s, x, y, m_reinf)


    def designInnerWall(material1):

        def calculateLh2Required(m_h2_flight):

            # calculate the required LH2 mass
            return m_h2_flight / (max_frac - empty_frac)  / (1. - c_3) #return m_h2_rq


        def calculateLh2MassToTankVolume(m_h2,p,y_max):

            # calculate H2 tank volume based on h2 mass, internal pressure and liquid volume fraction
            rho_lh2 = np.interp(p,[1.,10.98],[70.83,47.29])
            rho_gh2 = np.interp(p, [1., 10.73], [1.427, 15.898])
            rho_mean = y_max * rho_lh2 + (1.-y_max) * rho_gh2
            #print("v_spe =",1/(rho_mean/1000))
            return m_h2/rho_mean/(1 - c_1 - c_2) #return h2 volume = V_h2


        def calculateVolumeToRadius(V_h2,lamb):

            # calculate sphere radius
            return ( V_h2 / np.pi / (4. / 3. + lamb))**(1. / 3.) # = r_in_1 = inside (1) radius of inner wall


        def calcualteCylindricalWallThickness(p_vent,r_in_1,material):

            # calculate cylindrical thickness
            p_burst = material1['SF'] * p_vent  # burst pressure at which the tank fails according to Brewer book, w/o safety factor
            return (p_burst * r_in_1) / (weld_eff * material['stress_limit'] - 0.6 * p_burst) + c_4 + c_5 # = s_w_cyl


        def calcualteSphericalWallThickness(p_vent,r_in_1):

            # calculate sphere thickness
            p_burst = material1['SF'] * p_vent  # burst pressure at which the tank fails according to Brewer book, w/o safety factor
            return (p_burst * r_in_1) / (2. * weld_eff * material1['stress_limit']  - 0.2 * p_burst) + c_4 + c_5 # = s_w_sphe


        def calculateInnerWallMass(t_in_cyl,t_in_sphe,r_in_2_cyl,r_in_2_sphe):

            # calculate inner wall mass
            return (4. / 3. * np.pi * (r_in_2_sphe** 3 - r_in_1**3) + np.pi * (r_in_2_cyl**2. - r_in_1**2.) * lamb * r_in_1)* material1['rho_mat'] # = m_tank


        # combine all fucntions together
        # subscript 1 and 2 refer to interior and exterior dimensions respectively
        m_h2_rq = calculateLh2Required(m_h2_flight)
        r_in_1 = calculateVolumeToRadius(calculateLh2MassToTankVolume(m_h2_rq, p_fill, y_max), lamb)
        m_liner = 0
        if material1 == CFRP:

            t_liner_al = 0.635e-3 # https://apps.dtic.mil/dtic/tr/fulltext/u2/b211395.pdf
            r_in_1 = r_in_1 + t_liner_al
            m_liner = (4. / 3. * np.pi * (r_in_1** 3 - (r_in_1 - t_liner_al)**3) + np.pi * (r_in_1**2. - (r_in_1 - t_liner_al)**2.) * lamb * r_in_1)* alum_2024['rho_mat']


        t_in_cyl = calcualteCylindricalWallThickness(p_vent, r_in_1,material1)
        t_in_sphe = calcualteSphericalWallThickness(p_vent, r_in_1)
        r_in_2_cyl = r_in_1 + t_in_cyl #radius cylidner including thickness
        r_in_2_sphe = r_in_1 + t_in_sphe #radius spherical end including thickness
        m_in = calculateInnerWallMass(t_in_cyl,t_in_sphe,r_in_2_cyl,r_in_2_sphe) + m_liner

        # print("Wall cyl. thickness = ", t_in_cyl * 1000, "mm")
        # print("Wall sph. thickness = ", t_in_sphe * 1000, "mm")
        # print("Tank mass           = ", m_in, "kg")
        # print("Radius              = ", r_in_2_cyl, "m")
        # print("Rq. fuel mass       = ", m_h2_rq, "kg")

        return(m_in,r_in_1,r_in_2_cyl,r_in_2_sphe,m_h2_rq,t_in_cyl)


    def designThermalInsualtion(material1,material2,material3,material4,material5,material6,m_h2_rq,P_fc,vacuum,t_spaceneeded,t_added_vent, t_added_outlet, t_ins_old):

        # general inputs
        T_amb = 50. + 273.15  # K
        T_lh2 = 23.5  # K
        fc_eff = 0.4477
        SE_h2 = 141.8e6  # specific energy lh2
        h_f_h2 = 461.e3  # latent heat of vaporization in J/kg
        t_max_flight = 60 * 60.  # 1hr
        t_turnaround = 60 * 60.  # 1hr
        vent_ratio = c_3
        rho_lh2 = 70 #kg/m3
        rho_gh2 = 2.7  # kg/m3
        gamma = 1.405 #cp/cv from http://catalog.conveyorspneumatic.com/Asset/FLS%20Specific%20Heat%20Capacities%20of%20Gases.pdf
        gas_cte_h2 = 4124
        m_lh2_in = 60e-3  # kg/s from https://www.linde-engineering.com/en/images/DS_Cryo%20Pump_tcm19-523716.pdf
        v_lh2 = 20 # https://www.engineeringtoolbox.com/fluid-velocities-pipes-d_1885.html

        # geometrical inputs
        t_vac = t_ins_old + t_spaceneeded
        if t_ins_old[0] == 0:
            t_vac = 1e-2 + t_spaceneeded

        [m_in, r_in_1, r_in_2_cyl, r_in_2_sphe, m_h2_rq, t_in_cyl] = designInnerWall(material1)
        A_in_2 = r_in_2_cyl * 2. * np.pi * lamb * r_in_1 + 4. * np.pi * r_in_2_sphe**2
        A_in_2_cyl = r_in_2_cyl * np.pi * 2 * lamb * r_in_1
        A_out_1_cyl = (r_in_2_cyl + t_vac) * np.pi * 2 * lamb * r_in_1
        A_in_2_sphe = 4 * np.pi * r_in_2_sphe ** 2
        A_out_1_sphe = 4 * np.pi * (r_in_2_sphe + t_vac) ** 2
        A_out_2 = (r_in_2_cyl + t_vac) * 2. * np.pi * lamb * r_in_1 + 4. * np.pi * (r_in_2_sphe + t_vac) ** 2
        D_outlet = (1e-2 + t_added_outlet * 2) * np.ones(m_in.shape)  # diameter given by components, operations...

        # conduction inputs
        t_onelayer = (0.03 + 0.02) * 1e-2 #Cryogenic Heat Transfer
        t_onelayer_al = 0.006e-3 #https://books.google.nl/books?id=8fZ6-Nl491EC&printsec=frontcover&dq=Cryogenic+Heat+Transfer+De+Randall+F.+Barron&hl=fr&sa=X&ved=0ahUKEwj2iIjch-TiAhXKJVAKHTrrA90Q6AEIKDAA#v=onepage&q=aluminized%20polyester%20film&f=false
        t_onelayer_fg = 0.015e-3 #https://books.google.nl/books?id=8fZ6-Nl491EC&printsec=frontcover&dq=Cryogenic+Heat+Transfer+De+Randall+F.+Barron&hl=fr&sa=X&ved=0ahUKEwj2iIjch-TiAhXKJVAKHTrrA90Q6AEIKDAA#v=onepage&q=aluminized%20polyester%20film&f=false
        h_s = 0.0851 #solid conductance glass paper from "Cryogenic Engineering, Revised and Expanded"
        A_m_cyl = (A_out_1_cyl - A_in_2_cyl) / np.log(A_out_1_cyl / A_in_2_cyl)
        A_m_sphe = (A_in_2_sphe * A_out_1_sphe) ** (1. / 2.)
        n_support = 8.
        EnoughSpace = True

        # radiation inputs
        SB_cte = 5.6703e-8  # (W/m2K4)
        View_factor = 1.  # https://www.dspe.nl/knowledge-base/thermomechanics/chapter-1---basics/1-2-heat-transfer/radiation/

        # mission requirement calcualtions
        m_dot_lh2_fc = P_fc / SE_h2 / fc_eff
        Q_need = m_dot_lh2_fc * h_f_h2
        A_outlet = m_dot_lh2_fc / rho_lh2 / v_lh2
        D_outlet_rq = ((A_outlet / np.pi) ** 0.5 + t_pipe * 2) #diameter we need for the FC

        E_t = m_h2_rq * h_f_h2
        Q_allow = 0.025 * E_t / (t_turnaround)

        #print(Q_allow)
        m_dot_gh2_vent = Q_allow / h_f_h2
        #print(0.03 * m_h2_rq/m_dot_gh2_vent/60.)


        Mach = 0.1
        v_gh2 = (gamma * gas_cte_h2 * T_lh2 / Mach )**0.5  #speed of gh2 through vent
        A_vent = m_dot_gh2_vent / rho_gh2 / v_gh2
        D_vent_rq = (material6['SF']*(A_vent / np.pi)**0.5 + t_pipe*2)
        D_vent = (6e-3 + 2 * t_pipe + 2 * t_added_outlet) * np.ones(A_m_cyl.shape) #https://www.herose.com/eng/products/cryogenic-services/06001_dt.php

        #print(m_dot_gh2_vent)

        A_inlet = m_lh2_in / rho_lh2 / v_lh2
        D_inl = ((A_inlet / np.pi) ** 0.5 + t_pipe * 2)
        # Check D_inl is smaller than D_outlet

        # conduction initialization
        k_a = 1 / (material4["N_density"] * 100) * (h_s + (SB_cte * material4["emissivity"] * T_amb**3)/(2 - material4["emissivity"]) * (1 + (T_lh2 / T_amb)**2) * (1 + (T_lh2 / T_amb)))
        A_support = np.pi /4. * (0.04**2 - 0.025**2)

        A_pipe = np.pi * ((D_outlet/2)**2 - ((D_outlet/2) - t_pipe)**2) + np.pi * ((D_vent/2)**2 - ((D_vent/2) - t_pipe)**2)

        guess = 15

        n_layer = np.array([])
        for j in range(A_pipe.size):
            def findMinimumLayersMLI(x):  # minimum thickness finder

                Q_cond_pipe = material6["thermal_cond"] * (T_amb - T_lh2) * A_pipe[j] / (t_onelayer * x[0])

                Q_cond_support = material5["thermal_cond"] * (T_amb - T_lh2) * A_support * n_support / (
                        t_onelayer * x[0])
                Q_cond_MLI = (T_amb - T_lh2) * k_a * (A_m_cyl[j] + A_in_2_sphe[j]) / (t_onelayer * x[0])
                emissivity_fac = 1. / (1. / material1["emissivity"] + 1. / material4["emissivity"] - 1 + (x[0] - 1) * (
                        2. / material4["emissivity"]) + 1. / material4["emissivity"] + 1. / material2[
                                           "emissivity"] - 1)
                # according to "Cryogenic Heat Transfer" from De Randall F. Barron

                Q_rad = View_factor * SB_cte * A_in_2[j] * (T_amb ** 4. - T_lh2 ** 4.) * emissivity_fac

                Q_cond = Q_cond_support + Q_cond_pipe

                orig = Q_cond_support + Q_cond_MLI + Q_rad + Q_cond_pipe

                return (abs(orig - Q_allow[j]))

            n_layer = np.append(n_layer,math.ceil(fmin(findMinimumLayersMLI, guess, disp=False,ftol=0.1)[0]))

        t_ins = t_onelayer * n_layer

        m_pipe = A_pipe * t_ins * material6['rho_mat']
        m_support = A_support * t_ins * material5['rho_mat']
        m_MLI = (4. / 3. * np.pi * ((r_in_2_sphe + t_ins) ** 3 - (r_in_2_sphe)  ** 3.) + np.pi * ((r_in_2_cyl + t_ins) ** 2 - (r_in_2_cyl) ** 2) * lamb * r_in_1) * material4['rho_mat']


        t_ins += np.maximum(t_spaceneeded-t_ins_old,0)

        # print("n_layer =", n_layer)
        # print("t_ins =", t_ins)

        return(m_MLI,t_ins,m_pipe,m_support,D_vent,D_outlet)


    def designOuterWall(material1,material2,t_ins):

        p_c = material2['SF'] * p_atm  # collapse pressure

        # retrieve previous data
        # subscript 1 and 2 refer to interior and exterior dimensions respectively

        r_in_1 = designInnerWall(material1)[1]
        r_out_1_cyl = designInnerWall(material1)[2] + t_ins
        r_out_1_sphe = designInnerWall(material1)[3] + t_ins

        length_cyl = r_in_1 * lamb  # length of cylindrical part
        length_out_1 = r_out_1_sphe * 2. + length_cyl #length including spherical ends
        d_out_1_cyl = r_out_1_cyl * 2.


        def calculateCylindricalOuterWallThickness():

            # calculate cylindrical outer wall thickness
            t_out_cyl = np.array([])
            for j in range(d_out_1_cyl.size):

                def findMinimumCollapsePressureTchinessCyl(x): #minimum thickness finder

                    orig = (2.42 * material2['E_mod'] *  ((x[0] / d_out_1_cyl[j]) ** (5. / 2.))) / (((1. - material2['poisson'] ** 2.) ** (3. / 4.)) * ((length_cyl[j] / d_out_1_cyl[j]) - 0.45 * ((x[0] / d_out_1_cyl[j]) ** (1. / 2.))))
                    #according to https://pdfs.semanticscholar.org/1119/752b11d766cbd15ffb699ef89406110d52ff.pdf
                    #supported ends thanks to supports

                    #orig = 0.855 * material2['E_mod'] * 0.75 / ((1. - material2['poisson'] ** 2.) ** (3. / 4.) * (r_out_1_cyl[j] / x[0])**(5./2.) * (length_cyl[j] / r_out_1_cyl[j])) # accoridng to https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19690013955.pdf

                    return(abs(orig - p_c))

                guess = 0.003
                t_out_cyl = np.append(t_out_cyl,fmin(findMinimumCollapsePressureTchinessCyl, guess, disp=False,ftol=0.00005)[0] )

            return(t_out_cyl)


        def calculateSphericalOuterWallThickness():

            # calculate spherical outer wall thickness according to https://www.ijedr.org/papers/IJEDR1501040.pdf
            return r_out_1_sphe * (p_c * (3. * (1. - material2['poisson'] ** 2.)) ** (1. / 2.) / (0.5 * material2['E_mod'])) ** (1. / 2.) # = t_out_sphe


        def calculateOuterWallMass(r_out_2_cyl,r_out_2_sphe):

            #calculate outer wall mass
            return (4. / 3. * np.pi * (r_out_2_sphe** 3 - r_out_1_sphe**3.) + np.pi * (r_out_2_cyl**2 - r_out_1_cyl**2) * lamb * r_in_1)* material1['rho_mat'] # = m_tank


        t_out_cyl = calculateCylindricalOuterWallThickness()
        t_out_sphe = calculateSphericalOuterWallThickness()
        r_out_2_cyl = r_out_1_cyl + t_out_cyl
        r_out_2_sphe = r_out_1_sphe + t_out_sphe
        m_out = calculateOuterWallMass(r_out_2_cyl,r_out_2_sphe)

        T_amb = 50. + 273.15  # K
        T_lh2 = 23.5  # K
        coeff_expansion = 24e-6 #from http://hyperphysics.phy-astr.gsu.edu/hbase/Tables/thexp.html
        delta_expansion =  coeff_expansion * (T_amb -  T_lh2)

        # print("Wall cyl. thickness = ", t_out_cyl * 1000, "mm")
        # print("Wall sph. thickness = ", t_out_sphe * 1000, "mm")

        return(m_out,r_out_2_cyl,r_out_2_sphe,t_out_cyl)


    material1 = CFRP
    material2 = alum_2024
    material3 = None
    material4 = MLI
    material5 = G10
    material6 = stain_steel

    bar = progressbar.ProgressBar(maxval=100, widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
    bar.start()

    [m_in,r_in_1,r_in_2_cyl,r_in_2_sphe,m_h2_rq,t_in_cyl] = designInnerWall(alum_2024)
    bar.update(5)
    [m_ins,t_ins,m_pipe,m_support,D_vent,D_outlet] = designThermalInsualtion(material1, material2, material3, material4,material5, material6, m_h2_rq, P_fc, vacuum, np.zeros(m_in.shape),np.ones(m_in.shape)*0.0028,np.ones(m_in.shape)*0.0028,np.zeros(m_in.shape))
    bar.update(10)
    [m_out,r_out_2_cyl,r_out_2_sphe,t_out_cyl] = designOuterWall(material1, material2, t_ins)
    bar.update(20)


    [T_n_1, T_s_in_vent, x, y_in_vent, m_pipe_vent_in] = designOpeningReinforcement(D_vent,t_in_cyl,t_pipe,material1,material6)
    bar.update(25)
    [T_n_2, T_s_out_vent, x, y_out_vent, m_pipe_vent_out] = designOpeningReinforcement(D_vent,t_out_cyl,t_pipe,material1,material6)
    bar.update(35)
    [T_n_3, T_s_in_outlet, x, y_in_outlet, m_pipe_outlet_in] = designOpeningReinforcement(D_outlet,t_in_cyl,t_pipe,material1,material6)
    bar.update(5)
    [T_n_4, T_s_out_outlet, x, y_out_outlet, m_pipe_outlet_out] = designOpeningReinforcement(D_outlet,t_out_cyl,t_pipe,material1,material6)
    bar.update(5)

    t_spaceneeded_vent = (t_ins - y_in_vent - y_out_vent)
    t_spaceneeded_outlet = (t_ins - y_in_outlet - y_out_outlet)
    t_spaceneeded = -np.minimum(t_spaceneeded_vent,t_spaceneeded_outlet)

    t_added_vent = np.maximum(T_n_1,T_n_2)
    t_added_outlet = np.maximum(T_n_3,T_n_4)



    for i in range(2):
        bar.update(i*10+1)

        [m_ins, t_ins, m_pipe, m_support, D_vent,D_outlet] = designThermalInsualtion(material1, material2, material3, material4,material5, material6, m_h2_rq, P_fc, vacuum, t_spaceneeded, t_added_vent, t_added_outlet, t_ins)
        [m_out, r_out_2_cyl, r_out_2_sphe, t_out_cyl] = designOuterWall(material1, material2, t_ins)

        [T_n_2, T_s_out_vent, x, y_out_vent, m_pipe_vent_out] = designOpeningReinforcement(D_vent, t_out_cyl, t_pipe,material1, material6)
        [T_n_4, T_s_out_outlet, x, y_out_outlet, m_pipe_outlet_out] = designOpeningReinforcement(D_outlet, t_out_cyl, t_pipe, material1, material6)

        t_spaceneeded_vent = (t_ins - y_in_vent - y_out_vent)
        t_spaceneeded_outlet = (t_ins - y_in_outlet - y_out_outlet)
        t_spaceneeded = -np.minimum(t_spaceneeded_vent, t_spaceneeded_outlet)
        # print(t_spaceneeded)

    m_pipe =+ (m_pipe_vent_in + m_pipe_vent_out + m_pipe_outlet_in + m_pipe_outlet_out)
    m_tot = m_in + m_ins + m_out + m_support + m_pipe


    # def calculatePipeThickness(p_vent, r_in_1, material):
    #
    #     # calculate cylindrical thickness
    #     p_burst = material[
    #                   'SF'] * p_vent  # burst pressure at which the tank fails according to Brewer book, w/o safety factor
    #     return (p_burst * r_in_1) / (weld_eff * material['stress_limit'] - 0.6 * p_burst)
    # # Check
    # print("min_thick",calculatePipeThickness(p_vent,D_vent[i]/2,stain_steel))
    # print("min_thick", calculatePipeThickness(p_vent, D_outlet, stain_steel))

    # print("")
    # print("m_pipe =", m_pipe)
    # print("m_ins =", m_ins)
    # print("m_tot =", m_tot)
    #print("m_tank/m_h2 = ",m_tot/m_h2_rq)

    m_tot = np.array([m_tot])
    m_h2_rq = np.array([m_h2_rq])
    print(m_h2_rq/m_h2_flight)

    return(m_tot,m_h2_rq)


def main():
    pass

if __name__ == "__main__":

    # mission inputs

    m_h2_flight = np.arange(15,25,0.01)
    m_h2_flight = np.array([m_h2_flight])
    P_fc = [550e1]
    P_fc = np.array([P_fc])

    start_time = time.time()


    designMechanicalTank(m_h2_flight,P_fc,1)
    print("--- %s seconds ---" % (time.time() - start_time))


    def designFuelLine():


        return












    # m_tank = []
    # lamb_range = np.arange(0.,3.0,0.1)
    # for i in lamb_range:
    #     m_tank.append(designMechanicalTank(m_h2_flight,P_fc,i)[0]/designMechanicalTank(m_h2_flight,P_fc,i)[1])
    #
    # plt.plot(lamb_range,m_tank)
    # plt.show()


    #QUESTIONS
    #-welding
    #-fatigue for compression, external pressure vessel
    #t_liner
