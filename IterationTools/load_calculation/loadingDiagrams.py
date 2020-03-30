'''Code written by Gianluca Mancini'''

# import the necessary libraries
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import math
from sympy.solvers import solve
from sympy import Symbol
from scipy.interpolate import interp1d
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})


# establish the load limits: regulated limit loads, remote loads, maximum gust horizontal speed

class LimitLoadsCS29():

    """
    This class contains all the limits loads and
    horizontal gust speed that belong to the regulation
    CS29
    """

    def __init__(self):

        # establish limit loads

        self.unity_load = 1 # - (L/W)
        self.lower_bound_m_l =  -1  # - (L/W)
        self.upper_bound_m_l = 3.4   # - (L/W)  CS-23 limit based on weight relationship
        self.landing=1.15 # Roskam value for propeller aircraft

        self.upper_bound_m_l_cs_29=3.5 # - (L/W)  CS-29

        # establish remote load cases ( limits for loads which are reached only in extreme cases )

        self.lower_bound_r_l = -0.5   # - (L/W)
        self.upper_bound_r_l = 2   # - (L/W)

        # maximum gust horizontal speed

        self.gust_horizontal_speed = 9   # m/s

    def load_factors(self):

        return np.array([0,self.unity_load, self.landing, self.upper_bound_m_l,  self.upper_bound_m_l,0, self.lower_bound_m_l, self.lower_bound_m_l,0])


class Rotorcraft():

    """
    This class includes all the characteristics of the aicraft at MTOW for sea-level
    and cruise characteristics. It also includes the lift curve slope of the wing, the vtoss and take off distance
    """

    def __init__(self,MTOW,lateral_drag_area, front_drag_area, cruise_height, lift_curve, lt, lcg, wing_area, chord,cruise_speed):

        self.MTOW = MTOW  # kg
        self.lateral_drag_area = lateral_drag_area  # m^2
        self.front_drag_area = front_drag_area # m^2
        self.drag_coefficient = 1  # - equivalent flat plat area coefficient
        self.cruise_height = cruise_height  # kg/m^3
        self.lift_curve = lift_curve*180/math.pi #
        self.lt = lt # m (distance between main wing aerodynamic center and horizontal tail)
        self.lcg= lcg # m (distance between main wing aerodynamic center and center of gravity of plane)
        self.gamma=1.4 # gamma gas constant
        self.R = 287.05 # R gas constant
        self.temperature_gradient= -0.0065 # K/m
        self.temperature_sea_level = 288.15 # K
        self.density_sea_level = 1.225 # kg/m^3
        self.temparature_cruise= self.temperature_sea_level+self.temperature_gradient*cruise_height #K
        self.pressure_sea_level= 101325 # Pa
        self.g= 9.80665 #m/s^2
        self.pressure_cruise = self.pressure_sea_level*math.pow(self.temparature_cruise / self.temperature_sea_level, (-self.g / (self.temperature_gradient * self.R)))
        self.density_cruise= self.pressure_cruise/ (self.R * self.temparature_cruise)
        self.sound_speed_cruise=math.sqrt(self.gamma * self.R * self.temparature_cruise)
        self.cruise_speed=cruise_speed
        # self.lift_curve_cruise = lift_curve/(math.sqrt(1-self.cruise_speed/self.sound_speed_cruise))  # -
        self.vtoss = 15 #m/s speed derived from the constraint for vtoss at 61m altitude and horizontal take
        self.MTOW_weight=MTOW*self.g # N
        self.wing_area = wing_area #m^2
        self.wing_loading = self.MTOW_weight/wing_area # N/m^2
        self.chord=chord #m



class GustSpeedCS25():

    """
    This class contains all the gust speed loads
    for the operative altitude (below 6000m) at v cruise,
    v dive and v b (speed at high angle of attack). The relavent
    gust speed are expressed with the letter 'u'. They are needed to produce the gust loading diagram.
    """

    def __init__(self, Rotorcraft, vc_to, vc_cr, vs_to, vs_cr,lift_curve, wing_loading, chord, rho_sea_level, rho_cruise, safety_factor):

        # gust speed as set by cs25 regulations

        self.ub = 20 # m/s
        self.uc = 16 # m/s
        self.ud = 8 # m/s

        self.g = 9.80665  # m/s^s gravity constant
        self.kg_to_lbs=2.20462 # conversion from kg to lbs
        self.N_to_lbs=2.20462/9.80665 # conversion from N to lbs
        self.m_to_ft=1./0.3048 # conversion from m to ft
        self.kg_to_slug=0.0685218 # conversion from kg to slug
        self.m_s_to_knots=1.94384 # conversion from m/s to knots

        # calculation of gust loads parameters

        mu_sea_level_gust = Rotorcraft.MTOW / (
                    Rotorcraft.wing_area * 0.5 * Rotorcraft.density_sea_level * Rotorcraft.lift_curve * Rotorcraft.chord)
        Kg_sea_level_gust = (0.88 * mu_sea_level_gust) / (5.3 + mu_sea_level_gust)

        self.mu_sea_level_gust=mu_sea_level_gust
        self.Kg_sea_level_gust=Kg_sea_level_gust

        mu_cruise_gust = Rotorcraft.MTOW / (
                    Rotorcraft.wing_area * 0.5 * Rotorcraft.density_cruise * Rotorcraft.lift_curve * Rotorcraft.chord)
        Kg_cruise_gust = (0.88 * mu_cruise_gust) / (5.3 + mu_cruise_gust)

        self.mu_cruise_gust=mu_cruise_gust
        self.Kg_cruise_gust=Kg_cruise_gust


        # horizontal velocities at take off

        mu_to=2*(wing_loading*self.N_to_lbs/(self.m_to_ft)**2)/(chord*self.m_to_ft*(rho_sea_level*self.kg_to_slug/(self.m_to_ft)**3)*lift_curve*self.g*self.m_to_ft)

        Kg_to=0.88*mu_to/(5.3+mu_to)

        self.vb_to=(vs_to*math.sqrt(1+(Kg_to*self.ub*self.m_to_ft*vc_to*self.m_s_to_knots*lift_curve)/(498*(wing_loading*self.N_to_lbs/(self.m_to_ft)**2))))#/self.m_s_to_knots

        # self.vb_to = vc_to-22.12 # m/s # Roskam V relationship

        self.vc_to = vc_to # m/s
        self.vd_to = vc_to*1.25 # m/s # Roskam V relationship

        # horizontal velocities at cruise

        mu_cr=2*(wing_loading*self.N_to_lbs/(self.m_to_ft)**2)/(chord*self.m_to_ft*(rho_cruise*self.kg_to_slug/(self.m_to_ft)**3)*lift_curve*self.g*self.m_to_ft)

        Kg_cr=0.88*mu_cr/(5.3+mu_cr)

        self.vb_cr=(vs_to*math.sqrt(1+(Kg_cr*self.ub*self.m_to_ft*vc_cr*self.m_s_to_knots*lift_curve)/(498*(wing_loading*self.N_to_lbs/(self.m_to_ft)**2))))#/self.m_s_to_knots

        # self.vb_cr = vc_cr-22.12 # m/s # Roskam V relationship

        self.vc_cr = vc_cr # m/s
        self.vd_cr = vc_cr*1.25 # m/s # Roskam V relationship

        # establish a safety factor since only two altitude conditions at MTOW are considered to find the limit loads

        self.safety_factor = safety_factor # -

class LiftCoefficients():

    """
    This class contains the lift coefficient needed to find the different speeds to produce the maneuver load diagram.
    The lift coefficient are the same for take off and cruise diagrams
    """

    def __init__(self, cl_stall_sealevel, cl_cruise_sealevel, cl_dive_sealevel,cl_stall_cruise, cl_cruise_cruise,cl_dive_cruise, cl_stall_full_flaps):

        self.cl_stall_sealevel = cl_stall_sealevel  # -
        self.cl_cruise_sealevel = cl_cruise_sealevel  # -
        self.cl_dive_sealevel = cl_dive_sealevel  # -

        # in case the cl would not be different
        self.cl_stall_cruise = cl_stall_cruise  # -
        self.cl_cruise_cruise = cl_cruise_cruise  # -
        self.cl_dive_cruise = cl_dive_cruise  # -

        self.cl_stall_full_flaps=cl_stall_full_flaps # - cl stall in full flaps condition


class HorizontalSpeeds():

    """
    This class contains the horizontal speed needed to produce the maneuver load diagram. The dive speed
    """

    def __init__(self,vc_to,vc_cr):

        # take off diagram
        self.vs_to=1 # m/s, stall speed take off
        self.va_to = 1  # m/s, va take off
        self.vc_to = vc_to # m/s
        self.vd_to = vc_to*1.25 # m/s # Roskam V relationship
        self.H_to = 1 # m/s

        self.vs_full_flaps=1 #m/s

        # cruise diagram
        self.vs_cr= 1 # m/s, stall speed cruise
        self.va_cr = 1  # m/s, va cruise
        self.vc_cr = vc_cr # m/s
        self.vd_cr = vc_cr*1.25 # m/s # Roskam V relationship
        self.H_cr = 1  # m/s

    def to_speeds(self):
        return np.array([0,self.vs_to,self.vs_full_flaps, self.va_to, self.vd_to,self.vd_to, self.vc_to, self.H_to,0])

    def cr_speeds(self):
        return np.array([0,self.vs_cr, self.va_cr, self.vd_cr,self.vd_cr, self.vc_cr, self.H_cr,0])



class Checks():

    def __init__(self, v_sound_cruise, v_sound_dive, vs_to,va_to, vs_cr, va_cr, upper_bound_m_l):

        self.v_sound_cruise=v_sound_cruise # -, speed of sound at cruise
        self.v_sound_dive = v_sound_dive  # -, speed of sound at dive
        self.vs_to= vs_to # m/s, stall speed take off
        self.va_to = va_to # m/s, va take off speed
        self.vs_cr = vs_cr # m/s, stall speed cruise
        self.va_cr = va_cr # m/s, va cruise speed
        self.upper_bound_m_l = upper_bound_m_l # -

    def outcome_dive(self):

        if 0.8*self.v_sound_dive>self.v_sound_cruise:
            return True
        else:
            return False

    def stall_speed(self):

        if self.va_cr>self.vs_cr*math.sqrt(self.upper_bound_m_l) and self.va_to>self.vs_to*math.sqrt(self.upper_bound_m_l):

            return True

        else:

            return False

class TakeOffProcedure():

    """
    This class includes all the data and the calculations needed to find the required runway length
    """

    def __init__(self, Ab,front_drag_area,v_tip,A, MTOW_weight, density_sea_level, drag_coefficient_blades, induced_drag):

        self.Ab = Ab # m^2, Area of the blades of one rotor
        self.drag_coefficient_blades = drag_coefficient_blades # -, coefficient of drag of blades
        self.front_drag_area = front_drag_area # m^2, equivalent front area
        self.v_tip = v_tip # m/s, tip speed velocity of the rotor
        self.mu =0.15 #
        self.ct_sigma = 0.04 # ct/sigma, other values necessary to retrieve the value of 'e' are 0.08 and 0.12
        self.induced_drag = induced_drag # N
        self.alpha_tpp= -57.3*self.induced_drag/self.MTOW_weight# degrees, necessary to find the value of e
        self.A = A # m^2, Disk area of rotor
        self.MTOW_weight = MTOW_weight # kg
        self.density_sea_level = density_sea_level #kg/m^3
        self.lift_coefficient_sigma = 20 # other values as shown in page 356 of the green book are 15, 10, 5, 0

    def vmin(self,e):
        """
        This function allows to find velocity at minimum power for an helicopter (applied for a rotor) once 'e'(Induced Efficiency Factor) is found (page 128 green book)
        :return: v minimum speed
        """
        vmin=Symbol('vmin')
        vmin_solution=solve(vmin^4+0.5*(self.Ab)*(self.v_tip)*(self.drag_coefficient_blades)*(vmin^3)/(self.front_drag_area) -((self.MTOW_weight/self.density_sea_level)^2)/(3*e*self.front_drag_area*self.A),vmin)
        return vmin_solution

    def vcrit(self, v_critical):
        """
        Given the value of lift_coefficient_sigma, read at page 365 of the green book what is v critical and report it in m/s
        :param v_critical:
        :return: v_critical
        """
        return v_critical


class AicraftLandingProcedure():
    """
    This class contains all the necessary calculations and data needed to find the stall speed which is the minimum speed at landing.
    The calculations are based on an hypothetical horizontal landing given an airborne landing distance. No consideration is made on the ground run distance needed to stop the aircraft
    """
    def __init__(self):

        self.delta_n=0.15 # load factor of 1.15 at landing for propeller aircraft
        self.gamma= 3 # degree landing angle
        self.h_scr= 15.24 # m, landign screen height distance as based on ADSEE 1 slide 45 and NLR investigation (50ft)
        self.h_flare = 6.1 # m , based on NLR landing investigation for jet aircraft but I think it is ok (20ft-30ft)
        self.g = 9.80665 # m/s^2
        self.airborne_distance =410 #m, based on NLR landing investigation for jet aircraft but I think it is ok (400m-430m)
        self.runway_length=1508 #m, london city airport runway length, shortest among the selected airports

    def vapproach(self):
        """
        Calculate the approach speed based on statistical flare height and on airborne horizontal statistical length.
        Also the statistical estimation provided in ADSEE I is considered even if not accurate enough since it is based on runway length while we
        want our emergency landing condition to be dependent on the approach phase
        :return: v_approach
        """
        R=self.h_flare/(1-math.cos(math.radians(self.gamma)))
        v_approach=math.sqrt(R*(self.delta_n)*(self.g))

        R1=(self.airborne_distance-self.h_scr/math.tan(math.radians(self.gamma)))/(math.sin(math.radians(self.gamma))-1+math.cos(math.radians(self.gamma)))
        v_approach_1 = math.sqrt(R1 * (self.delta_n) * (self.g))

        v_approach_2=math.sqrt(self.runway_length/0.5847) #empirical relationship for cs-25 jet ADSEE I


        return v_approach,v_approach_1, v_approach_2

    def vmin(self,v_approach):
        """
        Calculate minimum speed

        :param v_approach:
        :return:v_min
        """
        v_min=v_approach/1.3
        return v_min

def maneuver_diagram(HorizontalSpeeds, Rotorcraft, LimitLoadsCS29, LiftCoefficients, plot=True):

    # find speeds for to based on loading

    vs_to=math.sqrt((2*Rotorcraft.MTOW_weight*LimitLoadsCS29.unity_load)/(Rotorcraft.density_sea_level*LiftCoefficients.cl_stall_sealevel*Rotorcraft.wing_area))
    va_to=math.sqrt((2*Rotorcraft.MTOW_weight*LimitLoadsCS29.upper_bound_m_l)/(Rotorcraft.density_sea_level*LiftCoefficients.cl_stall_sealevel*Rotorcraft.wing_area))

    # vc_to=math.sqrt((2*Rotorcraft.MTOW_weight*LimitLoadsCS29.lower_bound_m_l)/(Rotorcraft.density_sea_level*(-LiftCoefficients.cl_cruise_sealevel)*Rotorcraft.wing_area))
    # vd_to=math.sqrt((2*Rotorcraft.MTOW_weight*LimitLoadsCS29.upper_bound_m_l)/(Rotorcraft.density_sea_level*LiftCoefficients.cl_cruise_sealevel*Rotorcraft.wing_area))

    H_to= math.sqrt((2*Rotorcraft.MTOW_weight*LimitLoadsCS29.lower_bound_m_l)/(Rotorcraft.density_sea_level*(-LiftCoefficients.cl_stall_sealevel)*Rotorcraft.wing_area))

    vs_full_flaps=math.sqrt((2*Rotorcraft.MTOW_weight*LimitLoadsCS29.landing)/(Rotorcraft.density_sea_level*LiftCoefficients.cl_stall_full_flaps*Rotorcraft.wing_area))

    HorizontalSpeeds.vs_full_flaps = vs_full_flaps

    HorizontalSpeeds.vs_to =vs_to
    HorizontalSpeeds.va_to =va_to
    HorizontalSpeeds.H_to =H_to


    # find speeds for cruise based on loading

    vs_cr=math.sqrt((2*Rotorcraft.MTOW_weight*LimitLoadsCS29.unity_load)/(Rotorcraft.density_cruise*LiftCoefficients.cl_stall_sealevel*Rotorcraft.wing_area))
    va_cr=math.sqrt((2*Rotorcraft.MTOW_weight*LimitLoadsCS29.upper_bound_m_l)/(Rotorcraft.density_cruise*LiftCoefficients.cl_stall_sealevel*Rotorcraft.wing_area))

    # vc_cr=math.sqrt((2*Rotorcraft.MTOW_weight*LimitLoadsCS29.lower_bound_m_l)/(Rotorcraft.density_cruise*(-LiftCoefficients.cl_cruise_sealevel)*Rotorcraft.wing_area))
    # vd_cr=math.sqrt((2*Rotorcraft.MTOW_weight*LimitLoadsCS29.upper_bound_m_l)/(Rotorcraft.density_cruise*LiftCoefficients.cl_cruise_sealevel*Rotorcraft.wing_area))

    H_cr=math.sqrt((2*Rotorcraft.MTOW_weight*LimitLoadsCS29.lower_bound_m_l)/(Rotorcraft.density_cruise*(-LiftCoefficients.cl_stall_sealevel)*Rotorcraft.wing_area))


    HorizontalSpeeds.vs_cr =vs_cr
    HorizontalSpeeds.va_cr =va_cr
    HorizontalSpeeds.H_cr =H_cr

    # make arrays sea levels

    if plot==True:

        x=HorizontalSpeeds.to_speeds()
        y=LimitLoadsCS29.load_factors()

        # x_inter_quadratic = x[0:3]
        # y_inter_quadratic = y[0:3]
        # f_quadratic = interp1d(x_inter_quadratic, y_inter_quadratic, kind='quadratic')
        # x_new_quadratic = np.linspace(x[0], x[2], 100)

        # x1_inter_quadratic = x[6:]
        # y1_inter_quadratic = y[6:]
        # # f_quadratic = interp1d(x1_inter_quadratic, y1_inter_quadratic, kind='quadratic')
        # x1_new_quadratic = np.linspace(x[6], x[-1], 100)
        #
        # f_quadratic_1=np.poly1d(np.polyfit(x1_inter_quadratic, y1_inter_quadratic, 2))

        velocities_a=np.linspace(0, HorizontalSpeeds.va_to, 1000)
        velocities_a_minus=np.linspace(0, HorizontalSpeeds.vs_to, 1000)

        n_velocities_a=0.5*Rotorcraft.density_sea_level*LiftCoefficients.cl_stall_sealevel*Rotorcraft.wing_area*(velocities_a**2)/Rotorcraft.MTOW_weight

        n_velocities_a_minus = -0.5 * Rotorcraft.density_sea_level * LiftCoefficients.cl_stall_sealevel * Rotorcraft.wing_area * (
                velocities_a_minus ** 2) / Rotorcraft.MTOW_weight

        # plot maneuver loading sea level

        plt.title('Maneuver loading for sealevel condition')
        plt.scatter(x,y, color='b')

        plt.plot(velocities_a_minus,n_velocities_a_minus, 'b')
        plt.plot(velocities_a,n_velocities_a, 'b')
        plt.plot(x[3:8],y[3:8], 'b')

        plt.xlabel(r'True Airspeed [$\frac{m}{s}$]', fontsize=12)
        plt.ylabel(r'Load Factor [-]', fontsize=12)

        plt.show()

        # make arrays cruise

        x=HorizontalSpeeds.cr_speeds()
        y=LimitLoadsCS29.load_factors()

        # x_inter_quadratic = x[0:3]
        # y_inter_quadratic = y[0:3]
        # f_quadratic = interp1d(x_inter_quadratic, y_inter_quadratic, kind='quadratic')
        # x_new_quadratic = np.linspace(x[0], x[2], 100)
        #
        # x1_inter_quadratic = x[6:]
        # y1_inter_quadratic = y[6:]
        # # f_quadratic = interp1d(x1_inter_quadratic, y1_inter_quadratic, kind='quadratic')
        # x1_new_quadratic = np.linspace(x[6], x[-1], 100)
        #
        # f_quadratic_1=np.poly1d(np.polyfit(x1_inter_quadratic, y1_inter_quadratic, 2))

        velocities_a=np.linspace(0, HorizontalSpeeds.va_cr, 1000)
        velocities_a_minus=np.linspace(0, HorizontalSpeeds.vs_cr, 1000)

        n_velocities_a=0.5*Rotorcraft.density_cruise*LiftCoefficients.cl_stall_cruise*Rotorcraft.wing_area*(velocities_a**2)/Rotorcraft.MTOW_weight

        n_velocities_a_minus = -0.5 * Rotorcraft.density_cruise * LiftCoefficients.cl_stall_cruise * Rotorcraft.wing_area * (
                velocities_a_minus ** 2) / Rotorcraft.MTOW_weight


        # plot maneuver loading cruise

        plt.title('Maneuver loading for cruise condition')
        plt.scatter(x,y, color='b')

        plt.plot(velocities_a_minus,n_velocities_a_minus, 'b')
        plt.plot(velocities_a,n_velocities_a, 'b')
        plt.plot(x[3:8],y[3:8], 'b')

        plt.xlabel(r'True Airspeed [$\frac{m}{s}$]', fontsize=12)
        plt.ylabel(r'Load Factor [-]', fontsize=12)

        plt.show()


def gust_diagram(Rotorcraft,GustSpeedCS25, plot=True):

    #find delta n sea level

    delta_n_ub_to=1+GustSpeedCS25.ub*GustSpeedCS25.vb_to*Rotorcraft.lift_curve*Rotorcraft.density_sea_level*GustSpeedCS25.Kg_sea_level_gust/(2*Rotorcraft.wing_loading)
    delta_n_uc_to =1+ GustSpeedCS25.uc*GustSpeedCS25.vc_to*Rotorcraft.lift_curve*Rotorcraft.density_sea_level*GustSpeedCS25.Kg_sea_level_gust/(2*Rotorcraft.wing_loading)
    delta_n_ud_to = 1+GustSpeedCS25.ud*GustSpeedCS25.vd_to*Rotorcraft.lift_curve*Rotorcraft.density_sea_level*GustSpeedCS25.Kg_sea_level_gust/(2*Rotorcraft.wing_loading)

    # find delta n cruise

    delta_n_ub_cr=1+GustSpeedCS25.ub*GustSpeedCS25.vb_cr*Rotorcraft.lift_curve*Rotorcraft.density_sea_level*GustSpeedCS25.Kg_cruise_gust/(2*Rotorcraft.wing_loading)
    delta_n_uc_cr = 1+GustSpeedCS25.uc*GustSpeedCS25.vc_cr*Rotorcraft.lift_curve*Rotorcraft.density_sea_level*GustSpeedCS25.Kg_cruise_gust/(2*Rotorcraft.wing_loading)
    delta_n_ud_cr = 1+GustSpeedCS25.ud*GustSpeedCS25.vd_cr*Rotorcraft.lift_curve*Rotorcraft.density_sea_level*GustSpeedCS25.Kg_cruise_gust/(2*Rotorcraft.wing_loading)

    delta_n_u=[delta_n_ub_to,delta_n_uc_to,delta_n_ud_to,delta_n_ub_cr,delta_n_uc_cr,delta_n_ud_cr]



    #plot sea level gust loading

    if plot==True:

        plt.title('Gust loading for sealevel condition')
        plt.plot([0,GustSpeedCS25.vb_to],[1,delta_n_ub_to], 'b')
        plt.plot([0,GustSpeedCS25.vc_to], [1, delta_n_uc_to], 'b')
        plt.plot([0, GustSpeedCS25.vd_to], [1, delta_n_ud_to], 'b')
        plt.plot([0,GustSpeedCS25.vb_to],[1,1-delta_n_ub_to], 'b')
        plt.plot([0,GustSpeedCS25.vc_to], [1, 1-delta_n_uc_to], 'b')
        plt.plot([0, GustSpeedCS25.vd_to], [1, 1-delta_n_ud_to], 'b')
        plt.xlabel(r'True Airspeed [$\frac{m}{s}$]', fontsize=12)
        plt.ylabel(r'Load Factor [-]', fontsize=12)
        plt.show()


        #plot cruise gust loading

        plt.title('Gust loading for cruise condition')
        plt.plot([0,GustSpeedCS25.vb_cr],[1,delta_n_ub_cr], 'b')
        plt.plot([0,GustSpeedCS25.vc_cr], [1, delta_n_uc_cr], 'b')
        plt.plot([0, GustSpeedCS25.vd_cr], [1, delta_n_ud_cr], 'b')
        plt.plot([0,GustSpeedCS25.vb_cr],[1,1-delta_n_ub_cr], 'b')
        plt.plot([0,GustSpeedCS25.vc_cr], [1, 1-delta_n_uc_cr], 'b')
        plt.plot([0, GustSpeedCS25.vd_cr], [1, 1-delta_n_ud_cr], 'b')
        plt.xlabel(r'True Airspeed [$\frac{m}{s}$]', fontsize=12)
        plt.ylabel(r'Load Factor [-]', fontsize=12)
        plt.show()

    return delta_n_u

# AicraftLandingProcedure=AicraftLandingProcedure()
# print ('V landing/approach based on flare height',AicraftLandingProcedure.vapproach()[0])
# print ('V landing/approach based on airborne distance',AicraftLandingProcedure.vapproach()[1])
# print ('V landing/approach based on runway length',AicraftLandingProcedure.vapproach()[2])
#
# print ('V min based on flare height',AicraftLandingProcedure.vmin(AicraftLandingProcedure.vapproach()[0]))
# print ('V min based on airborne distance',AicraftLandingProcedure.vmin(AicraftLandingProcedure.vapproach()[1]))
# print ('V min based on runway length',AicraftLandingProcedure.vmin(AicraftLandingProcedure.vapproach()[2]))

def combined_load_diagram(HorizontalSpeeds, Rotorcraft, LimitLoadsCS29, LiftCoefficients, delta_n_u, GustSpeedCS25, plot=True, plot_final=False):

    """
    Take all the parameters needed to plot the combined loading diagram and return the maximum load
    :param HorizontalSpeeds:
    :param Rotorcraft:
    :param LimitLoadsCS29:
    :param LiftCoefficients:
    :param delta_n_u:
    :param GustSpeedCS25:
    :param plot: False or True according if plotting is desired to understand the limiting condition.
    This is essential to validate whether or not the maximum load returned is actually the limiting case.
    :return: maximum_load
    """

    # find speeds for to based on loading

    vs_to = math.sqrt((2 * Rotorcraft.MTOW_weight * LimitLoadsCS29.unity_load) / (
                Rotorcraft.density_sea_level * LiftCoefficients.cl_stall_sealevel * Rotorcraft.wing_area))
    va_to = math.sqrt((2 * Rotorcraft.MTOW_weight * LimitLoadsCS29.upper_bound_m_l) / (
                Rotorcraft.density_sea_level * LiftCoefficients.cl_stall_sealevel * Rotorcraft.wing_area))

    # vc_to = math.sqrt((2 * Rotorcraft.MTOW_weight * LimitLoadsCS29.lower_bound_m_l) / (
    #             Rotorcraft.density_sea_level * (-LiftCoefficients.cl_cruise_sealevel) * Rotorcraft.wing_area))
    # vd_to = math.sqrt((2 * Rotorcraft.MTOW_weight * LimitLoadsCS29.upper_bound_m_l) / (
    #             Rotorcraft.density_sea_level * LiftCoefficients.cl_cruise_sealevel * Rotorcraft.wing_area))

    H_to = math.sqrt((2 * Rotorcraft.MTOW_weight * LimitLoadsCS29.lower_bound_m_l) / (
                Rotorcraft.density_sea_level * (-LiftCoefficients.cl_stall_sealevel) * Rotorcraft.wing_area))

    vs_full_flaps=math.sqrt((2*Rotorcraft.MTOW_weight*LimitLoadsCS29.landing)/(Rotorcraft.density_sea_level*LiftCoefficients.cl_stall_full_flaps*Rotorcraft.wing_area))

    HorizontalSpeeds.vs_full_flaps = vs_full_flaps

    HorizontalSpeeds.vs_to = vs_to
    HorizontalSpeeds.va_to = va_to
    HorizontalSpeeds.H_to = H_to

    # find speeds for cruise based on loading

    vs_cr = math.sqrt((2 * Rotorcraft.MTOW_weight * LimitLoadsCS29.unity_load) / (
                Rotorcraft.density_cruise * LiftCoefficients.cl_stall_sealevel * Rotorcraft.wing_area))
    va_cr = math.sqrt((2 * Rotorcraft.MTOW_weight * LimitLoadsCS29.upper_bound_m_l) / (
                Rotorcraft.density_cruise * LiftCoefficients.cl_stall_sealevel * Rotorcraft.wing_area))

    # vc_cr = math.sqrt((2 * Rotorcraft.MTOW_weight * LimitLoadsCS29.lower_bound_m_l) / (
    #             Rotorcraft.density_cruise * (-LiftCoefficients.cl_cruise_sealevel) * Rotorcraft.wing_area))
    # vd_cr = math.sqrt((2 * Rotorcraft.MTOW_weight * LimitLoadsCS29.upper_bound_m_l) / (
    #             Rotorcraft.density_cruise * LiftCoefficients.cl_cruise_sealevel * Rotorcraft.wing_area))

    H_cr = math.sqrt((2 * Rotorcraft.MTOW_weight * LimitLoadsCS29.lower_bound_m_l) / (
                Rotorcraft.density_cruise * (-LiftCoefficients.cl_stall_sealevel) * Rotorcraft.wing_area))

    HorizontalSpeeds.vs_cr = vs_cr
    HorizontalSpeeds.va_cr = va_cr
    HorizontalSpeeds.H_cr = H_cr

    # make arrays sea levels

    x_to = HorizontalSpeeds.to_speeds()
    y_to = LimitLoadsCS29.load_factors()

    velocities_a_to = np.linspace(0, HorizontalSpeeds.va_to, 1000)
    velocities_full_flaps_to=np.linspace(0, HorizontalSpeeds.vs_full_flaps, 1000)
    velocities_a_minus_to = np.linspace(0, HorizontalSpeeds.vs_to, 1000)

    n_velocities_full_flaps_to = 0.5 * Rotorcraft.density_sea_level * LiftCoefficients.cl_stall_full_flaps * Rotorcraft.wing_area * (
            velocities_full_flaps_to ** 2) / Rotorcraft.MTOW_weight

    n_velocities_a_to = 0.5 * Rotorcraft.density_sea_level * LiftCoefficients.cl_stall_sealevel * Rotorcraft.wing_area * (
            velocities_a_to ** 2) / Rotorcraft.MTOW_weight

    n_velocities_a_minus_to = -0.5 * Rotorcraft.density_sea_level * LiftCoefficients.cl_stall_sealevel * Rotorcraft.wing_area * (
            velocities_a_minus_to ** 2) / Rotorcraft.MTOW_weight


    vb_to = Symbol('vb_to')
    vb_to_solution = solve(0.5 * Rotorcraft.density_sea_level * LiftCoefficients.cl_stall_sealevel * Rotorcraft.wing_area * (
            vb_to ** 2) / Rotorcraft.MTOW_weight-1-GustSpeedCS25.ub*vb_to*Rotorcraft.lift_curve*Rotorcraft.density_sea_level*GustSpeedCS25.Kg_sea_level_gust/(2*Rotorcraft.wing_loading), vb_to)

    vb_to_solution=min(max(vb_to_solution),GustSpeedCS25.vc_to)

    n_intersection_to=1+GustSpeedCS25.ub*vb_to_solution*Rotorcraft.lift_curve*Rotorcraft.density_sea_level*GustSpeedCS25.Kg_sea_level_gust/(2*Rotorcraft.wing_loading)

    n_intersection_to_final=min(n_intersection_to,3.8)
    vb_to_solution_final=(n_intersection_to_final-1)/(GustSpeedCS25.ub*Rotorcraft.lift_curve*Rotorcraft.density_sea_level*GustSpeedCS25.Kg_sea_level_gust/(2*Rotorcraft.wing_loading))


    if plot==True:
        # plot maneuver loading sea level

        fig, ax = plt.subplots()

        plt.title('Combined loading diagram for sealevel condition')
        plt.scatter(x_to, y_to, color='b')

        plt.plot(velocities_full_flaps_to, n_velocities_full_flaps_to, 'b')
        plt.plot(velocities_a_minus_to, n_velocities_a_minus_to, 'b')
        plt.plot(velocities_a_to, n_velocities_a_to, 'b')
        plt.plot(x_to[3:8], y_to[3:8], 'b')

        plt.plot([0, GustSpeedCS25.vb_to], [1, delta_n_u[0]], 'r--')
        plt.plot([0, GustSpeedCS25.vc_to], [1, delta_n_u[1]], 'r--')
        plt.plot([0, GustSpeedCS25.vd_to], [1, delta_n_u[2]], 'r--')
        plt.plot([0, GustSpeedCS25.vb_to], [1, 1 - delta_n_u[0]], 'r--')
        plt.plot([0, GustSpeedCS25.vc_to], [1, 1 - delta_n_u[1]], 'r--')
        plt.plot([0, GustSpeedCS25.vd_to], [1, 1 - delta_n_u[2]], 'r--')

        ax.axhline(y=0, color='black')
        plt.xlabel(r'True Airspeed [$\frac{m}{s}$]', fontsize=12)
        plt.ylabel(r'Load Factor [-]', fontsize=12)

        plt.scatter(vb_to_solution, n_intersection_to, color='black')
        plt.scatter(vb_to_solution_final, n_intersection_to_final, color='red')
        plt.show()

    # make arrays cruise

    x_cr = HorizontalSpeeds.cr_speeds()
    y_cr = LimitLoadsCS29.load_factors()

    velocities_a_cr = np.linspace(0, HorizontalSpeeds.va_cr, 1000)
    velocities_a_minus_cr = np.linspace(0, HorizontalSpeeds.vs_cr, 1000)

    n_velocities_a_cr = 0.5 * Rotorcraft.density_cruise * LiftCoefficients.cl_stall_cruise * Rotorcraft.wing_area * (
                velocities_a_cr ** 2) / Rotorcraft.MTOW_weight

    n_velocities_a_minus_cr = -0.5 * Rotorcraft.density_cruise * LiftCoefficients.cl_stall_cruise * Rotorcraft.wing_area * (
                velocities_a_minus_cr ** 2) / Rotorcraft.MTOW_weight

    # finding maximum load and intersection

    vb_cr = Symbol('vb_cr')
    vb_cr_solution = solve(0.5 * Rotorcraft.density_cruise * LiftCoefficients.cl_stall_cruise * Rotorcraft.wing_area * (
            vb_cr ** 2) / Rotorcraft.MTOW_weight-1-GustSpeedCS25.ub*vb_cr*Rotorcraft.lift_curve*Rotorcraft.density_cruise*GustSpeedCS25.Kg_cruise_gust/(2*Rotorcraft.wing_loading), vb_cr)

    vb_cr_solution = min(max(vb_cr_solution), GustSpeedCS25.vc_cr)

    n_intersection_cr=1+GustSpeedCS25.ub*vb_cr_solution*Rotorcraft.lift_curve*Rotorcraft.density_cruise*GustSpeedCS25.Kg_cruise_gust/(2*Rotorcraft.wing_loading)

    n_intersection_cr_final=min(n_intersection_cr,3.8)
    vb_cr_solution_final=(n_intersection_cr_final-1)/(GustSpeedCS25.ub*Rotorcraft.lift_curve*Rotorcraft.density_cruise*GustSpeedCS25.Kg_cruise_gust/(2*Rotorcraft.wing_loading))

    if plot == True:
        # plot maneuver loading cruise

        fig, ax = plt.subplots()

        plt.title('Combined loading diagram for cruise condition')
        plt.scatter(x_cr, y_cr, color='b')

        plt.plot(velocities_a_minus_cr, n_velocities_a_minus_cr, 'b')
        plt.plot(velocities_a_cr, n_velocities_a_cr, 'b')
        plt.plot(x_cr[3:8], y_cr[3:8], 'b')

        plt.plot([0, GustSpeedCS25.vb_cr], [1, delta_n_u[3]], 'r--')
        plt.plot([0, GustSpeedCS25.vc_cr], [1, delta_n_u[4]], 'r--')
        plt.plot([0, GustSpeedCS25.vd_cr], [1, delta_n_u[5]], 'r--')
        plt.plot([0, GustSpeedCS25.vb_cr], [1, 1 - delta_n_u[3]], 'r--')
        plt.plot([0, GustSpeedCS25.vc_cr], [1, 1 - delta_n_u[4]], 'r--')
        plt.plot([0, GustSpeedCS25.vd_cr], [1, 1 - delta_n_u[5]], 'r--')

        ax.axhline(y=0, color='black')
        plt.xlabel(r'True Airspeed [$\frac{m}{s}$]', fontsize=12)
        plt.ylabel(r'Load Factor [-]', fontsize=12)
        plt.scatter(vb_cr_solution, n_intersection_cr, color='black')
        plt.scatter(vb_cr_solution_final, n_intersection_cr_final, color='red')
        plt.show()

    maximum_load = max(n_intersection_to_final, 3.5, n_intersection_cr_final)

    if plot_final:

        # plot the limiting case of sea level condition

        fig, ax = plt.subplots()

        plt.title('Combined loading diagram for sealevel condition')
        plt.scatter(x_to, y_to, color='b')

        plt.plot(velocities_full_flaps_to, n_velocities_full_flaps_to, 'b')
        plt.plot(velocities_a_minus_to, n_velocities_a_minus_to, 'b')
        plt.plot(velocities_a_to, n_velocities_a_to, 'b')
        plt.plot(x_to[3:8], y_to[3:8], 'b')

        plt.plot([0, GustSpeedCS25.vb_to], [1, delta_n_u[0]], 'r--')
        plt.plot([0, GustSpeedCS25.vc_to], [1, delta_n_u[1]], 'r--')
        plt.plot([0, GustSpeedCS25.vd_to], [1, delta_n_u[2]], 'r--')
        plt.plot([0, GustSpeedCS25.vb_to], [1, 1 - delta_n_u[0]], 'r--')
        plt.plot([0, GustSpeedCS25.vc_to], [1, 1 - delta_n_u[1]], 'r--')
        plt.plot([0, GustSpeedCS25.vd_to], [1, 1 - delta_n_u[2]], 'r--')

        # final intersection


        v_1_solution=19.45 # magnitude

        n_1=0.5*Rotorcraft.density_sea_level*LiftCoefficients.cl_stall_sealevel*Rotorcraft.wing_area*(v_1_solution**2)/Rotorcraft.MTOW_weight

        v_2=vb_to_solution_final
        n_2=n_intersection_to_final

        v_3=Symbol('v_3')
        v_3_solution= solve(LimitLoadsCS29.upper_bound_m_l-(1+ GustSpeedCS25.uc*v_3*Rotorcraft.lift_curve*Rotorcraft.density_sea_level*GustSpeedCS25.Kg_sea_level_gust/(2*Rotorcraft.wing_loading)),v_3)

        v_3_solution=v_3_solution[0]

        n_3=1 + GustSpeedCS25.uc * v_3_solution * Rotorcraft.lift_curve * Rotorcraft.density_sea_level * GustSpeedCS25.Kg_sea_level_gust / (2 * Rotorcraft.wing_loading)

        v_4=HorizontalSpeeds.vd_to
        n_4=delta_n_u[2]

        v_5=HorizontalSpeeds.vd_to
        n_5=1-delta_n_u[2]

        v_6=HorizontalSpeeds.vd_to*0.99
        n_6=-0.4*n_intersection_to_final

        n_7=-0.4*n_intersection_to_final
        v_7=59.46 #m/s

        plt.scatter([v_1_solution,v_2,v_3_solution,v_4,v_5,v_6,v_7,v_1_solution],[n_1,n_2,n_3,n_4,n_5,n_6,n_7,n_1], color='black')
        plt.plot([v_1_solution,v_2,v_3_solution,v_4,v_5,v_6,v_7,v_1_solution],[n_1,n_2,n_3,n_4,n_5,n_6,n_7,n_1], color='black')

        ax.axhline(y=0, color='black')
        plt.xlabel(r'True Airspeed [$\frac{m}{s}$]', fontsize=12)
        plt.ylabel(r'Load Factor [-]', fontsize=12)
        plt.show()


    return maximum_load, vb_to_solution, vb_cr_solution

def main():
   pass

if __name__ == "__main__":
    LiftCoefficients=LiftCoefficients(1.386, 0.436, 0.436,1.386, 0.436,0.436,2) # assumption that cl dive and cruise are the same
    HorizontalSpeeds=HorizontalSpeeds(97.22,97.22)
    LimitLoadsCS29=LimitLoadsCS29()
    Rotorcraft=Rotorcraft(4000,2, 1.2, 2000, 0.12066, 3, 0.3, 21.48,2,97.22)
    maneuver_diagram(HorizontalSpeeds, Rotorcraft, LimitLoadsCS29, LiftCoefficients, plot=False)
    GustSpeedCS25=GustSpeedCS25(Rotorcraft,97.22,97.22, HorizontalSpeeds.vs_to,HorizontalSpeeds.vs_cr,Rotorcraft.lift_curve,Rotorcraft.wing_loading,Rotorcraft.chord,Rotorcraft.density_sea_level, Rotorcraft.density_cruise, 1.1) # speeds at sea level and cruise are the same
    print(GustSpeedCS25.vb_to, GustSpeedCS25.vc_to, GustSpeedCS25.vd_to)
    delta_n_u=gust_diagram(Rotorcraft,GustSpeedCS25, plot=False)
    print (combined_load_diagram(HorizontalSpeeds, Rotorcraft, LimitLoadsCS29, LiftCoefficients, delta_n_u, GustSpeedCS25, plot=False, plot_final=True))

    print (GustSpeedCS25.vb_to, GustSpeedCS25.vc_to, GustSpeedCS25.vd_to)
