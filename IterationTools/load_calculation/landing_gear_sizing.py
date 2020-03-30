import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import math


"""This code is used to size the landing gear according to CS-25 requirements since emergency landing will occur in aircraft configuration"""

class LandingGear():

    """ The initialization of the function takes the parameters of the landing gear
    such as the stroke length, the absorbtion coefficient of the stroke, the tire dimensions,the absorption coefficient of the tire"""

    def __init__(self,MTOW,landing_wing_loading, stroke_length,nt_tire_main,ns,st_main,nt_tire_nose,st_nose, number_strut_mlg,ln,lm,D0_main,d_main,b_main,loaded_radius_main,D0_nose,d_nose,b_nose,loaded_radius_nose):
        """

        :param MTOW:
        :param landing_wing_loading:
        :param stroke_length:
        :param nt_tire_main: tire absorbtion efficiency
        :param ns: shock absorbtion efficiency
        :param st_main: allowable tire deflection main landing gear
        :param nt_tire_nose: tire absorbtion efficiency
        :param st_nose: allowable tire deflection nose landing gear
        :param number_strut_mlg:
        :param ln: horizontal distance from the main landing gear to the cg
        :param lm: horizontal distance from the nose landing gear to the cg
        :param D0_main: big diameter of the main landing gear tire in inch
        :param d_main: small diameter of the main landing gear tire in inch
        :param b_main: width of the main landing gear tire in inch
        :param loaded_radius_main:
        :param D0_nose: big diameter of the nose landing gear tire in inch
        :param d_nose: small diameter of the nose landing gear tire in inch
        :param b_nose: width of the nose landing gear tire in inch
        :param loaded_radius_nose:
        """

        # relevant constant input and conversions
        self.MTOW= MTOW#kg
        self.g=9.80665 # m/s^s gravity constant
        self.kg_to_lbs=2.20462 # conversion from kg to lbs
        self.N_to_lbs=2.20462/9.80665 # conversion from N to lbs
        self.m_to_ft=1./0.3048 # conversion from m to ft
        self.inch_to_m=0.025
        self.ln=ln # distance from nose landing gear to cg
        self.lm=lm # distance from main landing gear to cg

        self.D0_main=D0_main*self.inch_to_m  #m
        self.d_main=d_main*self.inch_to_m #m
        self.b_main=b_main*self.inch_to_m #m
        self.loaded_radius_main=loaded_radius_main*self.inch_to_m #m
        self.D0_nose=D0_nose*self.inch_to_m #m
        self.d_nose=d_nose*self.inch_to_m #m
        self.b_nose=b_nose*self.inch_to_m #m
        self.loaded_radius_nose=loaded_radius_nose*self.inch_to_m #m

        #calculation of touch down speed

        self.landing_weight_proportion=0.88 # proportion of the total mass to use for landing gear load calculations (Roskam I)
        self.Wl=self.MTOW*self.landing_weight_proportion
        self.touch_down_speed=min(4.4*(landing_wing_loading*self.N_to_lbs/(self.m_to_ft)**2),10) #10 ft/s is the maximum touchdown speed allowable for CS-23


        # this part is not directly relevant for sizing
        self.LCN_m=0.0001847 # from regression
        self.LCN_i=20.80465 # from regression
        self.LCN=min(20,self.LCN_m*self.MTOW+self.LCN_i)
        self.inflation_pressure=680*math.log(self.LCN)-680 # kPa


        # landing gear tire and strut design
        self.load_proportion_ml=0.92
        self.load_proportion_nl=0.08
        self.static_nose_lg_load=self.load_proportion_nl*self.MTOW # kg
        self.static_main_lg_load=self.load_proportion_ml*self.MTOW/2 # kg

        self.Ng=3 # number of g which determine the dynamic loading
        self.dynamic_main_lg_load=self.Ng*self.static_main_lg_load

        #maximum energy to be absorbed based on touch down speed

        self.energy_touch_down=0.5*self.MTOW*((self.touch_down_speed*1/self.m_to_ft)**2)/self.g # J
        self.required_stroke_length_main=(0.5*self.Wl*((self.touch_down_speed*1/self.m_to_ft)**2)/(number_strut_mlg*self.static_main_lg_load*self.Ng*self.g)-nt_tire_main*st_main)/ns

        self.required_stroke_length_nose = (0.5 * self.static_nose_lg_load * ((self.touch_down_speed * 1 / self.m_to_ft) ** 2) / (ns * self.static_main_lg_load * self.Ng * self.g) - nt_tire_nose * st_nose) / ns

        self.maximum_energy_sustainable=number_strut_mlg*self.Ng*self.static_main_lg_load*(nt_tire_main*st_main+ns*self.required_stroke_length_main)

        self.maximum_energy_sustainable_design = number_strut_mlg * self.Ng * self.static_main_lg_load * (
                    nt_tire_main * st_main + ns * stroke_length)

        # landing gear dimensions

        self.diameter_landing_gear_main=(0.041+0.0025*((self.static_main_lg_load*self.kg_to_lbs)**0.5))*1/self.m_to_ft #m
        self.diameter_landing_gear_nose = (0.041 + 0.0025 * ((self.static_nose_lg_load * self.kg_to_lbs) ** 0.5))*1/self.m_to_ft  # m , Not accurate estimation



    def lateral_tip_over(self,tip_back_landing_angle,fuselage_width, span_engine,phi, plot=True):
        """
        Enter the psi angle and the height of the main landing gear and find the half track of the main landing gear.
        Enter psi in degrees
        :param: tip_back_landing_angle: tip back angle in degrees for landing purposes, standard approach angle is a bit less than the stall angle which from the lift curve slope, for an approach speed of 50.48m/s and a cl of 1.15 is 8.65 degrees
        :param: fuselage width in meters
        :param: span engines is the distance between the propellers of the wing
        :param phi is the ground engine clearance which is typically 5 degrees
        :return
        """
        tip_back_landing_angle=tip_back_landing_angle*math.pi/180
        lms=np.linspace(0,0.5,100)
        psi_constant = 55 * math.pi / 180
        heights_landing_gear = lms / np.tan(tip_back_landing_angle) - self.required_stroke_length_main
        y_mlgs = (lms + self.ln) / (np.sqrt(((self.ln) ** 2) * (math.tan(psi_constant) ** 2) / (heights_landing_gear ** 2) - 1))


        # height_landing_gear = self.lm / math.tan(tip_back_landing_angle) - self.required_stroke_length_main
        # y_mlg=np.linspace(0,fuselage_width/2,100)
        #
        # delta=np.arctan(y_mlg*2/(2*(self.ln+self.lm)))
        # psi=np.arctan(height_landing_gear/(self.ln*np.sin(delta)))
        # y_mlg_final_index=np.where(y_mlg>(self.lm+self.ln)/(np.sqrt(((self.ln)**2)*(np.tan(psi)**2)/(height_landing_gear**2)-1)))
        #
        # y_mlg_final=y_mlg[y_mlg_final_index]
        # psi_final=psi[y_mlg_final_index]
        # delta_final = delta[y_mlg_final_index]
        # index_optimized=np.argmin(psi_final)
        #
        # psi_optimized=psi_final[index_optimized]
        # delta_optimized = delta_final[index_optimized]
        # y_mlg_final_optimized = y_mlg_final[index_optimized]

        tip_engine_clearance=(y_mlgs-span_engine/2)*(-math.tan(phi*math.pi/180))

        if plot==True:

            plt.figure(figsize=(13,9))
            plt.subplot(221)
            plt.plot(lms,heights_landing_gear)
            plt.xlabel('lm', fontsize=17)
            plt.ylabel('height landing gear', fontsize=17)
            plt.xticks(fontsize=15)
            plt.yticks(fontsize=15)

            plt.subplot(222)
            plt.plot(lms,y_mlgs)
            plt.xlabel('lm', fontsize=17)
            plt.ylabel('ymlg', fontsize=17)
            plt.xticks(fontsize=15)
            plt.yticks(fontsize=15)

            plt.subplot(223)
            plt.plot(heights_landing_gear, y_mlgs)
            plt.xlabel('height landing gear', fontsize=17)
            plt.ylabel('ymlg', fontsize=17)
            plt.xticks(fontsize=15)
            plt.yticks(fontsize=15)

            plt.subplot(224)
            plt.plot(lms, tip_engine_clearance)
            plt.xlabel('lm', fontsize=17)
            plt.ylabel('tip engine clearance', fontsize=17)
            plt.xticks(fontsize=15)
            plt.yticks(fontsize=15)
            plt.show()




def main():
    pass

if __name__ == "__main__":
    LandingGear=LandingGear(4000,3158, 0.19,0.47,0.8,0.06,0.47,0.0178,2,3.4,0.28,17.5,9,6.75,9.1,13.6,6,4.25,6.1)

    print('Designed maximum sustainable energy main landing gear', LandingGear.maximum_energy_sustainable_design)
    print('Main landing gear max energy encountered', LandingGear.energy_touch_down)
    print('Main landing gear required stroke length',LandingGear.required_stroke_length_main)
    print('Main landing gear diameter',LandingGear.diameter_landing_gear_main)
    print('Nose landing gear required stroke length',LandingGear.required_stroke_length_nose)
    print('Nose landing gear diameter',LandingGear.diameter_landing_gear_nose)

    print (LandingGear.lateral_tip_over(8.65,1.480,10.74,5))

