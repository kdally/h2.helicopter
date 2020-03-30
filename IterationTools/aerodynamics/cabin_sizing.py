import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

class cockpit_nose():

    """
    This class contains all the masses and volumes of the interior components of the aicraft located in the cocckpit/nose and cabin part of the aircraft.
    All the cg distances are given with respect to the nose.
    """

    def __init__(self,MTOW,OEW,width_fuselage,height_fuselage):

        # Torenbeek Method

        self.MTOW=MTOW #kg
        self.OEW = OEW  # kg

        self.seat_mass_correction_factor_confortability = 1.5

        self.mass_reduction_factor = 0.6 # mass reduction factor due to cable shrinkage, AFDX and component weight reduction

        self.number_cockpit_crew_members=1 # - number of pilots
        self.number_pilot_seats = 1  # -
        self.mass_pilot = 100  # kg
        self.mass_pilot_seat = 10*self.seat_mass_correction_factor_confortability   #kg 10kg, Roskam III

        self.xcg_pilot =  0.975# m

        self.mass_avionics_instruments = (40+0.008*self.MTOW)*self.mass_reduction_factor# kg
        self.xcg_avionics_instruments = 0.61 # m

        self.mass_radar_antenna_equipment = self.mass_avionics_instruments*0.05   # antenna size assumption
        self.xcg_radar_antenna =  0.2# m

        self.mass_air_conditioner=0.081*self.OEW*0.05 #kg reduction of 95% since no pressurization is needed
        self.xcg_air_conditioner = 0.390 # m

        self.mass_flight_control=0.64*((self.MTOW)**(2/3))*self.mass_reduction_factor/2 #kg divided by two since components for one pilot only are needed
        self.xcg_flight_control = 0.61  # m

        self.mass_wordrobe = 15  # kg
        self.xcg_wordrobe = 1.5 # m


        # self.volume_fan= # m^3
        # self.mass_fan =  # kg
        # self.xcg_fan = 0.2  # m

class passenger_cabin():

    def __init__(self,MTOW,OEW,width_fuselage,height_fuselage):

        self.MTOW=MTOW #kg
        self.OEW = OEW  # kg

        self.seat_mass_correction_factor_confortability=1.5

        self.number_passengers = 6  # -
        self.mass_passengers = 90  # kg
        self.xcg_passenger1= 2.5 #m
        self.xcg_passenger_pitch=1.100 #m

        self.number_baggage = 6  # -
        self.mass_baggage = 15  # kg
        self.xcg_cargo = 5.2  # m

        self.number_passenger_seats = 6  # -
        self.mass_passenger_seats=20*self.seat_mass_correction_factor_confortability #kg 20kg, Roskam III
        self.xcg_seat_pitch = 1.100  # m

        self.mass_frigo_bar_table = 5 #  kg
        self.number_frigo_bar_table = 6  # -
        self.xcg_frigo_bar_table1 = 1.950  # m
        self.xcg_frigo_bar_pitch = 1.050  # m

    def landing_gear_calculation(self, xcg, lm, ln):

        self.kg_to_lbs = 2.20462  # conversion from kg to lbs

        #main landing gear
        Kg_r=1 # gears are in the fuselage
        Ag_main=40
        Bg_main=0.16
        Cg_main=0.019
        Dg_main=1.5*10**(-5)
        mass_main_landing_gear =  Kg_r*(Ag_main+Bg_main*(self.MTOW*self.kg_to_lbs)**0.75+Cg_main*self.MTOW+Dg_main*(self.MTOW*self.kg_to_lbs)**1.5)  #  lbs
        mass_main_landing_gear=mass_main_landing_gear/self.kg_to_lbs
        xcg_main_landing = xcg+lm # m

        # nose landing gear
        Ag_nose=20
        Bg_nose=0.10
        Cg_nose=0
        Dg_nose=2*10**(-6)

        mass_front_landing_gear = Kg_r*(Ag_nose+Bg_nose*(self.MTOW*self.kg_to_lbs)**0.75+Cg_nose*self.MTOW+Dg_nose*(self.MTOW*self.kg_to_lbs)**1.5)  #  lbs
        mass_front_landing_gear=mass_front_landing_gear/self.kg_to_lbs
        xcg_front_landing = xcg-ln  # m

        return mass_main_landing_gear, xcg_main_landing, mass_front_landing_gear,xcg_front_landing

def xcg_fixed_equipment(cockpit_nose,passenger_cabin):

    mass_cockpit_nose=(cockpit_nose.mass_pilot+\
                    cockpit_nose.mass_pilot_seat+\
                    cockpit_nose.mass_radar_antenna_equipment+\
                    cockpit_nose.mass_air_conditioner+\
                    cockpit_nose.mass_avionics_instruments+\
                    cockpit_nose.mass_flight_control+cockpit_nose.mass_wordrobe)

    mass_cabin=(passenger_cabin.number_passengers*passenger_cabin.mass_passengers+\
                          passenger_cabin.number_passenger_seats*passenger_cabin.mass_passenger_seats+\
                          passenger_cabin.mass_baggage+passenger_cabin.number_baggage+\
                          passenger_cabin.number_frigo_bar_table*passenger_cabin.mass_frigo_bar_table+ \
                          passenger_cabin.number_baggage*passenger_cabin.mass_baggage)

    mass_fixed_equipment_cabin_cockpit=(cockpit_nose.mass_pilot_seat+\
                    cockpit_nose.mass_pilot+\
                    cockpit_nose.mass_radar_antenna_equipment+\
                    cockpit_nose.mass_air_conditioner+\
                    cockpit_nose.mass_avionics_instruments+\
                    cockpit_nose.mass_flight_control+ \
                    passenger_cabin.number_passenger_seats * passenger_cabin.mass_passenger_seats + \
                    passenger_cabin.number_frigo_bar_table * passenger_cabin.mass_frigo_bar_table
                    )

    xcg=((cockpit_nose.mass_pilot*cockpit_nose.xcg_pilot+ \
          cockpit_nose.mass_pilot_seat * cockpit_nose.xcg_pilot + \
          cockpit_nose.mass_pilot_seat*cockpit_nose.xcg_pilot+\
          cockpit_nose.mass_radar_antenna_equipment*cockpit_nose.xcg_radar_antenna+\
          cockpit_nose.mass_air_conditioner*cockpit_nose.xcg_air_conditioner+\
          cockpit_nose.mass_avionics_instruments*cockpit_nose.xcg_avionics_instruments+\
          cockpit_nose.mass_flight_control*cockpit_nose.xcg_flight_control+\
          cockpit_nose.mass_wordrobe+cockpit_nose.xcg_wordrobe)+
         (passenger_cabin.mass_baggage*passenger_cabin.number_baggage*passenger_cabin.xcg_cargo+\
          passenger_cabin.mass_passenger_seats*2*passenger_cabin.xcg_passenger1+ \
          passenger_cabin.mass_passenger_seats * 2 * (passenger_cabin.xcg_passenger1+passenger_cabin.xcg_seat_pitch) + \
          passenger_cabin.mass_passenger_seats * 2 * (passenger_cabin.xcg_passenger1 + 2*passenger_cabin.xcg_seat_pitch) + \
          passenger_cabin.mass_frigo_bar_table * 2 * passenger_cabin.xcg_frigo_bar_table1 + \
          passenger_cabin.mass_frigo_bar_table * 2 * (passenger_cabin.xcg_frigo_bar_table1 + passenger_cabin.xcg_frigo_bar_pitch) + \
          passenger_cabin.mass_frigo_bar_table * 2 * (passenger_cabin.xcg_frigo_bar_table1 + 2 * passenger_cabin.xcg_frigo_bar_pitch)))\
        /(mass_fixed_equipment_cabin_cockpit)

    mass_fixed_equipment_cabin_cockpit_list=[cockpit_nose.mass_pilot_seat, \
                    cockpit_nose.mass_pilot,\
                    cockpit_nose.mass_radar_antenna_equipment,\
                    cockpit_nose.mass_air_conditioner,\
                    cockpit_nose.mass_avionics_instruments,\
                    cockpit_nose.mass_flight_control, \
                    passenger_cabin.number_passenger_seats * passenger_cabin.mass_passenger_seats , \
                    passenger_cabin.number_frigo_bar_table * passenger_cabin.mass_frigo_bar_table]


    return xcg,mass_cockpit_nose,mass_cabin, mass_fixed_equipment_cabin_cockpit, mass_fixed_equipment_cabin_cockpit_list

def main():
    pass

if __name__ == "__main__":
    #Reference from Mitsubishi Mu-2B

    passenger_cabin=passenger_cabin(4000, 2600,1.480,1.4)
    print (passenger_cabin.landing_gear_calculation(4.98, 0.2,3.5))
    cockpit_nose=cockpit_nose(4000, 2600,1.480,1.4)

    print (xcg_fixed_equipment(cockpit_nose,passenger_cabin))

    print (cockpit_nose.mass_avionics_instruments*0.25)

    print ((passenger_cabin.number_passenger_seats * passenger_cabin.mass_passenger_seats + passenger_cabin.number_frigo_bar_table * passenger_cabin.mass_frigo_bar_table)*0.3*1.5)





