import numpy as np
import pickle
import Unpacker as un
import fpp.WeightFPPfinal as pp

class MassBudget():

    def __init__(self,oew_mass_budget,mass_fixed_equipment_cabin_cockpit_list, tank_mass_list, mass_nacelle_list,prop_mass_list):


# ======================== Fuselage  =====================================================================================================

        components=[]

        self.fuselage_mass=oew_mass_budget[0]

        self.fuselage_frame_mass=self.fuselage_mass*0.5
        components.append(self.fuselage_frame_mass)
        self.fuselage_floor_mass=self.fuselage_mass*0.15
        components.append(self.fuselage_floor_mass)
        self.fuselage_panels_mass=self.fuselage_mass*0.20
        components.append(self.fuselage_panels_mass)
        self.fuselage_interior_mass=self.fuselage_mass*0.15
        components.append(self.fuselage_interior_mass)

#======================== Fixed Equipment =====================================================================================================

        self.fixed_equipment = sum(mass_fixed_equipment_cabin_cockpit_list)

        self.pilot_mass=mass_fixed_equipment_cabin_cockpit_list[1]
        components.append(self.pilot_mass)

        self.radar_antenna_mass =mass_fixed_equipment_cabin_cockpit_list[2]
        components.append(self.radar_antenna_mass)

        self.air_condition_mass = mass_fixed_equipment_cabin_cockpit_list[3]
        components.append(self.air_condition_mass)

        self.avionics_mass=mass_fixed_equipment_cabin_cockpit_list[4]
        components.append(self.avionics_mass)

        self.flight_control_mass = mass_fixed_equipment_cabin_cockpit_list[5]
        components.append(self.flight_control_mass)

        self.furnishing_mass= mass_fixed_equipment_cabin_cockpit_list[0]+mass_fixed_equipment_cabin_cockpit_list[6]+mass_fixed_equipment_cabin_cockpit_list[7]
        components.append(self.furnishing_mass)

        self.hydraulics=175 # kg mass of hydraulics system
        components.append(self.hydraulics)

        self.electric_cables_mass =self.avionics_mass*0.25 # estimation from typical electrical systems mass roskam

#======================== Empennage =====================================================================================================

        self.horizontal_tail_mass=oew_mass_budget[2]

        ####put weight breakdown of horizontal tail and append to components
        self.horizontal_tail_skin_mass=0.41*self.horizontal_tail_mass
        components.append(self.horizontal_tail_skin_mass)
        self.horizontal_tail_spar_mass=0.38*self.horizontal_tail_mass
        components.append(self.horizontal_tail_spar_mass)
        self.horizontal_tail_stiffner_mass=0.21*self.horizontal_tail_mass
        components.append(self.horizontal_tail_stiffner_mass)


        self.vertical_tail_mass=oew_mass_budget[3]
        self.vertical_tail_skin_mass=0.52*self.vertical_tail_mass
        components.append(self.vertical_tail_skin_mass)
        self.vertical_tail_spar_mass=0.34*self.vertical_tail_mass
        components.append(self.vertical_tail_spar_mass)
        self.vertical_tail_stiffner_mass=0.14*self.vertical_tail_mass
        components.append(self.vertical_tail_stiffner_mass)

        ####put weight breakdown of vertical tail and append to components

#======================== Propeller =====================================================================================================

        self.propeller_mass=oew_mass_budget[5]
        components.append(self.propeller_mass)
        self.propeller_external_cover_mass=prop_mass_list[0]*0.1 # assumption
        self.propeller_honeycomb_mass = prop_mass_list[0] * 0.9 # assumption
        self.hub_propeller_mass=prop_mass_list[1]

#======================== Nacelle =====================================================================================================

        self.nacelle_mass=sum(mass_nacelle_list)


        self.electric_motor_mass=mass_nacelle_list[0]
        components.append(self.electric_motor_mass)
        self.motor_controller_mass=mass_nacelle_list[1]
        components.append(self.motor_controller_mass)
        self.gearbox_mass = mass_nacelle_list[2]
        components.append(self.gearbox_mass)
        self.tilt_mechanism_mass = mass_nacelle_list[3]
        components.append(self.tilt_mechanism_mass)
        self.nacelle_cover_mass=self.nacelle_mass *0.1 # assumption
        # components.append(self.nacelle_cover_mass)


#======================== Main Wing =====================================================================================================

        self.wing_mass=oew_mass_budget[7]

        self.wing_ribs_mass=self.wing_mass*0.15
        components.append(self.wing_ribs_mass)
        self.wing_upper_panel_mass=self.wing_mass*0.065
        components.append(self.wing_upper_panel_mass)
        self.wing_lower_panel_mass = self.wing_mass * 0.065
        components.append(self.wing_lower_panel_mass)

        self.wing_structure_part_mass=self.wing_mass*0.72

        self.wing_spar_mass=self.wing_structure_part_mass*0.035
        components.append(self.wing_spar_mass)
        self.wing_skin_mass = self.wing_structure_part_mass *0.035
        components.append(self.wing_skin_mass)
        self.wing_stringers_mass=self.wing_structure_part_mass*0.93
        components.append(self.wing_stringers_mass)

# ======================== Landing Gear =====================================================================================================

        self.nose_landing_gear_mass=oew_mass_budget[4]
        components.append(self.nose_landing_gear_mass)

        self.nose_landing_gear_strut_mass = self.nose_landing_gear_mass * 0.8  # assumption

        self.main_landing_gear_mass=oew_mass_budget[8]
        components.append(self.main_landing_gear_mass)

        self.main_landing_gear_strut_mass=self.main_landing_gear_mass*0.8 #assumption

#======================== Power Plant =====================================================================================================

        self.battery_mass=oew_mass_budget[9]
        components.append(self.battery_mass)

        self.tank_mass=sum(tank_mass_list)

        self.tank_inner_wall_mass=tank_mass_list[0]
        components.append(self.tank_inner_wall_mass)
        self.tank_mli_radiation_shield_mass=tank_mass_list[1]
        components.append(self.tank_mli_radiation_shield_mass)
        self.tank_mli_spacer_mass=tank_mass_list[2]
        components.append(self.tank_mli_spacer_mass)
        self.tank_outer_wall_mass=tank_mass_list[3]
        components.append(self.tank_outer_wall_mass)
        self.support_mass=tank_mass_list[4]
        components.append(self.support_mass)
        # self.tank_liner_mass=tank_mass_list[]
        # components.append()

        self.fuel_cell_mass=oew_mass_budget[11]
        components.append(self.fuel_cell_mass)
        self.radiator_mass=oew_mass_budget[12]
        components.append(self.radiator_mass)

        self.oew_mass_budget_sum=sum(oew_mass_budget)
        self.oew_mass_budget_components=sum(components)


    def calculate_material_mass(self):

 # ======================== Fuselage  =====================================================================================================


        self.fuselage_frame_mass_material = self.fuselage_frame_mass*1.5
        self.fuselage_floor_mass_material = self.fuselage_floor_mass*1.5
        self.fuselage_panels_mass_material = self.fuselage_panels_mass*1.5
        self.fuselage_interior_mass_material = self.fuselage_interior_mass*1.5

# ======================== Fixed Equipment==============================================================================================

        self.furnishing_mass_material = self.furnishing_mass*0.3*1.5
        self.electric_cables_mass_material = self.electric_cables_mass*1.5

# ======================== Empennage =====================================================================================================

        self.horizontal_tail_mass_material = self.horizontal_tail_mass*1.5

        self.vertical_tail_mass_material = self.vertical_tail_mass*1.5

# ======================== Propeller =====================================================================================================

        self.propeller_external_cover_mass_material =  self.propeller_external_cover_mass*1.5# assumption
        self.propeller_honeycomb_mass_material = self.propeller_honeycomb_mass*1.5  # assumption

# ======================== Nacelle =====================================================================================================

        self.nacelle_cover_mass_material = self.nacelle_cover_mass*1.5

# ======================== Main Wing =====================================================================================================

        self.wing_ribs_mass_material = self.wing_ribs_mass*1.5
        self.wing_upper_panel_mass_material = self.wing_upper_panel_mass*1.5
        self.wing_lower_panel_mass_material = self.wing_lower_panel_mass*1.5
        self.wing_structure_part_mass_material = self.wing_structure_part_mass*1.5
        self.wing_spar_mass_material = self.wing_spar_mass*1.5
        self.wing_skin_mass_material = self.wing_skin_mass*1.5
        self.wing_stringers_mass_material = self.wing_stringers_mass*1.5

# ======================== Landing Gear =====================================================================================================

        self.nose_landing_gear_strut_mass_material = self.nose_landing_gear_strut_mass*1.5

        self.main_landing_gear_strut_mass_material = self.main_landing_gear_strut_mass*1.5

# ======================== Power Plant =====================================================================================================

        self.tank_inner_wall_mass_material =self.tank_inner_wall_mass*1.5
        self.tank_mli_radiation_shield_mass_material = self.tank_mli_radiation_shield_mass*1.5
        self.tank_mli_spacer_mass_material =self.tank_mli_spacer_mass*1.5
        self.tank_outer_wall_mass_material =self.tank_outer_wall_mass*1.5
        # self.tank_liner_mass_material = self.tank_liner_mass*1.5


def main():
    pass

if __name__ == "__main__":

    def unpack(pickle_name):
        with open(pickle_name, 'rb') as handle:
            session = pickle.load(handle)
        for key in session:
            dir().append(key)
            globals()[key] = session[key]


    unpack('finaldata.pickle')

    MassBudget=MassBudget(oew_mass_budget,mass_fixed_equipment_cabin_cockpit_list, tank_mass_list, mass_nacelle_list,prop_mass_list)
    print()
    print('components equipment mass',oew_mass_budget[1])
    print('components equipment', sum(mass_fixed_equipment_cabin_cockpit_list)+175)
    print()

    print('mass nacelle list', MassBudget.nacelle_mass)
    print('mass nacelle sum', sum([MassBudget.electric_motor_mass, MassBudget.motor_controller_mass, MassBudget.gearbox_mass, MassBudget.tilt_mechanism_mass ]))
    print()

    print('mass tank list', sum([MassBudget.tank_inner_wall_mass,MassBudget.tank_mli_radiation_shield_mass,MassBudget.tank_mli_spacer_mass,MassBudget.tank_outer_wall_mass,MassBudget.support_mass]))
    print ('mass tank sum',oew_mass_budget[10])
    print()

    print ('wing mass sum', MassBudget.wing_mass)
    print('wing mass list', sum([MassBudget.wing_ribs_mass,MassBudget.wing_upper_panel_mass,MassBudget.wing_lower_panel_mass,MassBudget.wing_spar_mass,MassBudget.wing_skin_mass, MassBudget.wing_stringers_mass]))

    print (MassBudget.oew_mass_budget_sum)
    print(MassBudget.oew_mass_budget_components)

    print('final mass', MTOW)
    print('joining', MTOW*0.04)


    print()
    print('Material breakdown')
    print (MassBudget.furnishing_mass)

    print ('Thibault Masses')
    print (MassBudget.wing_mass)
    print (MassBudget.vertical_tail_mass+MassBudget.horizontal_tail_mass)
    print(MassBudget.fuselage_mass)
    print(MassBudget.main_landing_gear_mass+MassBudget.nose_landing_gear_mass)
    print(MassBudget.battery_mass+MassBudget.fuel_cell_mass+MassBudget.tank_mass+MassBudget.radiator_mass+MassBudget.propeller_mass+MassBudget.nacelle_mass)




