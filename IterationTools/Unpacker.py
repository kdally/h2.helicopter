import pickle
import numpy as np
import fpp.WeightFPPfinal as pp
import load_calculation.wing_sizing as ws

def convertTuple(tup):
    tup = list(tup)
    for i in range(len(tup)):
        tup[i] = float(tup[i])
    return tuple(tup)


def unpack(pickle_name):
    with open(pickle_name, 'rb') as handle:
        session = pickle.load(handle)
    for key in session:
        dir().append(key)
        globals()[key] = session[key]

def change(a):
    a = 3
    return

unpack('finaldata_prop_changes.pickle')
print(battery)
# print(mission_profile[10000])
mass_powerplant, mass_fuelcell, mass_DCDC, mass_compressor, mass_radiator, mass_electric_motor, mass_motor_controller, mass_tank, mass_fuel, mass_battery, m_tankinnerwall, m_radshield, m_spacer, m_out, m_support = convertTuple(pp.calculateWeightFPP(mission_profile+11, mission_time, np.around(alex_data[7][4], decimals=1), np.around(alex_data[10][4], decimals=1), True))
print(mass_battery)
print(ws.wingstructure(M_h,Sy_h,T_h,b))
