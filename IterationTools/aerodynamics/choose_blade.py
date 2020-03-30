import pandas as pd
import numpy as np
import os
import pickle
import progressbar
from IterationTools.fpp.WeightFPPfinal import calculateWeightFPP
import matplotlib.pyplot as plt


pre = os.path.dirname(os.path.realpath(__file__))
# fname = 'mission_profiles2.xlsx'
fname = 'mission_profiles_fin.xlsx'
path = os.path.join(pre, fname)

mission_profiles = pd.read_excel(path, sheet_name='power').values

pre = os.path.dirname(os.path.realpath(__file__))
fname = 'mission_data.xlsx'
path = os.path.join(pre, fname)

mission_data = pd.read_excel(path, sheet_name='mission_data').values[:,-1]
# print(mission_data)


masses    = []
Mintots   = []
massFCs   = []
massDCDCs = []
massCs    = []
massRADs  = []
massEMs   = []
massMCs   = []
tanks     = []
fuel2s    = []
wbefinas  = []
idxs      = []



t_ho    = mission_data[1]
t_cl_v  = mission_data[3]-t_ho
t_cl_h  = mission_data[7]-t_ho-t_cl_v
t_cr    = mission_data[8]-t_ho-t_cl_v-t_cl_h
t_de    = 0

TOe = t_ho+t_cl_v+t_cl_h

Tendcruise = TOe+t_cr+t_de
switch = False

# time = np.arange(0, t_ho+t_cl_v+t_cl_h+t_cr+t_de+t_ho, 0.1)

# time = np.arange(0.1, mission_data[-1], 0.1)
time = np.arange(0.1, np.shape(mission_profiles)[0]/10+0.1, 0.1)



# plt.plot(time, mission_profiles[:, -12])
# plt.show()

bar = progressbar.ProgressBar(maxval=len(mission_profiles[0]), widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
bar.start()

# plt.plot(time, mission_profiles[:,20], label='linear, little twist', linewidth=3)
# plt.plot(time, mission_profiles[:, 30], label='linear, medium twist', linewidth=3)
# plt.plot(time, mission_profiles[:, 39], label='linear, high twist', linewidth=3)
# plt.ylabel('Power [kW]')
# plt.xlabel('Time [s]')


for idx, _ in enumerate(mission_profiles[0]):
    mission_profile = mission_profiles[:, idx]
    # print('shape time')
    # print(np.shape(time))
    # print('shape mission profile')
    # print(np.shape(mission_profile))
    # if idx ==0:
    #     label_i = 'ideal, little twist'
    # if idx == 10:
    #     lablel_i = 'ideal, medium twist'
    # if idx == 20:
    #     label_i = 'linear, little twist'
    # if idx == 30:
    #     label_i = 'linear, medium twist'
    # if idx%10==0:
    #     plt.plot(time, mission_profile, label= label_i, linewidth=3)
    # plt.show()

    bar.update(idx+1)
    #
    mass_powerplant, mass_fuelcell, mass_DCDC, mass_compressor, mass_radiator, mass_electric_motor, mass_motor_controller, mass_tank, mass_fuel, mass_battery, m_tankinnerwall, m_radshield, m_spacer, m_out, m_support = calculateWeightFPP(mission_profile, time, np.around(mission_data[7], decimals=1), np.around(mission_data[10], decimals=1),switch)

    Mintots.append(mass_powerplant)
    # massFCs.append(mass_fuelcell)
    # massDCDCs.append(mass_DCDC)
    # massCs.append(mass_compressor)
    # massRADs.append(mass_radiator)
    # massEMs.append(mass_electric_motor)
    # massMCs.append(mass_motor_controller)
    # tanks.append(mass_tank)
    # fuel2s.append(fuel2s)
    # wbefinas.append(wbefinas)
    idxs.append(idx)



# plt.plot(idxs, Mintots)
# plt.axvline(x=t_ho)#, label='ho')
# plt.axvline(x=t_cl_v+t_ho)#, label='cl_v')
# plt.axvline(x=t_cl_h+t_cl_v+t_ho)#, label='cl_h')
# plt.axvline(x=t_cr+t_cl_h+t_cl_h+t_ho-15)#, label='cr')
# plt.axvline(x=t_de+t_cr+t_cl_h+t_cl_h+t_ho-15)#, label='de')
# plt.axvline(x=TOe)#, label='TOe')
# plt.axvline(x=Tendcruise)#, label='end cruise')
# plt.legend()

# plt.figure(2)
plt.plot(idxs, Mintots)
plt.show()


# mass_dict = {'Mintots':Mintots, 'massFCs':massFCs, 'massDCDCs':massDCDCs,
#              'massCs':massCs, 'massRADs':massRADs, 'massEMs':massEMs,
#              'massMCs':massMCs, 'tanks':tanks, 'fuel2s':fuel2s,
#              'wbefinas':wbefinas, 'idxs': idxs}

mass_dict = {'Mintots':Mintots, 'idxs': idxs}

pre = os.path.dirname(os.path.realpath(__file__))
fname = 'masses_fin.pickle'
path = os.path.join(pre, fname)
with open(path, 'wb') as handle:
    pickle.dump(mass_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)

# print(min(masses))
# print(masses.index(min(masses)))
