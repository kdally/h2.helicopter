import math
import numpy as np
class Radiator():

    def __init__(self):

        self.t = 32 * 2  # [mm]
        self.l = 454  # [mm]
        self.w = 212  # [mm]
        self.m = 2.51 # [kg]
        self.coolant_volume = 508 * 1e-6 # [m3]
        self.coolant_flow = 20 * 1.6667e-5  # [m3/s]
        self.coolant_dp = 0.55  # [bar]
        self.max_v_flow = 0.25  # [m3/s]
        self.max_QT = 175  # [W/dC]

        self.T_cool_in    = 90   # [degC]
        self.cp_air       = 1006 # [J/kg K]@ ISA SL
        self.density_air  = 1.226# [kg/m3] @ ISA SL
        self.cp_egw       = 850  # [J/kg K]@ 70degC
        self.density_egw  = 1042 # [kg/m3] @ 80degC

    def HTcoeff(self,v_flow):  # heat transfer coefficient [W/dC] (volumetric flow rate [m3/s])
        # if v_flow < 0.25:
        #     return self.max_QT/self.max_v_flow * v_flow
        # else:
        #     return self.max_QT

        return  np.where(v_flow < 0.25, self.max_QT/self.max_v_flow * v_flow, self.max_QT)

    def v_flow(self,QT):
        # if QT >= self.max_QT:
        #     raise AssertionError
        # else:
        #     return QT/self.max_QT * self.max_v_flow

        return np.where(QT < self.max_QT, QT/self.max_QT * self.max_v_flow, QT)

    def air_dp(self,v_flow):  # Air side pressure drop [Pa]  (volumetric flow rate [m3/s])
            return 45/0.1333 * v_flow

    def n_rad_calc(self, Q, T_air_in):
        heat_needed = Q/2
        QT = self.HTcoeff(self.max_v_flow)

        # heat transfer per unit area of radiator core
        heat_per_rad = QT * (self.T_cool_in - T_air_in)/1000
        n_rad = np.ceil(heat_needed / heat_per_rad)

        tot_area = n_rad *(self.w/1000 * self.l/1000)

        return n_rad, tot_area

    def v_flow_calc(self, Q, T_air_in, n_rad):
        heat_needed = Q/2
        heat_out_per_rad = heat_needed *1000 / n_rad
        QT_needed = heat_out_per_rad / (self.T_cool_in - T_air_in)
        v_flow_needed = self.v_flow(QT_needed)

        return v_flow_needed * n_rad

    def in_out_calc(self, v_flow, P_inf, v_e, air_density, n_rad):

        # calculate total pressure of rotor wake
        P_tot_rot = P_inf + (1/2)*air_density*(v_e*v_e)
        # print('ptotrot', P_tot_rot)
        # print(P_inf + DP_rot/2)

        # calculate inlet velocity and area
        v_in = v_e / 2
        area_in = v_flow / v_in

        # calculate total pressure after radiator
        P_tot_2 = P_tot_rot - self.air_dp(v_flow/n_rad)

        # calculate outlet velocity and area
        v_out = np.sqrt((P_tot_2 - P_inf) * 2 / air_density)
        area_out = v_flow / v_out

        return area_in, area_out

    def final_weight(self, n_rad):
        mass_rad = self.m * n_rad
        mass_cool = self.coolant_volume * self.density_egw * n_rad

        return mass_cool + mass_rad

    def final_dimensions(self, n_rad):
        length = n_rad*self.w
        width  = self.l
        return length, width

    def n_rad_calc_2(self, Q, T_air_in):
        heat_needed = Q / 2
        nrad = np.ceil(heat_needed/75)
        tot_area = nrad*(1.45*0.375)
        return nrad, tot_area

    def final_weight_2(self, n_rad):
       mass_rad = n_rad*33.45
       return mass_rad

if __name__ == '__main__':
    rad2 = Radiator()
    # # calculate radiator area at take off
    # n_rad, area_tot = rad2.n_rad_calc(np.array([300, 200]), 50)
    # print('number radiator',n_rad)
    # print('total rad area', area_tot)
    # # Calculate volumetric flow needed at cruise
    # v_flow_tot = rad2.v_flow_calc(300, 30, n_rad)
    # print('total v_flow', v_flow_tot)
    #
    # a_in, a_out = rad2.in_out_calc(v_flow_tot, 101325, 100, 1.000, n_rad)
    #
    # print('area_in', a_in, 'area_out', a_out)
    # mass = rad2.final_weight(n_rad) * 2
    # print('total mass', mass)
    print(rad2.air_dp(15))