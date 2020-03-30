import pandas as pd
import numpy as np
from math import pi, floor
import matplotlib.pyplot as plt
import os

'''
Each of the airfoils contais:
    ref_deg, ref_cl, cl_0, cl_max, cd_0
'''


def cl_alpha(deg, cl):
    return cl/deg*(180/pi)

def max_alpha(cl_alpha, cl_max, cl_0):
    return (cl_max-cl_0)/cl_alpha

def alpha_0(cl_alpha, cl_0):
    return cl_0/cl_alpha


airfoildata = {'0009': [10, 0.9, 0, 1, 0.008],
               '0011': [5, 0.65, 0, 1, 0.013],
               '0012': [10, 1.1, 0, 1.2, 0.008],
               '0018': [10, 1.1, 0, 1.3, 0.015],
               '23015': [10, 1.15, 0.01, 1.5, 0.015],
               '64210': [5, 0.45, 0.2, 0.9, 0.01]}


class Airfoil():
    def __init__(self, name):
        self.t_c        = int(name)
        self.cl_alpha   = cl_alpha(airfoildata[name][0], airfoildata[name][1])
        self.cl_0       = airfoildata[name][2]
        self.cl_max     = airfoildata[name][3]
        self.cd_0       = airfoildata[name][4]
        self.alpha_max  = max_alpha(self.cl_alpha, self.cl_max, self.cl_0)
        self.alpha0     = alpha_0(self.cl_alpha, self.cl_0)
        pre = os.path.dirname(os.path.realpath(__file__))
        fname = 'airfoils.xlsx'
        path = os.path.join(pre, fname)
        self.info = pd.read_excel(path, sheet_name=name).values

        pre2 = os.path.dirname(os.path.realpath(__file__))
        fname2 = 'validation_viterna.xlsx'
        path2 = os.path.join(pre2, fname2)
        self.alfavit = pd.read_excel(path2).values[:,0]
        self.cdvit = pd.read_excel(path2).values[:,1]



class Rotor():
    def __init__(self, N, B, R, cutout, taper, solidity, theta_tip, airfoil, v_inf, be, h):
        self.N      = N                 # [-]   Number of rotors
        self.B      = B                 # [-]   Number of blades
        self.cutout = cutout            # [-]   percentage of blade at the root that does not produce lift
        self.sigma  = solidity          # [-]   Rotor solidity
        self.R0     = self.cutout*R     # [m]   Begin of blade
        self.R      = R                 # [m]   Blade radius
        self.taper  = taper             # [-]   Taper ratio
        # self.c0     = self.sigma*pi*R/N # [m]   Blade chord
        self.ctip   = 2*pi*self.R*self.sigma/(self.B*(self.taper+1)) # [m]  Tip blade chord
        self.croot  = self.ctip*self.taper                           # [m]  Root blade chord
        # self.ctip   = 0.2358
        # self.croot  = 0.2358
        # self.theta_0 = theta_0          # [rad] Root twist (at R0)
        self.theta_tip = theta_tip      # [rad] Tip twist (at R)
        self.M_tip  = 0.75              # [-]   Maximum tip Mach number
        self.v_inf  = v_inf             # [m/s] Undisturbed flow speed
        self.be     = be               # [-]   Number of blade elements

        self.area   = self.N*np.pi*(self.R**2-self.R0**2)        # [m2]   rotor area
        self.area   = self.N*np.pi*self.R**2
        self.airfoil = Airfoil(airfoil)

        self.h = h
        self.ISA()
        self.nonDimensionalCoeff()
        self.calcTwist('flat', 0)
        # print('ctip', self.ctip)
        # print('croot', self.croot)

    def ISA(self):
        self.temp = 288.15-0.0065*self.h
        self.pressure = 101325*pow(self.temp/288.15,(-9.80665/(-0.0065*287.05)))
        self.rho = self.pressure/(287.05*self.temp)
        self.v_sound = np.sqrt(1.4 * 287.05 * self.temp)
        self.v_tip_max = self.M_tip * self.v_sound
        self.maxOmega = self.v_tip_max/self.R
        print('v_tip', self.v_tip_max)

    def nonDimensionalCoeff(self):
        self.r = np.linspace(self.cutout, 1, self.be)[:-1]       # Non dimensional stations at blade elements
        # print(self.r)
        # self.r = np.linspace(self.cutout, 1, self.be+1)[1:]       # Non dimensional stations at blade elements
        # print(self.r)
        self.dr = (self.r[1]-self.r[0])
        self.c = self.croot + (self.croot-self.ctip)/(self.R0-self.R)*(self.R*self.r-self.R0)
        # self.c = (self.c0 - (self.c0 - self.ctip) / (self.cutout*self.R)*self.r).self.R # Non dimensional chords at blade elements
        # self.M = self.r*self.omega/self.v_sound

    def calcTwist(self, name, *args):
        if name == 'ideal':
            self.twist = args[0]/self.r*args[1]

        elif name == 'linear':
            theta_root = args[0]
            theta_tip  = args[1]
            # slope = (theta_tip-theta_root)/(self.R-self.R0-self.dr*self.R)
            slope = (theta_tip-theta_root)/(1-self.cutout)
            self.twist = theta_root + slope*(self.r-self.cutout)
            # self.twist = theta_root + slope *(self.R*self.r - self.R0)

        elif name == 'poly':
            if len(args) == 3:
                self.twist = args[0] + args[1]*self.r + args[2]*self.r**2
            if len(args) == 4:
                self.twist = args[0] + args[1]*self.r + args[2]*self.r**2 + args[3]*self.r**3

        elif name == 'coords':
            self.twist = args[0]

        elif name == 'flat':
            self.twist = np.array([args[0]]*len(self.r))
        # plt.plot(self.r, self.twist)
        # plt.show()


    # def calcCollective(self, thrust, omega):
    #     v_i = np.sqrt(thrust / (2*self.rho*self.area))
    #
    #     n1 = 2*thrust/(self.airfoil.cl_alpha*self.rho*self.c0*self.N*self.B)
    #     n2 = (self.theta_tip*self.R*v_i**2 - v_i**3/omega)/4*np.log(1/self.cutout)
    #     n3 = (self.theta_tip*self.R*omega**2 - v_i*omega)*(1-self.cutout**2)*self.R**2/2
    #     d1 = v_i**2/4*(1-self.cutout)*self.R + omega**2*(1-self.cutout**3)*self.R**3/3
    #
    #     alpha_collective = (n1-n2-n3-d1*self.airfoil.cl_0)/d1
    #
    #     alphas = alpha_collective + self.airfoil.alpha0 + self.theta_tip/self.r - v_i/(omega*self.r*self.R)
    #
    #     # print(self.r)
    #     # print(alphas)
    #     # print(alphas*self.airfoil.cl_alpha)
    #
    #
    #     condition = self.theta_tip*self.R - v_i/omega
    #     alpha_max = np.where(condition>0, alpha_collective+self.airfoil.alpha0+ \
    #                          (self.theta_tip*self.R - v_i/omega)/self.R0,
    #              alpha_collective+self.airfoil.alpha0+ (self.theta_tip*self.R - v_i/omega)/self.R)
    #
    #     # print()
    #
    #     # print(alpha_collective*180/pi)
    #     # print(alpha_max*self.airfoil.cl_alpha)
    #     # print(alpha_max)
    #
    #     alpha_collective = np.where(alpha_max > self.airfoil.alpha_max, 1000, alpha_collective)
    #
    #     # print(alpha_collective)
    #
    #     return alpha_collective
    #
    # def powerInduced(self):
    #     v_induced = np.sqrt()
    #
    # def calcPProfile(self, h, thrust):
    #     temp, rho, v_sound = self.ISA(h)
    #     v_i = np.sqrt(thrust / (2*rho*self.area))


    def calcF(self, infl):
        phi = np.arctan(infl / self.r)
        f = self.B/2*(1-self.r)/(self.r*np.sin(phi))
        # f = self.B/2*(1/(self.r*phi)-1/phi)
        F = 2/pi*(np.arccos(np.e**(-f)))
        return F

    def calcF2(self, infl):
        phi = np.arctan(infl / self.r)
        f_tip = self.B/2*(1-self.r)/(self.r*np.sin(phi))
        F_tip = 2/pi*(np.arccos(np.e**(-f_tip)))

        f_hub = self.B*(self.r - self.cutout + 0.000001)/(2*self.cutout*np.sin(phi))
        F_hub = 2/pi*(np.arccos(np.e**(-f_hub)))
        # F_hub[0] = 0.0001
        # if F_hub[0] == 0:
        #     F_hub[0] = 0.00001
        # print(F_tip*F_hub)
        return F_tip*F_hub

    def calcF3(self, infl, idx):
        new_r = self.r[idx]
        phi = np.arctan(infl / new_r)
        f_tip = self.B/2*(1-new_r)/(new_r*np.sin(phi))
        F_tip = 2/pi*(np.arccos(np.e**(-f_tip)))

        f_hub = self.B*(new_r - new_r[0])/(2*new_r[0]*np.sin(phi))
        F_hub = 2/pi*(np.arccos(np.e**(-f_hub)))
        # F_hub[0] = 0.0001
        # if F_hub[0] == 0:
        #     F_hub[0] = 0.00001
        # print(F_tip*F_hub)
        return F_tip*F_hub

    def inflow(self, F, inf, col):
        theta = self.twist + col
        n1 = np.sqrt((self.sigma*self.airfoil.cl_alpha/16/F - inf/2)**2 + \
                     self.sigma*self.airfoil.cl_alpha/8/F*theta*self.r)
        n2 = -(self.sigma*self.airfoil.cl_alpha/16/F)+inf/2

        return n1 + n2

    def inflow2(self, F, inf, inflow, idx, col):
        theta = self.twist[idx] + col
        n1 = inf/2
        aoa_r = theta-np.arctan(inflow[idx]/self.r[idx])
        phi = theta-aoa_r
        closest_aoa_r = np.round(aoa_r*180/pi *4)/4
        index = ((closest_aoa_r - self.airfoil.info[0, 0])*4-1).astype(int)
        cls = self.airfoil.info[index, 1]*np.cos(phi) - self.airfoil.info[index, 2]*np.sin(phi)
        # print('cls')
        # print(cls)
        if inf == 0:
            idx = np.where(cls > 0)
            new_r = self.r[idx]
        else:
            inside = inf**2/4+self.sigma*cls[idx]*self.r[idx]/(8*F[idx])
            idx = np.where(inside > 0)
            # r_min = -inf**2*2/(self.sigma*cls)
            # idx = np.arange(idx[-1]+1, self.be, 1)
            new_r = self.r[idx]
            # new_r = self.r[idx]
            # idx = np.arange(0, self.be-1)
            # new_r = self.

        # new_cutout = self.r[new_cutout[-1]]
        # print((np.isnan(closest_aoa_r)))
        # print('collective')
        # print(col)
        # if min(closest_aoa_r) > 10:
        #     print('col problem')
        # if np.any(np.isnan(closest_aoa_r)) or min(closest_aoa_r) > 10:
        #     print('col and nan problem')
        #     return 1000

        # phi = theta-aoa_r
        # index = ((closest_aoa_r - self.airfoil.info[0, 0])*4-1).astype(int)

        # cls = self.airfoil.info[index, 1]*np.cos(phi) - self.airfoil.info[index, 2]*np.sin(phi)
        # cls = self.airfoil.info[index, 1]

        n2 = np.sqrt(inf**2/4+self.sigma*cls[idx]*new_r/(8*F[idx]))
        # print('shape end')
        # print(np.shape(n2))
        return  n1+n2, idx

    def iterations(self, omega, collective):
        F1 = np.array([1]*len(self.r))

        inf = self.v_inf / omega / self.R

        iterating = True
        # print('here')
        while iterating:
            infl = self.inflow(F1, inf, collective)
            F2   = self.calcF(infl)
            # print('inflow')
            # print(infl)
            # print('F2')
            # print(F2)
            diff = F1-F2
            if abs(np.average(diff)) > 0.002:
                # print('print diff rotor F2')
                # print(np.average(diff))
                # print(count)
                F1 = (F2+F1)/2
            else:
                iterating = False

        return infl, F2

    def iterations2(self, omega, collective):
        F1 = np.array([1]*len(self.r))
        infl1 = np.array([0]*len(self.r))

        inf = self.v_inf / omega / self.R
        idx1 = F1

        iterating = True
        while iterating:
            # print('inflow')
            # print(infl1)
            # infl2, idx2 = self.inflow2(F1, inf, infl1, idx1, collective)
            infl2 = self.inflow(F1, inf, collective)
            # print(infl2)
            # if type(infl2) == int:
            #     return 1000, 1000
            # print('idx')
            # print(idx2)
            # F2   = self.calcF3(infl2, idx2)
            F2   = self.calcF2(infl2)
            # print('inflow2')
            # print(infl2)
            # print('F2')
            # print(F2)
            if len(F1) == len(F2):
                diff = (F1-F2)
                # if np.any

                if abs(np.average(diff)) > 0.002:
                    print(np.average((F1-F2)))
                    F1 = (F2+diff/2)
                    infl1 = infl2

                else:
                    iterating = False
            else:
                F1 = F1[-len(F2):]
                infl1 = infl1[-len(F2):]
            # print(self.r[idx2])
        # print(F2)
        # print(infl2)
        # plt.plot(self.r, F2)
        # plt.show()
        return infl2, F2
        # return infl2, F2, idx2


    def calcCtCp(self, col, inflow, F, omega):
        theta = self.twist + col
        # dct = self.sigma*self.airfoil.cl_alpha/2*(theta*self.r**2 - inflow*self.r)
        dct = 4*F*inflow*(inflow-self.v_inf/omega/self.R)*self.r
        ct = np.trapz(dct, self.r)+dct[-1]*self.dr

        # if 0.0042 < ct <0.044:
        #     plt.plot(self.r, dct)#*self.dr)

        aoa_r = theta-inflow/self.r
        closest_aoa_r = np.round(aoa_r*180/pi *4)/4
        index = ((closest_aoa_r - self.airfoil.info[0, 0])*4-1).astype(int)
        # if np.any(np.abs(index) > len(self.airfoil.info[:, 0])-1):

        # np.where(index)
        # plt.plot(self.r, aoa_r*180/pi)
        # plt.show()
        # if np.any(index > len(self.airfoil.info[:, 0])-1) or np.any(index < 0):
        if np.any(index > len(self.airfoil.info[:, 0])-1) or np.any(index == -2147483648) or np.any(index < -200):
            # print('index', index)
            # print('twist')
            # print(self.twist)
            # print('col')
            # print(col)
            # print('aoa')
            # print(aoa_r*180/pi)
            # print('')
            # print(np.max(index), np.min(index), len(self.airfoil.info[:, 0]-1))
            return 1000, 1000, 1000, 1000
        # print(index)
        phi = theta-aoa_r
        Mach = np.sqrt((omega*self.r*self.R)**2+self.v_inf**2)/self.v_sound
        # print(Mach)
        # if self.cond:
        #     ext = np.sqrt(1-Mach**2)
        # else:
        #     ext = 1
        cd = np.array(self.airfoil.info[index, 2])/np.sqrt(1-Mach**2)
        # cd = self.airfoil.info[index, 1]*np.sin(phi) - self.airfoil.info[index, 2]*np.cos(phi)
        # cd = 0.011
        # idx = 180+aoa_r
        # print(self.airfoil.cdvit[np.floor(idx)])
        # print(np.floor(idx))
        # cd = self.airfoil.cdvit[np.floor(idx)] + idx-np.floor(idx)*(self.airfoil.cdvit[np.floor(idx)+1]-self.airfoil.cdvit[np.floor(idx)])



        dcp = inflow*dct*self.dr + 0.5*self.sigma*cd*self.r**3*self.dr
        cp = np.sum(dcp)
        # cp = np.trapz(inflow, dct) + 0.5*self.sigma*cd*((self.R/self.R)**4-(self.R0/self.R)**4)/4
        # print('ct', ct)
        # print('cp', cp)
        # print('dct', dct)
        return ct, cp, dct, dcp

    def thrust(self, ct, omega):
        # return ct*self.rho*pi*self.R**2*(omega*self*R)**2
        return ct*self.rho*pi*self.R**2*(omega*self.R)**2

    def power(self, cp, omega):
        # return cp*self.rho*pi*self.R**2*(omega*self.R)**2*self.R
        return cp*self.rho*pi*self.R**2*(omega*self.R)**3

    def simulation(self, omega, col):
        inflow, F = self.iterations(omega, col)
        # print(self.calcCtCp(col, inflow, F, omega))
        # plt.plot(self.r, inflow)
        # plt.show()
        # print('inflow', inflow)
        # print('F', F)
        ct, cp, dct, dcp = self.calcCtCp(col, inflow, F, omega)
        # print('ct, cp', ct, cp)
        if ct == 1000:
            return 1000, 1000, 1000, 1000, 1000, 1000, 1000
        thrust = self.thrust(ct, omega)
        P_indu = self.power(cp, omega)
        return thrust, P_indu, ct, cp, inflow, dct, dcp

    def simulation2(self, omega, col):
        # inflow, F, idx = self.iterations2(omega, col)
        inflow, F = self.iterations2(omega, col)
        # print('col')
        # print(col)
        if type(inflow) == int:
            return 1000, 1000, 1000, 1000, 1000
        # print(self.calcCtCp(col, inflow, F, omega))
        # plt.plot(self.r, inflow)
        # plt.show()
        # print('inflow', inflow)
        # print('F', F)
        ct, cp, dct, dcp = self.calcCtCp(col, inflow, F, omega)
        # print('ct, cp', ct, cp)
        if ct == 1000:
            return 1000, 1000, 1000, 1000, 1000
        thrust = self.thrust(ct, omega)
        P_indu = self.power(cp, omega)
        # return thrust, P_indu, ct, cp, inflow, idx
        return thrust, P_indu, ct, cp, inflow, F


    # def spin(self, h, MTOW, omega):
    #     self.ISA(h)
    #     self.nonDimensionalCoeff()
    #     self.calcCollective(MTOW, omega)
    #     # print(self.calcCollective(MTOW, omega))


