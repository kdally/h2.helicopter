##
import numpy as np
#from IterationTools import Mission
from fpp import Componentslibrary
from fpp import Componentspicker

from fpp import Lh2Tank
import matplotlib.pyplot as plt
from fpp import Radiator
##
def main():
    pass
if __name__ == "__main__":
    import numpy as np
    import warnings

    warnings.filterwarnings("ignore")
    import math
    from IterationTools import Mission
    import IterationTools.aerodynamics.RotorPerformance  as rp
    from IterationTools.aerodynamics.BET import Rotor
    from IterationTools.fpp import Componentslibrary
    from IterationTools.fpp import Componentspicker
    from IterationTools.fpp import Lh2Tank
    import matplotlib.pyplot as plt
    from IterationTools.fpp import Radiator
##
# OUTPUT= TOT MASS, mass fuel cell, mass DCDC , mass compressor, mass radiator, mass electric motor, mass motor controller, tank mass , fuel mass, mass battery
def calculateWeightFPP(powert, time, TOe, tendcruise,switch):
    print(max(powert))
    powert = np.around(powert, decimals=1)
    print(powert[11000])
    time = np.around(time, decimals=1)
    Pmin = 1
    Toend = np.where(time == TOe)[0]
    Pmax = np.amax(powert)
    time = time - 0.1
    dc = 10
    end2 = np.where(time == tendcruise)[0]
    teto = np.array(time[::dc])
    Pmfc = np.linspace(Pmin, Pmax, len(teto))
    DataFC, weightFC, nameFC = Componentslibrary.libfuelcell()
    name2FC, numberFC, massFC, ratingFC = Componentspicker.choosecomponentsarray(Pmfc, DataFC, weightFC, nameFC)
    DataEM, weightEM, nameEM = Componentslibrary.libelectricmotor()
    nameEM2, numberEM, massEM, ratingEM = Componentspicker.choosecomponents(Pmax/2, DataEM, weightEM, nameEM)
    DataMC, weightMC, nameMC = Componentslibrary.libmotorcontroller()
    nameMC2, numberMC, massMC, ratingMC = Componentspicker.choosecomponents(ratingEM, DataMC, weightMC, nameMC)
    DataDCDC, weightDCDC, nameDCDC = Componentslibrary.libdcdcconverter()
    # nameDCDC2, numberDCDC, massDCDC, ratingDCDC = Componentspicker.choosecomponentsarray(ratingFC, DataDCDC, weightDCDC, nameDCDC)
    massDCDC=3*weightDCDC[0]+16*2
    #print(numberDCDC)
    numberDCDC=3
    #massDCDC=massDCDC*numberDCDC
    if switch:
        plt.plot(Pmfc, ratingFC)
        plt.plot(Pmfc, np.ones(len(Pmfc))*343.0)
        plt.xlabel('Power to be delivered by fuel cell [kW]')
        plt.ylabel('Actual fuel cell stack power output [kW]')
        plt.show()
    Pmfc=ratingFC
    airmassflow=3.57E-4*1.6*Pmfc/0.7
    DataC, weightC, nameC = Componentslibrary.libcompressor()
    name2C, numberC, massC, ratingC = Componentspicker.choosecomponentsarray(airmassflow, DataC, weightC, nameC)
    massC = massC *(numberC+1)/numberC
    numberC=numberC+1

    Pmfc = np.array([Pmfc])
    powertTO = np.array([powert[::dc]])
    teto2, Pmfc2 = np.meshgrid(teto, Pmfc)
    RHS = (powertTO - Pmfc2) * (teto[2] - teto[1])
    Toend = np.around(Toend / 10, decimals=0)
    Toend = Toend.astype(int)[0]
    end2 = np.around(end2 / 10, decimals=0)
    end2 = end2.astype(int)[0]
    RHS[:, Toend:end2] = RHS[:, Toend:end2] + 0.3 * Pmfc.T
    RHS[RHS == 0] = 1e-10
    RHS[RHS < 0] = 1e-10
    RHS2 = RHS
    RHS = np.cumsum(RHS, axis=1)
    RHS = RHS[:, -1]
    res = RHS
    activitytime = np.ones(len(res))
    for i in range(len(res) - 1):
        activitytime[i] = teto[[j for j, e in enumerate(RHS2[i, :]) if e != 0][-1]]
    #activitytime[-1] = 0
    activitytime[activitytime==time[-1]] = 0
    # activitytime = np.array([activitytime])
    # teto = np.array([teto])
    # res = np.array([res])
    nmotor=0.96
    nbattery=0.96
    eff = 0.48#0.477 #0.54
    SEh2 = 141.8  # MJ/kg
    mjtokwh = 0.277778
    SEmp = 1.73 + 1  # kg/kg
    SPmp = 0.1#0.1102#0.1262 #0.5682  # kg/kw
    rad = Radiator.Radiator()
    n_rad, area_tot= rad.n_rad_calc_2(Pmfc*(1-eff)/eff +Pmax * np.ones(len(teto))*(1-nmotor)/nmotor+(Pmax - Pmfc)*(1-nbattery)/nbattery,  50)
    massRAD = rad.final_weight_2(n_rad)
    #SPcc = 0.195  # kg/kw
    SEb = 0.448#0.12#0.460
    SPb = 4.48#3
    numberFCcorrected=numberFC
    #PSF = [ 2 if k==1 else k/(k-1) for k in numberFC ]#numberFC/(numberFC-1) if numberFC>1 else 2 #1.3 #numberFC/(numberFC-1)
    #PSF=1.3
    wbe1 = (RHS) * mjtokwh * 1E-3 / SEb
    wbe2 = (Pmax - Pmfc) / SPb
    wbe1 = np.array([wbe1])
    #wbe22 = wbe2 * PSF
    wbefinal = np.maximum(wbe1, wbe2)
    #wbefinal2 = np.maximum(wbe1, wbe22)
    energymp = Pmfc * (TOe + time[-1] - tendcruise + 0.7 * (tendcruise - TOe))

    fuel= energymp / eff / (SEh2 * 1E3)
    #tank, fuel2=[Lh2Tank_beta.designMechanicalTank(fuel[0,j], Pmfc[0, j]) for j in range(len(fuel.T)-1)]
    # tank, fuel2 = [Lh2Tank_beta.designMechanicalTank(fuel[0, j], Pmfc[0, j])[0] for j in range(len(fuel.T) )], [
    #     Lh2Tank_beta.designMechanicalTank(fuel[0, j], Pmfc[0, j])[1] for j in range(len(fuel.T))]
    # tank=np.array([tank])
    tank, fuel2, m_tankinnerwall, m_radshield, m_spacer, m_tankouterwall, m_tanksupport = Lh2Tank.designMechanicalTank(fuel, Pmfc) #,
    #tank = np.array([tank])
    #fuel2=np.array([fuel2])
    wemp=tank+fuel2
    #wemp=np.array([wemp])
    #wemp = energymp / eff / (SEh2 * 1E3) * SEmp
    wpmp = massRAD*2 + massFC + massDCDC + massC #Pmfc * SPmp
    #wpmp2 = wpmp * PSF
    wpcc= (massEM+ massMC) * np.ones(Pmfc.size)*2
    #wpcc = Pmax * SPcc * np.ones(Pmfc.size)
    #wpcc2 = wpcc * PSF
    # batterypack
    energybatterykwh = RHS * mjtokwh * 1E-3
    Vfccell = 0.65
    Ifccell = 450
    nfccell = 167  # 335
    nseries = 2
    print('check2')
    nparallel = np.ceil((numberFC+1)/nseries)
    print(nparallel)# 2
    Vfc = nfccell * Vfccell
    Vfcstack = Vfc * nseries
    Ifcstack = Ifccell * nparallel
    # batteries
    space = 2
    Vbatt = 3.4
    pmax = (Pmax - Pmfc)*1E3  # 851E3
    rho = 1000
    Ispec = 1400
    minseries = np.ceil(Vfcstack / Vbatt)
    print('check1',minseries)
    Ireq = pmax / (minseries * Vbatt)
    diametro = 21
    storage = energybatterykwh*1E3
    l = 70.5
    V = (diametro / 2) ** 2 * np.pi * l
    Vl = V * 10 ** (-6)
    ec = Vl * rho

    spe = 0.448
    rows = np.ceil(storage / (ec * minseries))
    print('check3',ec)
    width = ((rows + 1) / 2) * diametro + space * (rows / 2 + 1)
    lenght = (minseries / 2) * diametro + space * (minseries / 2 + 1)
    crow = 1400 * (ec / spe / 1000)
    currenttotb = crow * rows
    rows=np.array(rows)
    mcell = Ireq / Ispec * 1 / rows
    print('check4',np.where(wbe1>wbe2)[1] )
    print('check5',wbefinal.shape)
    print('check6', wbe1.shape)
    wbefinal[0,np.where(wbe1>wbe2)[1]]=wbefinal[0,np.where(wbe1>wbe2)[1]]*(2+rows[np.where(wbe1>wbe2)[1]])/rows[np.where(wbe1>wbe2)[1]]
    #wgearbox=252
    Wtot = wpmp + wbefinal + wemp + wpcc #+wgearbox
    #Wtot2 = wpmp2 + wbefinal2 + wemp + wpcc2
    M = np.amin(Wtot)
    I = np.where(Wtot == np.amin(Wtot))[1][0]
    print('energy fuelc cell',energymp[0,I]*mjtokwh*1E-3)
    #pesopsf = Wtot2[0, I]
    POWERRRRR = Pmfc[0, I]
    print('hello',Pmfc[0,I+1])
    energybatterykwh = RHS[I] * mjtokwh * 1E-3
    volumedensity = 620
    volumeb = (RHS[I]) * mjtokwh * 1E-3 / volumedensity
    #fuelweight = 1 / 2.73 * wemp[0, I]
    fuelweight=fuel
    if switch:
        plt.plot(Pmfc.T, wpmp.T, linewidth=5, label='Fuel cell and converter weight')
        plt.plot(Pmfc.T, wemp.T, linewidth=5, label='Tank and hydrogen weight')
        plt.plot(Pmfc.T, wbefinal.T, linewidth=5, label='Weight of batteries')
        plt.plot(Pmfc.T, wpcc.T, linewidth=5, label='E-motors')
        plt.plot(Pmfc.T, wbe1.T)
        plt.plot(Pmfc.T, Wtot.T, linewidth=5, label='Total weight')
        plt.xlabel('Max rated power fuel cell [kW]')
        plt.ylabel('Weight [kg]')
        plt.grid()
        plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc='lower left',
               ncol=3, mode="expand", borderaxespad=0.)
        plt.show()
    if switch:
        plt.plot(Pmfc.T, wbe1.T, linewidth=5, label='Battery mass sized on energy')
        plt.plot(Pmfc.T, wbe2.T,linewidth=5, label='Battery mass sized on power')
        plt.plot(Pmfc.T,wbefinal.T,linewidth=5, label='tot')
        plt.xlabel('Max rated power fuel cell [kW]')
        plt.ylabel('Battery weight [kg]')
        plt.grid()
        plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc='lower left',ncol=3, mode="expand", borderaxespad=0.)
        plt.show()
    if switch:

        plotforme = np.concatenate((ratingFC[I] * np.ones(len(time[1:Toend * 10:dc])),0.7 * ratingFC[I] * np.ones(len(time[Toend * 10:end2 * 10:dc])),ratingFC[I] * np.ones(len(time[end2 * 10::dc]))), axis=None)
        plotforme = np.array([plotforme])
        correct= np.ones(len(teto.T))*ratingFC[I]*0.7
        c2=plotforme-correct
        print(c2.shape)
        print(teto.shape)
        print(plotforme.shape)
        print(c2)

        plt.plot(time, powert)
        plt.plot(teto.T, plotforme.T)
        plt.fill_between(teto.T,plotforme[:,0], powert[::dc], where=powert[::dc]>plotforme[:,0], color='lightskyblue' )

        plt.fill_between(teto.T, plotforme[:, 0], 0,where=c2[0,:]>10,color='peachpuff')
        plt.fill_between(teto.T, correct, 0, where=correct < powert[::dc], color='peachpuff')
        plt.fill_between(teto.T, correct, powert[::dc], where=c2[0, :] < 10, color='lightskyblue')
        plt.ylim((0, 1300))
        plt.show()
    # plt.plot(teto, powertTO[0,:])
    # plt.plot(time, powert)
    # plt.show()


    if switch:
        print('results')
        print("Fuel cell=", name2FC[I])
        print("Fuel cell rating=", ratingFC[I], numberFC[I]+1)
        print("Fuel cell mass=", massFC[I])
        print("Radiator mass=", massRAD[0, I]*2)
        print("Radiator area=", area_tot[0,I])
        print("DCDC converter=", nameDCDC)
        print("DCDC converter rating=", DataDCDC, numberDCDC)
        print("DCDC converter mass=", massDCDC)
        print("Compressor =", name2C[I])
        print("Compressor rating=", ratingC[I], numberC[I]+1)
        print("Compressor mass=", massC[I])
        print("Electric motor=", nameEM2)
        print("Electric motor rating=", ratingEM*2)
        print("Electric motor mass=", massEM*2)
        print("Electric motor controller=", nameMC2)
        print("Electric motor controller rating=", ratingMC*2,(numberMC+1)*2)
        print("Electric motor controller mass=", massMC*2)
        print("BATTERY MASS=", wbefinal[0,I])
        print("Battery storage=", energybatterykwh, (Pmax - Pmfc)[0,I])
        print("Tank weight=",tank[0, I])
        print("Hydrogen fuel weight=", fuel2[0, I])
        print("Tot mass power plant=", M)
        print('Vfcstack', Vfcstack)
        print('Ifcstack', Ifcstack[I])
        print('Ireq', Ireq[0,I])
        print('series', minseries)
        print('volume', Vl,ec,energybatterykwh)
        print('parallel', rows[I]+2)
        print('width', width[I])
        print('lenght', lenght)
        print('mass cell', mcell[0,I])
        print('vbatt', minseries * Vbatt)
        print(Pmax/2)
        #print('current batt for motor', Ireq)


    return M, massFC[I], massDCDC, massC[I], massRAD[0, I]*2, massEM*2, massMC*2, tank[0, I], fuel2[0, I], wbefinal[0,I], m_tankinnerwall[I], m_radshield[I], m_spacer[I], m_tankouterwall[I], m_tanksupport[I]

if __name__ == "__main__":

    N, B, R, cutout, solidity, theta_tip, taper, airfoil, be = rp.initVariables(3.25)
    rotor = Rotor(N, B, R, cutout, taper, solidity, theta_tip, airfoil, 0, be, 0)
    rotor.calcTwist('linear', 10 * math.pi / 180, math.pi / 180)
    rp.calcPowersFast(rotor)

    cruise_altitude = 4000  # m
    aircraft = Mission.TiltRotor(cruise_altitude, 350, 270, 220, 4000., 9.80665, 8, 45, 6.17, 0.9, 60, 1.15, 0.75)
    aircraft.performFlight(1.225)
    data = aircraft.data
    timestamps = aircraft.timestamps
    powert = data['powert']
    energy = data['energy']
    time = data['time']
    TOe = np.around(timestamps['cruise speed'], decimals=1)
    tendcruise = np.around(timestamps['Min Height'], decimals=1)
    Mintot, massFC, massDCDC, massC, massRAD, massEM, massMC, tank, fuel2, wbefinal, m_tankinnerwall, m_radshield, m_spacer, m_tankouterwall, m_tanksupport = calculateWeightFPP(
        powert, time, TOe, tendcruise, True)

# import sympy as sy
# x = sy.symbols('x')
# x=np.linspace(0,1,101)
# #RHSsc = (powertTO - Pmfc2) * (teto[2] - teto[1])
# #Energy=(RHS) * mjtokwh * 1E-3
# Power=(Pmax - Pmfc.T)
# x2, Energy2= np.meshgrid(x, Energy)
# x2, Power2= np.meshgrid(x, Power)
# ##
# f=(x2*0.448+(1-x2)*15)*Energy2
# g=(x2*1.3+(1-x2)*14)*Power2
# final= np.maximum(f, g)
# minimum= np.amin(final, axis=1)
# minimumindex=np.argmin(final, axis=1)
# plt.plot(Pmfc[0, :], minimumindex)
# plt.show()
#
#
# ##
# def scbatt(Pmax, Pmfc, SEsc, SPsc, SEb, SPb, RHS, powertTO, teto):
#     x = np.linspace(0, 1, 101)
#     Power = (Pmax - Pmfc)
#     teto2, pb = np.meshgrid(teto, (Power*x).T)
#     energyb = (powertTO - Pmfc - pb) * (teto[2] - teto[1])
#     energyb = energyb * mjtokwh * 1E-3
#     energyb[energyb <= 0] = 0
#     energyb = np.cumsum(energyb, axis=1)
#     energyb = energyb[:, -1]
#     energysc= RHS-energyb
#     #energyb[energyb < 0] = 1e-10
#     powermass = (x*SPb+(1-x)*SPsc) *Power
#     energymass = energyb * SEb + energysc * SEsc
#     final = np.maximum(powermass, energymass)
#
#     minimum = np.amin(final)
#     minimumindex = np.argmin(final)
#
#     return minimum, np.array(x[minimumindex])
# ##
# minimum= np.ones(len(Pmfc.T))
# xminimum=np.ones(len(Pmfc.T))
# i=0
# for i in range(len(Pmfc.T) - 1):
#     minimum[i], xminimum[i]=scbatt(Pmax, Pmfc[0,i], 0.2, 14, 0.466, 1, RHS[i], powertTO, teto)
# xminimum=np.array(xminimum)
# ##
# plt.plot(Pmfc[0, :], xminimum)
# plt.show()
#