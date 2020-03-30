

def main():
    pass   
if __name__ == "__main__":
    from AircraftLoading import findEmptyCG, findCGRange, constructLoadingDiagram, convertToMAC
    from LoadingScissorOptimisation import checkSwitch, findScissorOptimum, findScissorCG, plotScissorOptimum, plotCustomCG
    import Scissor as ScissorCode
    import matplotlib.pyplot as plt
    from math import sqrt,pi,tan
    import numpy as np
else:
    from control_stability.AircraftLoading import findEmptyCG, findCGRange, constructLoadingDiagram, convertToMAC
    from control_stability.LoadingScissorOptimisation import checkSwitch, findScissorOptimum, findScissorCG, plotScissorOptimum, plotCustomCG
    import control_stability.Scissor as ScissorCode
    import matplotlib.pyplot as plt
    from math import sqrt,pi,tan
    import numpy as np


def runControlStabilitySizing(MAC,fuselage,fus,fixed,horizontal,vertical,nose,prop,nacelle,wing,landing,battery,tank,fuelcell,radiator,payload,passengermass,fuel,cargo,firstperson,seatpitch,fuelpos,c,Ah,Lambda_ch,Cl_alphaw,Sw,Cm_acf,Cm_acw,lh,b,Cl_f,Cl_alphaf,xac_w,SM,Cl_w,Sf,lf):
    """
    This function unifies all the control and stability functions into one function to run on the Main.py file.
    Inputs: OEW (float), fuel weight (float), fuel position (float), cargo position (float), location of first person (float), seat pitch (float), 
    Aircraft components distance from nose and mass in lists [distance, mass] {fus,horizontal,vertical,nose,prop,wing,landing},  MAC (float), xLEMAC approximation (float), 
    aerodynamics characteristics (float) {}.
    Ouputs: cg range [lowercg, uppercg], Sh/S (float), xLEMAC (float), tail arm (float), control surface sizes (floats)

    Written by: Benjamin
    """
    #set up tail arm iteration values
    # lh_list = np.linspace(5,5.1,1)
    # cg_aft = []
    # cg_forward = []
    # Sh_S_list = []
    # xLEMAC_list = []



    # #iterate on xLEMAC and lh
    # for i in lh_list:

    ##construct loading curves 
    m_f, b_f, m_a, b_a, OEW, MTOW = constructLoadingDiagram(MAC,fuselage,fus,fixed,horizontal,vertical,nose,prop,nacelle,wing,landing,battery,tank,fuelcell,radiator,payload,passengermass,fuel,cargo,firstperson,seatpitch,fuelpos,'no')
    ##construct scissor diagram
    m, q = ScissorCode.Scissor(MAC,Ah,Lambda_ch,Cl_alphaw,Sw,Cm_acf,Cm_acw,lh,b,Cl_f,Cl_alphaf,xac_w,SM,Cl_w,Sf,lf)

    ##find optimum
    m = [m[0], m[1], m_a, m_f]
    b = [q[0], q[1], b_a, b_f]
    # cg, Sh_S, xLEMAC, c = findScissorCG(m, b, 0.49)
    # plotCustomCG(m,b,cg,Sh_S,c)

    cg, Sh_S, xLEMAC, scale = findScissorOptimum(m, b) 

    if abs(0.5 - cg[1]) > 0.001:
        cg, Sh_S, xLEMAC, c = findScissorCG(m,b,0.4999)
        # plotCustomCG(m,b,cg,Sh_S,c)
    # else:
        # plotScissorOptimum(m,b,scale,cg,Sh_S)

        
    # cg_forward.append(cg[0])
    # cg_aft.append(cg[1])
    # Sh_S_list.append(Sh_S)
    # xLEMAC_list.append(xLEMAC)

    ##find control surface areas
    # plt.plot(lh_list,Sh_S_list)
    # plt.xlabel('lh [m]')
    # plt.ylabel('Sh/S [-]')
    # plt.show()
    # tail_arm = sum(lh_list)/len(lhl_ist)
    tail_arm = lh

    return cg, Sh_S, xLEMAC, tail_arm, OEW, MTOW#, aileron, elevator, rudder




if  __name__ == "__main__":
    fuselage = 9.73
    MAC = 2
    xLEMAC = 0.5*fuselage
    wing_pos = xLEMAC*1.05
    fus = [fuselage*0.5, 277.6]
    fixed = [fuselage*0.38, 625.8]
    nacelle = [fuselage*0.17, 100]
    horizontal = [fuselage*0.2, 72]
    vertical = [fuselage, 72]
    nose = [fuselage*0.13, 80]
    prop = [wing_pos, 400]
    wing = [wing_pos, 318]
    landing = [wing_pos*1.1, 200]
    battery = [fuselage*0.7, 360]
    tank = [fuselage*0.95, 25]
    fuelcell = [fuselage*0.8,180]
    radiator = [wing_pos, 90]
    fuelpos = fuelcell[0]
    firstperson = fuselage*0.25
    seatpitch = 1
    cargo = (fuselage*0.2+3)*1.03

    c=3
    Ah=30
    Lambda_ch=20*pi/180
    Cl_alphaw=6.0
    Cm_acw=-0.05 # moment coefficient aerodynamic center 
    Cm_acf=-0.15# moment coefficient aerodynamic center
    lh=5.8
    lf=-0.3
    b=14
    xac_w= 0.25 # position of the ac of the main wing as percentage of the MAC
    SM=0.04 # static margin
    Cl_w=0.4
    Cl_f=0.2
    Cl_alphaf=3
    Sw=20
    Sf=32

    emptycg, OEW = findEmptyCG(MAC,xLEMAC,fus,fixed,horizontal,vertical,nose,prop,nacelle,wing,landing,battery,tank,fuelcell,radiator)
    cg, Sh_S, xLEMAC, tail_arm, OEW, MTOW = runControlStabilitySizing(MAC,fuselage,fus,fixed,horizontal,vertical,nose,prop,nacelle,wing,landing,battery,tank,fuelcell,radiator,900,90,24,cargo,firstperson,seatpitch,fuelpos,c,Ah,Lambda_ch,Cl_alphaw,Sw,Cm_acf,Cm_acw,lh,b,Cl_f,Cl_alphaf,xac_w,SM,Cl_w,Sf,lf)
    print('OEW',OEW)
    print('cg',cg)
    print('Sh_S',Sh_S)
    print('xLEMAC', xLEMAC)
    print('tail arm', tail_arm)

    emptycg, OEW = findEmptyCG(MAC,xLEMAC,fus,fixed,horizontal,vertical,nose,prop,nacelle,wing,landing,battery,tank,fuelcell,radiator)
    m_f, b_f, m_a, b_a, OEW, MTOW = constructLoadingDiagram(MAC,fuselage,fus,fixed,horizontal,vertical,nose,prop,nacelle,wing,landing,battery,tank,fuelcell,radiator,900,90,24,cargo,firstperson,seatpitch,fuelpos,'no')
    ##construct scissor diagram
    
    m, q = ScissorCode.Scissor(c,Ah,Lambda_ch,Cl_alphaw,Sw,Cm_acf,Cm_acw,lh,b,Cl_f,Cl_alphaf,xac_w,SM,Cl_w,Sf,lf)
    mi = [m[0], m[1], m_a, m_f]
    bi = [q[0], q[1], b_a, b_f]
    cg, Sh_S, xLEMAC, c = findScissorCG(mi, bi, 0.5) 
    print('OEW',OEW)
    print('cg',cg)
    print('Sh_S',Sh_S)
    print('xLEMAC', xLEMAC)
    print('tail arm', tail_arm)
    plotCustomCG(mi,bi,cg,Sh_S,c)
    