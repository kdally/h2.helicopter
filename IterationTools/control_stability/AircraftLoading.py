import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

def convertToMAC(x,xLEMAC,MAC):
    return (x-xLEMAC)/MAC

def findEmptyCG(MAC,xLEMAC,fus,fixed,horizontal,vertical,nose,prop,nacelle,wing,landing,battery,tank,fuelcell,radiator):
    """
    Given the different components in the aircraft the c.g. of the empty tilt rotor is estimated as a distance from the cg as a fraction of MAC
    input: MAC [mean aerodynamic chord] (float), xLEMAC [position of the Leading Edge of the MAC from nose] (float), fus --> radiator are lists of distance from nose and mass [distance, mass]
    output: the OEW cg (float), the OEW mass (float)
    """

    wing_pos = xLEMAC*1.05
    prop[0] = xLEMAC
    wing[0] = wing_pos
    landing[0] = wing_pos
    radiator[0] = xLEMAC + 0.5*MAC
    nacelle[0] = xLEMAC + 0.5*MAC

    #put lists in list for iteration
    components = [fus,fixed,horizontal,vertical,nose,prop,nacelle,wing,landing,battery,tank,fuelcell,radiator]

    #set counters
    mass = 0
    weighted_distance = 0

    #iterate through components
    for i in components:
        mass += i[1]
        weighted_distance += i[0]*i[1]

    #find cg
    nosecg = weighted_distance/mass
    cg = (nosecg-xLEMAC)/MAC    

    return cg, mass

def findCGRange(OEW,MAC,xLEMAC,emptycg,payload,passengermass,fuel,cargo,firstperson,seatpitch,fuelpos):
    """
    Given the OEW and empty cg find the cg range of loading the aircraft (baggage ---> passengers ----> fuel)
    Inputs: OEW mass (float), OEW cg (float), ammount of payload (float), mass of a passanger (float), fuel mass (float), distance of first person from lemac in terms of MAC (float), 
    seat pitch in terms of MAC (float), fuel position from xLEMAC in terms of MAC (float)
    Output: list of cg range [lowercg, uppercg]
    """
    #convert values to in terms of LEMAC
    cargo = convertToMAC(cargo,xLEMAC,MAC)
    firstperson = convertToMAC(firstperson,xLEMAC,MAC)
    seatpitch = seatpitch/MAC
    fuelpos = convertToMAC(fuelpos,xLEMAC,MAC)


    #set up lists for different loading cases for potato diagram
    rightstreamx = []
    rightstreamy = []
    leftstreamx = []
    leftstreamy = []

    #Baggage
    #find end point and create line to it
    baggagemass = payload - 6 * passengermass
    mass = OEW+baggagemass
    baggagecg = ((OEW*emptycg)+(cargo*baggagemass))/mass

    for i in range(5):
        rightstreamx.append(emptycg+(i+1)*(baggagecg-emptycg)/5)
        rightstreamy.append(OEW + (i+1)*(mass-OEW)/5)

    #fronttoback passenger loading
    tempmass = mass
    passengercg1 = baggagecg
    for i in range(3):
        passengercg1 = (passengercg1*tempmass+(firstperson+i*seatpitch)*passengermass*2)/(tempmass+passengermass*2)
        tempmass += passengermass*2
        rightstreamx.append(passengercg1)
        rightstreamy.append(tempmass)

    #backtofront passenger loading
    passengercg2 = baggagecg
    leftstreamx.append(passengercg2)
    leftstreamy.append(mass)
    for i in range(3):
        passengercg2 = (passengercg2*mass+(firstperson+(2-i)*seatpitch)*passengermass*2)/(mass+passengermass*2)
        mass += passengermass*2
        leftstreamx.append(passengercg2)
        leftstreamy.append(mass)

    #fuel loading with straight line
    fuelcg = (mass*passengercg2+fuel*fuelpos)/(mass+fuel)

    for i in range(5):
        rightstreamx.append(passengercg2 + (i+1)*(fuelcg-passengercg2)/5)
        rightstreamy.append(mass + (i+1)*(fuel)/5)

    # leftstreamxplot = leftstreamx.copy()
    # leftstreamxplot[1] = 0.499
    # rightstreamxplot = rightstreamx.copy()
    # rightstreamxplot[-1] = 0.498
    # idx = rightstreamxplot.index(min(rightstreamxplot))
    # rightstreamxplot[idx] = 0.443
    # plt.plot(rightstreamxplot,rightstreamy)
    # plt.plot(leftstreamxplot,leftstreamy)
    # plt.xlabel('xcg/MAC [-]', size =14)
    # plt.ylabel('Mass [kg]', size = 14)
    # plt.xticks(fontsize=12)
    # plt.yticks(fontsize=12)
    # plt.show()
    MTOW = rightstreamy[-1]

    return [min(rightstreamx),max(max(rightstreamx),max(leftstreamx))], MTOW



def constructLoadingDiagram(MAC,fuselage,fus,fixed,horizontal,vertical,nose,prop,nacelle,wing,landing,battery,tank,fuelcell,radiator,payload,passengermass,fuel,cargo,firstperson,seatpitch,fuelpos,plot):
    """
    Generate cg range curves over a range of xLEMAC positions
    Inputs: an approximate LEMAC value (float), MAC (float), fuselage length (flaot), lists for distance from nose and mass {fus ---> radiator} [distance from nose, mass], OEW mass (float), payload mass (float), 
    mass of a passenger (float),  fuel mass (float), distance of first person from lemac in terms of MAC (float), seat pitch in terms of MAC (float), fuel position from xLEMAC in terms of MAC (float), 'yes' if plotting (string)
    Outputs: slope and intercepts for aft and forward cg curves (floats)
    """

    #Create xLEMAC range for cg range calculation
    xLEMAC = np.linspace(0.3,0.8,20)*fuselage
    xLEMACplot = xLEMAC*(1/fuselage)
    #ists for cg range values
    lowercgs = []
    uppercgs = []

    #iterate through different xLEMACs and append c.g. range for each to lists
    for i in xLEMAC:
        emptycg, OEW = findEmptyCG(MAC,i,fus,fixed,horizontal,vertical,nose,prop,nacelle,wing,landing,battery,tank,fuelcell,radiator)
        cgrange, MTOW = findCGRange(OEW,MAC,i,emptycg,payload,passengermass,fuel,cargo,firstperson,seatpitch,fuelpos)
        lowercgs.append(cgrange[0])
        uppercgs.append(cgrange[1])

    #plotting
    if plot == 'yes':     
        plt.plot(lowercgs,list(xLEMACplot), label = 'foward')
        plt.plot(uppercgs,list(xLEMACplot), label = 'aft')
        plt.ylabel('xLEMAC/Fuselage Length')
        plt.xlabel('xcg/MAC')
        plt.legend()
        plt.show()   

    #use linear regression to find slope and intercepts of each curve
    m_f, b_f, _,_,_ = stats.linregress(lowercgs,list(xLEMACplot))
    m_a, b_a, _,_,_ = stats.linregress(uppercgs,list(xLEMACplot))

    return m_f, b_f, m_a, b_a, OEW, MTOW


def main():
    pass

if __name__ == "__main__":
    # cgrange = findCGRange(2861,0.35,900,90,20,1.8,-2.88,0.22,2.5)
    """
    #cg esimation


    fuselage = 11.6
    MAC = 2
    xLEMAC = 0.5*fuselage
    wing_pos = xLEMAC*1.05
    fus = [fuselage*0.42, 277.6]
    fixed = [fuselage*0.38, 625.8]
    pilot = [fuselage*0.17, 100]
    horizontal = [fuselage*0.2, 72]
    vertical = [fuselage, 72]
    nose = [fuselage*0.13, 60]
    prop = [wing_pos, 300]
    wing = [wing_pos, 318]
    landing = [wing_pos*1.1, 120]
    battery = [wing_pos, 360]
    tank = [fuselage*0.95, 25]
    fuelcell = [fuselage*0.8,180]
    radiator = [wing_pos, 90]

    emptycg, OEW = findEmptyCG(MAC,xLEMAC,fus,fixed,pilot,horizontal,vertical,nose,prop,wing,landing,battery,tank,fuelcell,radiator)
    print('emptycg',emptycg)
    print('OEW', OEW)

    fuelpos = convertToMAC(fuelcell[0],xLEMAC,MAC)
    firstperson = convertToMAC(fuselage*0.2,xLEMAC,MAC)
    seatpitch = convertToMAC(1,xLEMAC,MAC)
    cargo = convertToMAC((fuselage*0.2+3)*1.03,xLEMAC,MAC)
    #findCGRange(OEW,emptycg,payload,passengermass,fuel,cargo,firstperson,seatpitch,fuelpos)
    # cg = findCGRange(OEW,emptycg,900.,90.,24.,cargo,firstperson,seatpitch,fuelpos)
    # print('cgrange', cg)

    m_f, b_f, m_a, b_a = constructLoadingDiagram(xLEMAC,MAC,fuselage,fus,fixed,pilot,horizontal,vertical,nose,prop,wing,landing,battery,tank,fuelcell,radiator,OEW,900.,90.,24.,cargo,firstperson,seatpitch,fuelpos,'yes')

    """