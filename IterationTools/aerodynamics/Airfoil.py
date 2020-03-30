import matplotlib.pyplot as plt
import numpy as np


# Reynolds Number
def ReynoldNumber(c, Vcruise, rho, mu):
    Re = (c*Vcruise*rho)/(mu)
    return Re


# Mach Number at Cruise
def Mcruise_calc(Vcruise, Tcruise, gamma, R):
    Mcruise = Vcruise/np.sqrt(Tcruise*gamma*R)
    return Mcruise


# Dynamic Pressure at Cruise
def Dynam_press(rho, Vcruise):
    q = 0.5*rho*Vcruise**2
    return q


# Thickness over Chord
def t_over_c(Mcruise, MTO, q, WArea, g):
    Mdd = Mcruise + 0.03
    tc = (0.935 - Mdd)-0.115*(MTO*g/(q*WArea))**1.5
    tc = np.minimum(tc, 0.18)
    return tc


# Lift of the Wing
def Lwing_calc(MTO, g):
    Lwing = 1.1*MTO*g
    return Lwing


# Desired Lift Coefficient of the Wing
def CLdes_calc(q, WoverS):
    CLdes = 1.1*(1/q)*(1/2)*(WoverS+1853.646839)
    return CLdes


# Desired Lift Coefficient of the Airfoil
def Cldes_calc(Lwing, q, WArea):
    Cldes = Lwing/(q*WArea)
    return Cldes


# Airfoil and Wing Dimensions
def AreaAirfoil(k, tc, c):
    AA = ((k+3)/6)*tc*c**2
    return AA


def ThicknessAirfoil(tc, c):
    tairf = tc * (c)
    return tairf


def SpanWing(WArea, c):
    b = WArea/(c)
    return b


def WingVolume(AA, b):
    VV = AA * (2*b-1.48)
    return VV


def main():
    pass


if __name__ == "__main__":

    # Variables
    Vcruise = 97.222
    Tcruise = 275.5
    gamma = 1.4
    R = 287.15
    rho0 = 1.225
    MTO = 4000
    g = 9.806
    WArea = 21.034880850919492
    WoverS = 1864.71225
    rho = 1.007
    mu = 1.726 * 10 ** (-5)
    c = 2  # m
    k = c*0.3


    # Results
    print("Reynolds Number (Re): " + str(ReynoldNumber(c, Vcruise, rho, mu)) + " -")
    print("Mach Number at Cruise (M): " + str(Mcruise_calc(Vcruise, Tcruise, gamma, R)) + " -")
    print("Dynamic Pressure (q): " + str(Dynam_press(rho, Vcruise)) + " Pa")
    print("Thickness over Chord (t/c): " + str(t_over_c(Mcruise_calc(Vcruise, Tcruise, gamma, R), MTO, Dynam_press(rho, Vcruise), WArea, g)) + " -")
    print("Wing Lift (L): " + str(Lwing_calc(MTO, g)) + " N")
    print("Design Lift Coefficient Wing (CL): " + str(CLdes_calc(Dynam_press(rho, Vcruise), WoverS)) + " -")
    print("Design Lift Coefficient Airfoil (Cl): " + str(Cldes_calc(Lwing_calc(MTO, g), Dynam_press(rho, Vcruise), WArea)) + " -")

    # Airfoil Plot
    xcoord = [0, 0.0125, 0.025, 0.05, 0.075, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 1]
    new_xcoord = [i * c for i in xcoord]
    ycoordup = [0, 0.0409, 0.0529, 0.0692, 0.0801, 0.0883, 0.0986, 0.1036, 0.1056, 0.1055, 0.1004, 0.0905, 0.0775, 0.0618, 0.044, 0.0239, 0.0132, 0]
    new_ycoordup = [i * c for i in ycoordup]
    ycoorddown = [0, -0.0183, -0.0271, -0.038, -0.046, -0.0522, -0.0618, -0.0686, -0.0727, -0.0747, -0.0737, -0.0681, -0.0594, -0.0482, -0.0348, -0.0194, -0.0109, 0]
    new_ycoorddown = [i * c for i in ycoorddown]
    chordlinex = [0, c]
    chordliney = [0, 0]
    camber = [0, 0.0113, 0.0129, 0.0156, 0.01705, 0.01805, 0.0184, 0.0175, 0.01645, 0.0154, 0.01335, 0.0112, 0.00905, 0.0068, 0.0046, 0.00225, 0.00115, 0]
    new_camber = [i * c for i in camber]

    plt.plot(new_xcoord, new_ycoordup, c='r', label='Airfoil Surface')
    plt.plot(new_xcoord, new_ycoorddown, c='r')
    plt.plot(chordlinex, chordliney, label='Chord Line', c='m')
    plt.plot(new_xcoord, new_camber, label='Camber Line', c='g')

    plt.title('Airfoil NACA 23018')
    plt.xlabel('Length [m]')
    plt.ylabel('Height [m]')
    plt.grid()
    plt.axis('equal')
    plt.legend(loc='upper left')
    #plt.savefig('airfoil.png')
    plt.show()

    # Cl vs alpha & Cd vs alpha plot
    alpha = [-19.75, -19.5, -19.25, -19, -18.75, -18.5, -18.25, -18, -17.75, -17.5, -17.25, -17, -16.75, -16.5, -16.25, -16, -15.75, -15.5, -15.25, -15, -14.75, -14.5, -14.25, -14, -13.75, -13.5, -13.25, -13, -12.75, -12.5, -12.25, -12, -11.75, -11.5, -11.25, -11, -10.75, -10.5, -10.25, -10, -9.75, -9.5, -9.25, -9, -8.75, -8.5, -8.25, -8, -7.75, -7.5, -7.25, -7, -6.75, -6.5, -6.25, -6, -5.75, -5.5, -5.25, -5, -4.75, -4.5, -4.25, -4, -3.75, -3.5, -3.25, -3, -2.75, -2.5, -2.25, -2, -1.75, -1.5, -1.25, -1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 2, 2.25, 2.5, 2.75, 3, 3.25, 3.5, 3.75, 4, 4.25, 4.5, 4.75, 5, 5.25, 5.5, 5.75, 6, 6.25, 6.5, 6.75, 7, 7.25, 7.5, 7.75, 8, 8.25, 8.5, 8.75, 9, 9.25, 9.5, 9.75, 10, 10.25, 10.5, 10.75, 11, 11.25, 11.5, 11.75, 12, 12.25, 12.5, 12.75, 13, 13.25, 13.5, 13.75, 14, 14.25, 14.5, 14.75, 15, 15.25, 15.5, 15.75, 16, 16.25, 16.5, 16.75, 17, 17.25, 17.5, 17.75, 18, 18.25, 18.5, 18.75, 19, 19.25]
    Cl = [-1.3032, -1.3096, -1.3135, -1.3154, -1.3155, -1.3141, -1.3119, -1.3084, -1.3044, -1.2988, -1.2917, -1.294, -1.2923, -1.2869, -1.2789, -1.2697, -1.2593, -1.248, -1.2357, -1.2221, -1.2131, -1.2064, -1.1951, -1.1819, -1.1673, -1.1523, -1.1371, -1.1213, -1.1088, -1.102, -1.0939, -1.0885, -1.0905, -1.095, -1.0649, -1.0328, -1.0014, -0.9684, -0.9374, -0.9049, -0.8712, -0.8413, -0.8154, -0.7877, -0.7619, -0.7397, -0.7168, -0.6929, -0.6694, -0.6459, -0.6211, -0.5972, -0.5727, -0.5483, -0.5242, -0.4997, -0.4761, -0.4539, -0.4317, -0.4076, -0.3828, -0.3573, -0.3312, -0.3048, -0.2783, -0.2514, -0.2244, -0.1974, -0.1704, -0.1436, -0.1167, -0.0896, -0.0625, -0.0356, -0.0082, 0.0191, 0.0461, 0.0736, 0.1013, 0.1285, 0.1559, 0.1832, 0.2103, 0.2377, 0.2646, 0.2915, 0.3456, 0.3724, 0.3989, 0.4259, 0.4527, 0.4791, 0.5062, 0.5322, 0.5591, 0.5848, 0.6118, 0.6377, 0.6634, 0.6895, 0.715, 0.7395, 0.7657, 0.7921, 0.8174, 0.8466, 0.8767, 0.9073, 0.9421, 0.9773, 1.0131, 1.047, 1.0775, 1.1125, 1.1472, 1.1803, 1.2143, 1.2464, 1.2687, 1.2758, 1.2775, 1.279, 1.2852, 1.2951, 1.3098, 1.323, 1.339, 1.3535, 1.3678, 1.3832, 1.3967, 1.4093, 1.4245, 1.4377, 1.449, 1.4586, 1.4724, 1.484, 1.4934, 1.4999, 1.5079, 1.5181, 1.5253, 1.5291, 1.529, 1.5368, 1.5402, 1.539, 1.5348, 1.5359, 1.5317, 1.5215, 1.5126, 1.5036, 1.489, 1.4698]
    Cd = [0.07356, 0.06953, 0.06583, 0.0624, 0.05925, 0.05632, 0.05353, 0.05092, 0.04842, 0.04613, 0.04404, 0.04111, 0.03864, 0.0366, 0.03482, 0.03321, 0.03175, 0.03042, 0.0292, 0.02811, 0.02674, 0.02528, 0.02419, 0.02329, 0.02251, 0.0218, 0.02115, 0.02057, 0.01989, 0.01907, 0.01854, 0.01815, 0.01789, 0.01769, 0.01725, 0.01685, 0.01617, 0.01564, 0.01522, 0.01483, 0.0145, 0.0139, 0.01354, 0.01321, 0.01294, 0.01255, 0.01224, 0.012, 0.01173, 0.01144, 0.01122, 0.01094, 0.01067, 0.0104, 0.01008, 0.0098, 0.00946, 0.00899, 0.00853, 0.00824, 0.00802, 0.00787, 0.00777, 0.00768, 0.00759, 0.00752, 0.00748, 0.00744, 0.00742, 0.00737, 0.00733, 0.00731, 0.00729, 0.00729, 0.00728, 0.00727, 0.00727, 0.00727, 0.00729, 0.00733, 0.00733, 0.00736, 0.00742, 0.00747, 0.00751, 0.00758, 0.00776, 0.00785, 0.00797, 0.00809, 0.00819, 0.00834, 0.00845, 0.00862, 0.00874, 0.00893, 0.00903, 0.0092, 0.00938, 0.00951, 0.00969, 0.00992, 0.01004, 0.01022, 0.01043, 0.01071, 0.01087, 0.01108, 0.01135, 0.01166, 0.01185, 0.0121, 0.01245, 0.01265, 0.01291, 0.01327, 0.01351, 0.01389, 0.01414, 0.01436, 0.01462, 0.01478, 0.01506, 0.01547, 0.01586, 0.0164, 0.01688, 0.01747, 0.01811, 0.01873, 0.01947, 0.02031, 0.02103, 0.02189, 0.02291, 0.0241, 0.02503, 0.02616, 0.0275, 0.02913, 0.03068, 0.0321, 0.03383, 0.03592, 0.03845, 0.04029, 0.04262, 0.04549, 0.04875, 0.05151, 0.05493, 0.05909, 0.0632, 0.0674, 0.07235, 0.07799]

    plt.plot(alpha, Cl)
    plt.title('Lift Coefficient (Cl) vs. Angle of Attack')
    plt.xlabel('Angle of Attack [°]')
    plt.ylabel('Lift Coefficient (Cl) [-]')
    plt.grid()
    #plt.savefig('Clvsalpha.png')
    plt.show()

    plt.plot(alpha, Cd)
    plt.title('Drag Coefficient (Cd) vs. Angle of Attack')
    plt.xlabel('Angle of Attack [°]')
    plt.ylabel('Drag Coefficient (Cd) [-]')
    plt.grid()
    #plt.savefig('Cdvsalpha.png')
    plt.show()


    print("Maximum Thickness Airfoil: " + str(ThicknessAirfoil(t_over_c(Mcruise_calc(Vcruise, Tcruise, gamma, R), MTO, Dynam_press(rho, Vcruise), WArea, g), c)) + " m")

    print("Area Airfoil: " + str(AreaAirfoil(k, t_over_c(Mcruise_calc(Vcruise, Tcruise, gamma, R), MTO, Dynam_press(rho, Vcruise), WArea, g), c)) + " m^2")

    print("Wing Span: " + str(SpanWing(WArea, c)) + " m")

    print("Wing Volume: " + str(WingVolume(AreaAirfoil(k, t_over_c(Mcruise_calc(Vcruise, Tcruise, gamma, R), MTO, Dynam_press(rho, Vcruise), WArea, g), c), SpanWing(WArea, c))) + " m^3")