import matplotlib.pyplot as plt
import numpy as np


# Wing Maximum Lift Coefficient
def CLmaxwing_calc(CLmaxoverClmax):
    Clmax = 1.646  # Airfoil
    CLmax = CLmaxoverClmax*Clmax
    return CLmax


# Mach Number at Cruise
def Mcruise_calc(Vcruise, Tcruise, gamma, R):
    Mcruise = Vcruise/np.sqrt(Tcruise*gamma*R)
    return Mcruise


# Oswald Efficiency Factor
def OswaldEfficiencyFActor(A):
    e = 1/(1.05+0.007*np.pi*A)
    return e


# Zero Lift Drag Coefficient
def ZeroLiftDragCoeff():
    CD0 = (0.007*0.15)+0.007
    return CD0


# Graph Line for Stall Speed (no flaps)
def SizingforStallSpeed(rho, Vs, CLmax):
    WS1 = 0.5*rho*Vs**2*CLmax
    return WS1


# Graph Line for Landing
def SizingforLanding(rho0, Vs, CLmax_land, WLWTO):
    WS2 = (CLmax_land*rho0*Vs**2)/(2*WLWTO)
    return WS2


# Graph Line for Cruise Speed
def SizingforCruiseSpeed(rho, CD0, Vcruise, A, e, n_p, x0, x1):
    x = np.arange(x0, x1)
    WP3 = n_p*(rho/rho0)**(3/4)*((CD0*0.5*rho*Vcruise**3)/x+x/(np.pi*A*e*0.5*rho*Vcruise))**(-1)
    return WP3


# Graph Line for Climb Rate
def SizingforClimbRate(rho0, n_p, c, A, e, CD0, x0, x1):
    x = np.arange(x0, x1)
    WP4 = n_p/(c+(np.sqrt(x)*np.sqrt(2/rho0))/((1.345*(A*e)**(3/4))/(CD0**(1/4))))
    return WP4


# Graph Line for Climb Gradient
def SizingforClimbGradient(c_V, CD0, CLclimb, rho0, n_p, A, e, x0, x1):
    x = np.arange(x0, x1)
    CDclimb = CD0+CLclimb**2/(np.pi*A*e)
    WP5 = n_p/(np.sqrt(x)*(c_V+CDclimb/CLclimb)*np.sqrt(2/(rho0*CLclimb)))
    return WP5


# Graph Line for Manoeuvring
def SizingforManoeuvring(nmax, Vcruise, A, e, rho, CD0, n_p, x0, x1):
    x = np.arange(x0, x1)
    WP6 = (2*A*e*n_p*rho*np.pi*Vcruise*x)/(A*e*CD0*rho**2*np.pi*Vcruise**4+4*nmax**2*x**2)
    return WP6


# Calculations for the Design Point
def ypoint1(rho, Vs, CLmax, rho0, n_p, c, A, e, CD0):
    xx = SizingforStallSpeed(rho, Vs, CLmax)
    y1 = n_p/(c+(np.sqrt(xx)*np.sqrt(2/rho0))/((1.345*(A*e)**(3/4))/(CD0**(1/4))))
    return y1


def ypoint2(Vs, CLmax, nmax, Vcruise, A, e, rho, CD0, n_p):
    xx = SizingforStallSpeed(rho, Vs, CLmax)
    y2 = (2*A*e*n_p*rho*np.pi*Vcruise*xx)/(A*e*CD0*rho**2*np.pi*Vcruise**4+4*nmax**2*xx**2)
    return y2


def dp_calc(rho, Vs, CLmax, rho0, n_p, c, A, nmax, Vcruise):
    xx = SizingforStallSpeed(rho, Vs, CLmax)
    dp = [xx, np.minimum(ypoint1(rho0, Vs, CLmax, rho0, n_p, c, A, OswaldEfficiencyFActor(A), ZeroLiftDragCoeff()), ypoint2(Vs, CLmax, nmax, Vcruise, A, OswaldEfficiencyFActor(A), rho0, ZeroLiftDragCoeff(), n_p))]
    return dp


# Optimum Wing Area
def OptimumWingArea(MTO, g, rho0, n_p, c, A, nmax, Vcruise, rho, Vs, CLmax):
    S = (MTO * g) / dp_calc(rho, Vs, CLmax, rho0, n_p, c, A, nmax, Vcruise)[0]
    return S


# Power at Climb
def PowerClimb(MTO, g, rho0, n_p, c, A, nmax, Vcruise, rho, Vs, CLmax):
    ClimbPower = (MTO * g) / dp_calc(rho, Vs, CLmax, rho0, n_p, c, A, nmax, Vcruise)[1]
    return ClimbPower


# Calculation of the Wing Lift Coefficient (DATCOM method)
def WingLiftCurveSlope(A, Mcruise, eta):
    beta = np.sqrt(1-Mcruise**2)
    CLalpha = (2*np.pi*A)/(2+np.sqrt(4+(A*beta)/eta))*(np.pi/180)
    return CLalpha


def CL_calc(CLalpha, alpha0, alpha1, alpha0L):
    alpha = np.arange(alpha0, alpha1)
    CL = CLalpha * (alpha - alpha0L)
    return CL


# Calculation of the Wing Drag Coefficient
def CD_calc(CD0, A, e, Mcruise, eta, alpha0, alpha1, alpha0L):
    CL = CL_calc(WingLiftCurveSlope(A, Mcruise, eta), alpha0, alpha1, alpha0L)
    CD = CD0+CL**2/(np.pi*A*e)
    return CD


# Trim Angle
def TrimAngleCalc(CLdes, CLalpha, alpha0L):
    TrimAngle = CLdes/CLalpha+alpha0L
    return TrimAngle


# Stall Angle
def StallAngleWing_calc(CLmax, CLalpha, alpha0L, DeltaalphaCLmax):
    StallAngle = CLmax/CLalpha + alpha0L + DeltaalphaCLmax
    return StallAngle


# FLAP
def DeltaCLmax_calc(CLmax_land, CLmax):
    DeltaCLmax = CLmax_land - (CLmax + 0.08704651108)
    return DeltaCLmax


def SwfoverS_calc(DeltaCLmax, DeltaClmax):
    SwfoverS = DeltaCLmax/(0.9*DeltaClmax)
    return SwfoverS


def cflap_calc(cfoverc, chord):
    cflap = cfoverc * chord
    return cflap


def lflap_calc(SwfoverS, S):
    lflap = ((SwfoverS * S) / 2)/2
    return lflap


def alphastall_flap_calc(StallAngle, SwfoverS):
    alphastallflap = StallAngle - 15 * SwfoverS
    return alphastallflap


def dragcoeff_flap(cflap, chord, SwfoverS):
    DeltaCdflap = 0.0144 * (cflap/chord) * SwfoverS * (60 - 10)
    return DeltaCdflap


def main():
    pass


if __name__ == "__main__":


    # Variables
    MTO = 4000
    FW = 24
    rho0 = 1.225
    Vs = 50
    CLmax_land = 2
    WLWTO = (MTO-FW)/MTO
    CD0 = 0.017
    Vcruise = 97.222
    A = 5.258
    n_p = 0.75
    c = 8
    c_V = 0.132
    CLclimb = 1.25
    nmax = 3.5
    x0 = 0.1
    x1 = 3500
    g = 9.806
    rho = 1.007
    eta = 0.95
    Tcruise = 275.5
    gamma = 1.4
    R = 287.058
    CLdes = 0.4297198751607246
    alpha0 = -19
    alpha1 = 20
    chord = 2
    alpha0L = -1.6
    CLmaxoverClmax = 0.9
    DeltaalphaCLmax = 2.2
    DeltaClmax = 0.9  # plain flap
    cfoverc = 0.25  # plain flap

    # Generate the line
    x = np.arange(x0, x1)
    line1 = SizingforCruiseSpeed(rho, ZeroLiftDragCoeff(), Vcruise, A, OswaldEfficiencyFActor(A), n_p, x0, x1)
    line2 = SizingforClimbRate(rho0, n_p, c, A, OswaldEfficiencyFActor(A), ZeroLiftDragCoeff(), x0, x1)
    line3 = SizingforClimbGradient(c_V, ZeroLiftDragCoeff(), CLclimb, rho0, n_p, A, OswaldEfficiencyFActor(A), x0, x1)
    line4 = SizingforManoeuvring(nmax, Vcruise, A, OswaldEfficiencyFActor(A), rho, ZeroLiftDragCoeff(), n_p, x0, x1)
    line5 = SizingforStallSpeed(rho, Vs, CLmaxwing_calc(CLmaxoverClmax))
    line6 = SizingforLanding(rho0, Vs, CLmax_land, WLWTO)

    # Plots
    plt.plot(x, line1, 'b', label='Cruise Speed')
    plt.plot(x, line2, 'g', label='Climb Rate')
    plt.plot(x, line3, 'r', label='Climb Gradient')
    plt.plot(x, line4, 'c', label='Manoeuvring Perf.')

    plt.axvline(x=line5, c='m', label='Stall Speed')
    plt.axvline(x=line6, c='y', label='Landing')
    plt.xlim([0, 3300])
    plt.ylim([0, 0.11])
    plt.title('Wing Loading vs. Power Loading')
    plt.xlabel('Wing Loading (W/S) [N/m^2]')
    plt.ylabel('Power Loading (W/P) [N/W]')
    plt.grid()
    plt.legend(loc='upper left')
    #plt.savefig('WingLoadingvsPowerLoading.png')
    plt.show()

    print("Design Point: " + "(" + str(dp_calc(rho, Vs, CLmaxwing_calc(CLmaxoverClmax), rho0, n_p, c, A, nmax, Vcruise)[0]) + ", " + str(dp_calc(rho, Vs, CLmaxwing_calc(CLmaxoverClmax), rho0, n_p, c, A, nmax, Vcruise)[1]) + ")")


    print("Optimum Wing Area: " + str(OptimumWingArea(MTO, g, rho0, n_p, c, A, nmax, Vcruise, rho, Vs, CLmaxwing_calc(CLmaxoverClmax))) + " m^2")


    print("Power at Climb (Pclimb): " + str(PowerClimb(MTO, g, rho0, n_p, c, A, nmax, Vcruise, rho, Vs, CLmaxwing_calc(CLmaxoverClmax))) + " W")


    print("Wing Lift Curve Slope: " + str(WingLiftCurveSlope(A, Mcruise_calc(Vcruise, Tcruise, gamma, R), eta)) + " deg")


    CLalphaline = CL_calc(WingLiftCurveSlope(A, Mcruise_calc(Vcruise, Tcruise, gamma, R), eta), alpha0, alpha1, alpha0L)
    plt.plot(np.arange(alpha0, alpha1), CLalphaline)
    plt.title('Lift Coefficient (CL) vs. Angle of Attack (Wing)')
    plt.xlabel('Angle of Attack [°]')
    plt.ylabel('Lift Coefficient (CL) [-]')
    #plt.savefig('LiftCoefficientAngleofAttack.png')
    plt.grid()
    #plt.show()

    CDalphaline = CD_calc(CD0, A, OswaldEfficiencyFActor(A), Mcruise_calc(Vcruise, Tcruise, gamma, R), eta, alpha0, alpha1, alpha0L)
    plt.plot(np.arange(alpha0, alpha1), CDalphaline)
    plt.title('Drag Coefficient (CD) vs. Angle of Attack (Wing)')
    plt.xlabel('Angle of Attack [°]')
    plt.ylabel('Drag Coefficient (CD) [-]')
    #plt.savefig('DragCoefficientAngleofAttack.png')
    plt.grid()
    #plt.show()

    print('Trim Angle: ' + str(TrimAngleCalc(CLdes, WingLiftCurveSlope(A, Mcruise_calc(Vcruise, Tcruise, gamma, R), eta), alpha0L)) + ' deg')

    print('Wing Maximum Lift Coefficient Clean (CLmax): ' + str(CLmaxwing_calc(CLmaxoverClmax)) + ' -')


    print('Stall Angle: ' + str(StallAngleWing_calc(CLmaxwing_calc(CLmaxoverClmax), WingLiftCurveSlope(A, Mcruise_calc(Vcruise, Tcruise, gamma, R), eta), alpha0L, DeltaalphaCLmax)) + ' deg')


    print('Flap Chord: ' + str(cflap_calc(cfoverc, chord)) + ' m')


    print('Flap Length: ' + str(lflap_calc(SwfoverS_calc(DeltaCLmax_calc(CLmax_land, CLmaxwing_calc(CLmaxoverClmax)), DeltaClmax), OptimumWingArea(MTO, g, rho0, n_p, c, A, nmax, Vcruise, rho, Vs, CLmaxwing_calc(CLmaxoverClmax)))) + ' m')


    print('Stall Angle With Flap: ' + str(alphastall_flap_calc(StallAngleWing_calc(CLmaxwing_calc(CLmaxoverClmax), WingLiftCurveSlope(A, Mcruise_calc(Vcruise, Tcruise, gamma, R), eta), alpha0L, DeltaalphaCLmax), SwfoverS_calc(DeltaCLmax_calc(CLmax_land, CLmaxwing_calc(CLmaxoverClmax)), DeltaClmax))) + ' deg')


    print(str(dragcoeff_flap(cflap_calc(cfoverc, chord), chord, SwfoverS_calc(DeltaCLmax_calc(CLmax_land, CLmaxwing_calc(CLmaxoverClmax)), DeltaClmax))))