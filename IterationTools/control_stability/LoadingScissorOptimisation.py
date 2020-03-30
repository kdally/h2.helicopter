from scipy.optimize import fmin
import matplotlib.pyplot as plt
import numpy as np

def checkSwitch(m,b):
    """
    If the slope of the stability curve is negatvive then switch the order of the slopes lists so that the correct intersections can be found
    Input is slope and intercept lists
    Return updated slope and intercept lists
    """

    #check if the slope of the stability curve is negative and switch stability and controllability curve
    switch = False
    if m[1] < 0:
        m[0], m[1] = m[1], m[0]
        b[0], b[1] = b[1], b[0]
        switch = True

    return m, b, switch

def findScissorOptimum(m,b):
    """ 
    Take in list of slopes (m) and list of y intercepts (b) from scissor plot in order [control,stability,aft,forward]
    Returns, in order: lower xc.g. range, higher xc.g. range, Sh/S, xLEMAC (All floats) and a list of the 2 optimum scaling factors

    Written by Benjamin
    """
    # assert m[0]<m[1], "Stability Curve more negative than Control Curve"
    m, b, switch = checkSwitch(m, b)

    def minimizeForKandC(x):
        """
        Function to scale load diagrams curves, find intercepts with stability and control axis and return difference between y intercepts (absolute value)
        Function is within a function so it can access lists m and b, whilst only taking argument x for the scipy.optimize fmin function
        Takes in list of scaling factors [multiplier,intercept]
        Returns difference between y intercepts (float)
        """

        #scale load diagram plots 
        m3 = x[0]*m[2]
        b3 = x[0]*b[2]+x[1]
        m4 = x[0]*m[3]
        b4 = x[0]*b[3]+x[1]

        #intercepts
        x1 = (b4-b[0])/(m[0]-m4)
        y1 = x1*m4+b4
        x2 = (b3-b[1])/(m[1]-m3)
        y2 = x2*m3+b3

        return abs(y2-y1)
    
    minimum  = fmin(minimizeForKandC,[3,-1],disp = False)[0:2] #list of scaling factors [k,c]

    #scale loading diagram by factors to the stability diagram
    m3 = minimum[0]*m[2]
    b3 = minimum[0]*b[2]+minimum[1]
    m4 = minimum[0]*m[3]
    b4 = minimum[0]*b[3]+minimum[1]

    #find intercepts between loading diagram curves (scaled) and control and stability curves
    x1 = (b4-b[0])/(m[0]-m4)
    y1 = x1*m4+b4
    x2 = (b3-b[1])/(m[1]-m3)
    y2 = x2*m3+b3
        
    yavg = (y1+y2)/2. #average y intercept values for one value of Sh/S
    xLEMAC = (yavg-minimum[1])/minimum[0] #convert back to loading diagram scale to find equivalent Sh/S

    # assert yavg>0, 'Sh/s not positive'

    return [min(x1,x2),max(x1,x2)],yavg,xLEMAC,minimum


def findScissorCG(m,b,uppercg):
    """
    Find the control characteristics given a boundary on the cg range
    input: list of slopes (m) and list of y intercepts (b) from scissor plot in order [control,stability,aft,forward], upper bound (most aft) on cg (float)
    output: cgrange [lower bound, upperbound], Sh_s, xLEMAC
    """
    xLEMAC = m[2]*uppercg+b[2]
    lowercg = (xLEMAC - b[3])/m[3]

    #case 1 go along control with forward cg value
    ycaseI = m[0]*lowercg + b[0]
    xtestI = (ycaseI - b[1])/m[1]
    caseI = xtestI > uppercg

    #case 2 go along stability with aft cg value
    ycaseII = m[1]*uppercg+b[1]
    xtestII = (ycaseII-b[0])/m[0]
    caseII = xtestII < lowercg

    cII = ycaseII - xLEMAC
    cI = ycaseI - xLEMAC


    #conditions
    if caseII:
        return [lowercg,uppercg], ycaseII, xLEMAC, cII
    elif caseI:
        return [lowercg,uppercg], ycaseI, xLEMAC, cI
    return None
    
def plotCustomCG(m,b,cg,Sh_S,c):
    labels = ['Control','Stability','Aft','Forward']
    graphx = np.linspace(-2,2,100)
    graphy = []
    b[3] = b[3] + c
    b[2] = b[2] + c
    for i in range(4):
        graphy.append(m[i]*graphx+b[i])

    for i in range(len(graphy)):
        plt.plot(graphx,graphy[i],label = labels[i])

    plt.hlines(y = Sh_S, xmin = cg[0], xmax = cg[1], linestyle = '--', label = 'Optimum') #horizontal line to show optimum point in scissor plot
    plt.ylim([-0.3,0.4])
    plt.ylabel('Sh/S', size = 14)
    plt.xlabel('xcg/MAC', size = 14)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.legend()
    plt.show()
    return

def plotScissorOptimum(m,b,scale,cg,Sh_S):
    """ 
    Plot superimposed optimum graph for scissor plots
    Takes list of slopes and intercepts from scissor and loading diagram plots [control,stability,aft,forward], a list of scaling factors [k,c] for the loading diagram plots
    the lower and upper bounds of c.g. location (floats) and Sh_S value (float).

    Written by Benjamin
    """
    
    m, b, switch = checkSwitch(m, b)
    if not switch:
        labels = ['Control','Stability','Aft','Forward'] #labels used for legend in graphing
    else:
        labels = ['Stability','Control','Aft','Forward']

    graphx = np.linspace(-2,2,100) #xrange of values
    graphy = [] #empty list for graphing
    

    #scale loading diagram slope and intercept
    m3 = scale[0]*m[2]
    b3 = scale[0]*b[2]+scale[1]
    m4 = scale[0]*m[3]
    b4 = scale[0]*b[3]+scale[1]

    #change values in list of slopes
    m[2:] = [m3,m4]
    b[2:] = [b3,b4]

    #for loop to calculate y values for each curve
    for i in range(4):
        graphy.append(m[i]*graphx+b[i])

    #plot each curve
    for i in range(len(graphy)):
        plt.plot(graphx,graphy[i],label = labels[i])

    plt.hlines(y = Sh_S, xmin = cg[0], xmax = cg[1], linestyle = '--', label = 'Optimum') #horizontal line to show optimum point in scissor plot
    plt.ylim([-0.3,0.4])
    plt.ylabel('Sh/S', size = 14)
    plt.xlabel('xcg/MAC', size = 14)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.legend()
    plt.show()
    return

def main():
    pass

if __name__ == "__main__":   
    # m = [-1.988377394,0.281167532,-0.116798638,-0.118237114] #test m values from ADSEE excel 
    # b = [0.640175101,-0.047689749,0.562063014,0.527569664] #test b values from ADSEE excel

    m = [-1.6,-0.4,-0.116798638,-0.118237114] #test m values for canard 
    b = [0.640175101,0.25,0.562063014,0.527569664] #test b values for canard

    # cg,Sh_S,xLEMAC,scale = findScissorOptimum(m,b)

    # plotScissorOptimum(m,b,scale,cg,Sh_S)