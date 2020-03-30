import math
from scipy.optimize import fmin, fsolve
import numpy as np
import matplotlib.pyplot as plt
import structural_loads as st
import Areastiffener as AreaSt
import Torque_tail as tt
# import control_stability.structural_loads as st
# import control_stability.Areastiffener as AreaSt
# import control_stability.Torque_tail as tt

#materials
AL7075 = {'E':72.7e9, 'nu':0.33, 'sigmay':424e6, 'price':3.98, 'density': 2810, 'G':27.4e9}
AL2024 = {'E':73.1e9, 'nu':0.33, 'sigmay':360e6, 'price':2.16, 'density': 2780, 'G':28.5e9}


def bendingMomentStress(t,b,M,h,l,tspar,Istif,Astif):
    n = math.ceil(l/b)-1
    Istring = n*2*(math.pow(h/2,2)*Astif+Istif)
    Iskin = 2*(l*t*math.pow(h/2,2)+l*math.pow(t,3)/12)
    Ispar = 2*tspar*math.pow(h,3)/12
    I = Istring + Iskin + Ispar
    return M*(h/2)/I, n

def shearFlows(Sy,T,t,tspar,n,l,h,Astif,G):
    # assume positive counter-clockwise direction, positive upward y value
    # calculate moment of inertia
    Istring = n*2*(math.pow(h/2,2)*Astif+Istif)
    Iskin = 2*(l*t*math.pow(h/2,2)+l*math.pow(t,3)/12)
    Ispar = tspar*math.pow(h,3)/12
    Ixx = Istring + Iskin + Ispar


    Bnormal = t*(l/n)+Astif # stringer booms
    Bcorner = Astif + (t*(l/n)/2) + (tspar*h/6) # spar cap booms
    # Bnormal = 0.0012
    # Bcorner = 0.0009
    B = [Bcorner] + [Bnormal]*n + [Bcorner]*2 + [Bnormal]*n # list of booms starting from bottom left corner to last before top left
    y = [-h/2]*(2+n) + [h/2]*(1+n) # distances in y direction from centroid to booms
    q_base = [] # initialise list for base shear flows
    temp = 0 # temporary holder for shear flow from previous step

    # loop calculating shear flow in each gap between booms
    for i in range(len(B)): 
        qb = (-Sy/Ixx)*B[i]*y[i]
        q_base.append(qb+temp)
        temp += qb
    # print('q_base', q_base)
    
    # set up for finding qs0
    s = [l/(n+1)]*(n+1) + [h] + [l/(n+1)]*(n+1) # length of each interval between booms
    l_list = [h/2]*(n+1) + [l/2] + [h/2]*(n+1) # moment arms of different segments between booms 
    sum_l = 0 # sum counter for calculating torque shear flow around area

    # sum up shear flows
    for i in range(len(l_list)):
        sum_l += q_base[i]*l_list[i]*s[i]

    qs0 = (T-sum_l)/(2*l*h) #shear flow at cut
    q = q_base.copy() #create copy of base shear flows to calculate actual shear flows in segments
    q.append(0) # base shear flow on cut segment
    s.append(h) # length of cut segment

    # add shear flow at cut to each base shear flow
    for i in range(len(q)):
        q[i] += qs0
    
    #intialise rate of twist
    rate_of_twist = 0
    ts = [t] * (n+1) + [tspar] + [t]*(n+1) + [tspar] #thickness of different segments

    #sum up components of loop for rate of twist
    for i in range(len(l_list)):
        rate_of_twist += q[i]*s[i]/ts[i]

    #divide by 2AG to find rate of twist
    rate_of_twist = rate_of_twist/(2*l*h*G)

    #find stress in each section
    stress = []
    for i in range(len(q)):
        stress.append(q[i]/ts[i])
    return q, rate_of_twist, stress

def equations(x, M, C, h, l, tspar, tstr, Istif, material):
    
    sigma_cr_else,sigma_cc,two_welse,Astif = st.Stress(C,material['nu'],material['E'],tstr,material['sigmay'])
    bendingstress = bendingMomentStress(x[0],x[1],M,h,l,tspar,Istif,Astif)
    cripplingstress = st.Sigma_panel(x[0],x[1],sigma_cr_else,sigma_cc,two_welse,Astif)
    f1 =  bendingstress - material['sigmay']
    f2 =  cripplingstress - material['sigmay']
    # print('simga else',sigma_cr_else)
    # print('sigma cc', sigma_cc*1e-6)
    # print('two_welse', two_welse)
    if x[0] <= 0 or x[0] >= 0.01:
        f1 = 1000000000000000
        f2 = 1000000000000000
    if x[1] <= 0 or x[1] >= 0.5:
        f1 = 1000000000000000
        f2 = 1000000000000000
    print('bendingstress', bendingstress*1e-6)
    print('crippling stress', cripplingstress*1e-6)
    print('f1', f1*1e-6)
    print('f2', f2*1e-6)
    return [f1,f2]

def area(t,n,h,l,tspar,Astif):
    Askin = l*t*2
    Aspar = h*tspar*2
    Astif_tot = Astif*n*2
    return Askin, Aspar, Astif_tot, Askin+Aspar+Astif_tot

def cost_per_length(area,material):
    return area*material['price']*material['density']

def vonMisesStress(bendingstress,stress):
    return math.sqrt(bendingstress**2+6*max(stress)**2)

def main():
    pass

if __name__ == "__main__":
    C = 4
    tstr = 0.002
    lh = 6.41
    lw = 0.1*2
    T, L = tt.Torque_crit(lh,lw)
    M = L/2
    h = (0.0511+0.08912)*0.819/2
    l = (0.725-0.25)*0.819
    print('h',h,'l',l)
    print('T',T,'L',L)
    tspar = 0.0075
    Istif = 0
    material = AL2024
    Astif = AreaSt.Area_stiff(tstr)

    # print(fsolve(equations,[0.002, 0.05],(M,C,h,l,tspar,tstr,Istif,material)))
    # #test shear flow calculation with textbook case
    # q, rate_of_twist = shearFlows(400*9.81, -400*9.81*0.3, 0.003, 0.075, 4, l, h, Astif, material['G'])
    # print('q', q)
    # print('rate of twist', rate_of_twist)

    xlow, xhigh = 0.0008, 0.003
    yhigh, ylow = 0.2, 0.05
    t = np.linspace(xhigh,xlow,200)
    b = np.linspace(yhigh, ylow,200)

    tlist = []
    blist = []
    bs = []
    cs = []
    vonmis = []
    shears = []
    sol = []
    aarr = []
    mass = []
    sigma_cr_else, sigma_cc, two_welse, Astif = st.Stress(C, material['nu'], material['E'], tstr, material['sigmay'])
    Askin, Aspar, Astif_tot, a = area(0.00125, 2, h, l, tspar, Astif)

    for i in t:
        for j in b:
            bendingstress, n = bendingMomentStress(i, j, M, h, l, tspar, Istif, Astif)
            _,_,stress = shearFlows(L/2,T,i,tspar,n,l,h,Astif,material['G'])
            vonmises = vonMisesStress(bendingstress, stress)
            cripplingstress = st.Sigma_panel(i, j, sigma_cr_else, sigma_cc, two_welse, Astif)
            Askin, Aspar, Astif_tot, a = area(i, n, h, l, tspar, Astif)
            m = a*material['density']*4
            shears.append(max(stress)/(1e6))
            vonmis.append(vonmises/(1e6))
            tlist.append(i)
            blist.append(j)
            bs.append(bendingstress/(1e6))
            cs.append(cripplingstress/(1e6))
            sol.append(cripplingstress/(1e6)-bendingstress/(1e6))
            aarr.append(a)
            mass.append(m)

    plt.subplot(1,3,1)
    plt.scatter(tlist,blist, c=cs, cmap='Oranges')
    cb = plt.colorbar()
    cb.ax.tick_params(labelsize=12)
    plt.xlabel('Skin Thickness [m]', size = 14)
    plt.ylabel('Stiffener Pitch [m]', size = 14)
    plt.title('Buckling Stress [MPa]', size = 14)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.xlim(xlow,xhigh)
    plt.ylim(ylow,yhigh)
    # plt.show()
    plt.subplot(1,3, 2)
    plt.scatter(tlist, blist, c=vonmis, cmap='Oranges')
    cb = plt.colorbar()
    cb.ax.tick_params(labelsize=12)
    plt.xlabel('Skin Thickness [m]', size = 14)
    plt.title('Von Mises Stress [MPa]', size = 14)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.xlim(xlow, xhigh)
    plt.ylim(ylow, yhigh)
    # plt.show()
    # plt.subplot(2,3, 3)
    # plt.scatter(tlist, blist, c=aarr, cmap='Blues')
    # plt.colorbar()
    # plt.xlim(xlow, xhigh)
    # plt.ylim(ylow, yhigh)

    plt.subplot(1, 3, 3)
    plt.scatter(tlist, blist, c=mass, cmap='Blues')
    cb = plt.colorbar()
    cb.ax.tick_params(labelsize=12)
    plt.xlabel('Skin Thickness [m]', size = 14)
    plt.title('Mass [kg]', size = 14)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.xlim(xlow, xhigh)
    plt.ylim(ylow, yhigh)

    plt.show()
    
    
    Askin, Aspar, Astif_tot, a = area(0.0015, 2, h, l, tspar, Astif)
    print(Askin, Aspar, Astif_tot, a)
    # print('skin', Askin*material['density']*1.5)
    # print('spar', Aspar*material['density']*1.5)
    # print('stif', Astif_tot*material['density']*1.5)
    m = a*material['density']*1.638*2*1.5
    print('mass', m)