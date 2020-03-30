# -*- coding: utf-8 -*-
"""
Created on Wed Jun 12 17:45:15 2019

@author: Matteo
"""

def Torque_crit(lh,lw):
    
    n=3.8
    MTOW=4000
    rho=1.225
    Vb=91.65
    S=2.7924
    g=9.81
    ##### the values below are valid for an aspect ratio of 4 and airfoil NACA0018
    arm=0.169
    Cmac=-0.017
    c=0.709
    g=9.81
    ############################################################
   
    L=n*(lw/(lw+lh))*MTOW*g+Cmac*(0.5*rho*Vb**2)*S*c/(lw+lh)
    Mac=Cmac*0.5*rho*Vb**2*S*c
    T=L*arm+Mac
    
    return T,L
    
    