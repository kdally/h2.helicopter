# -*- coding: utf-8 -*-
"""
Created on Tue Jun 18 14:10:31 2019

@author: Matteo
"""
from math import sqrt
Sw=20.77
def Hstab(Sw):
    
    AR=4
    Sh=0.13*Sw
    c=sqrt((Sh)/AR)
    Se=0.275*Sh # surface area elevator
    b=AR*c
    ce=Se/b # elevator cord
    l=1-ce/c
    b05=b/2
    
    return c,b05,l,ce,Se

span=5.19
lv=6.09

def Vstab(span,lv):
    
    T=2187.77 # thrust
    rho=1.1
    cl=0.764
    V=50
    AR=1.2
    
    
    Sv=T*span/(0.5*rho*cl*lv*V**2) # vertical fin surface area
    c=sqrt(Sv/AR)
    b=AR*c
        
    cr=0.415*c # ruffer cord
    l=(1-cr/c)
    arm=l/2
    
    L=0.5*rho*V**2*cl*Sv
    To=L*arm # torque
    
    return c,b,l,cr,Sv,L,To

Hstab(Sw)
    