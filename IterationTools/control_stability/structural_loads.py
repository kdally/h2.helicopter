# -*- coding: utf-8 -*-
"""
Created on Fri Jun  7 10:35:58 2019

@author: Matteo
"""
from math import pi,sqrt
import numpy as np

def Stress(C,nu,E,tstr,sigmay):   
    sigma_cr_else= (C*pi**2*E)/(12*(1-nu**2)) 
    



## stiffener strenght 
    
    alpha=0.8
    n=0.6
    
    Cstr=np.array([0.425,4.0,4.0])
    bstr=np.array([0.015,0.030,0.020])
    A1=0.015*tstr
    A2=0.030*tstr
    A3=0.015*tstr
    sigma_cc=alpha*(((Cstr/sigmay*((pi**2)*E)/(12*(1-nu**2))*(tstr/bstr)**2))**(1-n))*sigmay
    
    count=0
    for i in sigma_cc:
        
        if i/sigmay>1:
            sigma_cc[count]=sigmay
        
        count=count+1
    
    sigma_cc=(2*sigma_cc[0]*A1+2*sigma_cc[1]*A2+sigma_cc[2]*A3)/(2*A1+2*A2+A3)
    
    
    Astiff=A1+A2+A3
    
    two_welse=sqrt((C*pi**2)/(12*(1-nu**2)))*sqrt(E/sigma_cc) 
    
    return sigma_cr_else,sigma_cc,two_welse,Astiff

def  Sigma_panel(tskin,bskin,sigma_cr_else,sigma_cc,two_welse,Astiff):
    
    sigma_panel=(((Astiff+two_welse*tskin**2)*sigma_cc+(bskin-two_welse*tskin)*(tskin**3/bskin**2)*sigma_cr_else)/(Astiff+bskin*tskin))
    
    return sigma_panel
    
    
    
    
    
    
    
    
    
