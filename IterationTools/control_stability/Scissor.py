# -*- coding: utf-8 -*-
"""
Created on Wed Jun  5 14:39:06 2019

@author: Matteo
"""

from math import sqrt,pi,tan
import numpy as np
import matplotlib.pyplot as plt 

def Scissor(c,Ah,Lambda_ch,Cl_alphaw,Sw,Cm_acf,Cm_acw,lh,b,Cl_f,Cl_alphaf,xac_w,SM,Cl_w,Sf,lf):
    
    # this function returns the scissor plot diagram
    # c= MAC cord
    # Ah= aspect ratio tail 
    # Lambda_ch= sweep angle taken at 0.5 of the tail root cord
    # Cl_alphaw=  cl-alpha curve gradient  main wing 
    # Sw= wing surface area 
    # Cmacw= moment coefficient at the areodinamic center given the contribution of the wing 
    # Cmacf= moment coefficient at the areodinamic center given the contribution of the fuselage
    # lh= tail arm 
    # b=span
    # Cl_alphaf = cl-alpha curve gradient  fuselage
    
    position = lh/abs(lh)
    
    eta=0.95
    beta=0.957
    Cl_h=-0.8*position
    
    
    xcg= np.arange(-3,3,0.01)
    Cl_alphah=(2*pi*Ah)/(2+sqrt(4+((Ah*beta)/eta)**2+(1+(tan(Lambda_ch)**2/beta**2))))
   
    ms= 1/((Cl_alphah/Cl_alphaw)*lh/c)
    #qs=-(xac_w-SM+(Cl_alphaf/Cl_alphaw)*(Sf/Sw)*lf/c)/((Cl_alphah/Cl_alphaw)*lh/c)
  
    

    
    mc=1/((Cl_h/Cl_w)*lh/c)
    #qc=((Cm_acw/Cl_h)-xac_w-(Cm_acf/Cl_w)-(Cl_f/Cl_w)*(Sf/Sw)*lf/c)/((Cl_h/Cl_w)*lh/c)
    

    qs=-(xac_w-SM-position*(Cl_alphaf/Cl_alphaw)*(Sf/Sw)*lf/c)/((Cl_alphah/Cl_alphaw)*lh/c)
    qc=((Cm_acw/Cl_h)-xac_w+position*(Cm_acf/Cl_w)+position*(Cl_f/Cl_w)*(Sf/Sw)*lf/c)/((Cl_h/Cl_w)*lh/c)
    
    
    
# =============================================================================
#     elif lf>0:
#         qs=-(xac_w-SM-(Cl_alphaf/Cl_alphaw)*(Sf/Sw)*lf/c)/((Cl_alphah/Cl_alphaw)*lh/c)
#         qc=((Cm_acw/Cl_h)-xac_w+(Cm_acf/Cl_w)+(Cl_f/Cl_w)*(Sf/Sw)*lf/c)/((Cl_h/Cl_w)*lh/c)
# =============================================================================
        
   
    
    
    # print(mc)
    # print(ms)
    # print(Cl_alphah)
    # print(Cl_alphaw)
    # print(Cl_h)
    # print(cl_w)
    
    m=[mc,ms]
    q=[qc,qs]


    # Stability curve
    # ShS_stab=ms*xcg+qs
    # Controllability curve 
    # ShS_contr=mc*xcg+qc

    # plt.plot(xcg,ShS_stab,label="stability")
    # plt.plot(xcg,ShS_contr,label="controllability")
    # plt.ylim([0,0.5])
    # plt.xlim([-2,0.5])
    # plt.legend()
    # plt.show()
    
    return(m,q)


def main():
    pass
    
if __name__ == "__main__":
    c=3
    Ah=30
    Lambda_ch=20*pi/180
    Cl_alphaw=6.0
    Cm_acw=-0.05 # moment coefficient aerodynamic center 
    Cm_acf=-0.15# moment coefficient aerodynamic center
    lh=-5
    lf=-0.5
    b=14
    xac_w= 0.25 # position of the ac of the main wing as percentage of the MAC
    SM=0.04 # static margin
    Cl_w=0.4
    Cl_f=0.2
    Cl_alphaf=pi
    Sw=20
    Sf=32
    # m,q = Scissor(c,Ah,Lambda_ch,Cl_alphaw,Sw,Cm_acf,Cm_acw,lh,b,Cl_f,Cl_alphaf,xac_w,SM,Cl_w,Sf,lf)
    # print('m',m)
    # print('q',q)

    # Scissor(c,Ah,Lambda_ch,Cl_alphaw,Sw,Cm_acf,Cm_acw,lh,b,Cl_f,Cl_alphaf,xac_w,SM,Cl_w,Sf,lf)



    
    
        