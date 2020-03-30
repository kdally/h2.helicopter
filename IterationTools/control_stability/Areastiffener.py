# -*- coding: utf-8 -*-
"""
Created on Thu Jun 13 08:44:38 2019

@author: Matteo
"""

def Area_stiff(tstr):
    
    A1=0.0075*tstr
    A2=0.015*tstr*2
    A3=0.0075*tstr*2
    
    Astiff=A1+A2+A3
    
    return Astiff