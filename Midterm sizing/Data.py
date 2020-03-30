import math
import numpy as np


def density(height):
    temp = 288.15-0.0065*height
    pressure = 101325*math.pow(temp/288.15,(-9.80665/(-0.0065*287.05)))
    return pressure/(287.05*temp)

S    = 12   # m2
C_d  = 0.04 # -
C_dr = 0.3  # -
MTOM = 4000 # kg
g    = 9,80665 # m/s2
MTOW = MTOM * g # N
C_l  = 0.6 # -
