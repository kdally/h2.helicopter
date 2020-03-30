import numpy as np
import matplotlib.pyplot as plt
rad_data = np.loadtxt('radiator.txt')

# print(rad_data[:,0])


def htc(v): # function of speed (m/s)
    # Calculated @ 20 degC
    h = (25.87/1000/0.064)*0.0037*(v*0.064/15.06e-06)**(2/5)*0.708**(1/3)
    return h


v = rad_data[:, 0]/60/(0.387*0.212)
# v2 = np.linspace(0, 20, 100)
print(v[-1])
h1 = htc(v)
Ae_arr = rad_data[:, 1]/h1/(0.212*0.387)
vnew = 15
rad2 = 300000/40
Ae_new = rad2/htc(vnew)/(1.5*1.45)
print(Ae_arr[-1])
print(Ae_new)
print(Ae_arr[-1]*htc(15)*(1.5*1.45)*40)
# Ae = np.average(Ae_arr)
# Ae = Ae_arr[-1]
# plt.plot(rad_data[:, 0], rad_data[:, 1])
# plt.plot(v2*60*(0.387*0.212), htc(v2)*Ae)
# plt.show()
