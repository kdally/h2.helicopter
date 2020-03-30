import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

data1 = pd.read_excel('airfoils_short.xlsx', sheet_name='0012').values
alfa1 = data1[:, 0]
cl1 = data1[:, 1]
cd1 = data1[:, 2]

data2 = pd.read_excel('airfoils.xlsx', sheet_name='0012').values
alfa2 = data2[:, 0]
cl2 = data2[:, 1]
cd2 = data2[:, 2]

data3 = pd.read_excel('validation_viterna.xlsx').values
alfa3 = data3[180:, 0]
# cl3 = data2[:, 1]
cd3 = data3[180:, 1]

data4 = pd.read_excel('validation_viterna_cl.xlsx').values
alfa4 = data4[:, 0]
# cl3 = data2[:, 1]
cl4 = data4[:, 1]


data5 = pd.read_excel('airfoils_try.xlsx', sheet_name='0012').values
alfa5 = data5[180*4:, 0]
cl5 = data5[180*4:, 1]
cd5 = data5[180*4:, 2]
print(alfa5)

ticksize = 12
axislabels = 14


# plt.rc('text', usetex=True)
# plt.rc('font', family='serif')

# font = {'fontname':'Book Antiqua'}
# mpl.rc('font', family = 'Book Antiqua')

fs = 34

plt.figure(0)
plt.ylabel(r'$C_{d}$', fontsize=fs)#, **font)
plt.xlabel(r'$\alpha$', fontsize=fs)#, **font)
plt.xticks(fontsize=fs)
plt.yticks(fontsize=fs)

# plt.plot(alfa2, cd2, color='r', label= '')
plt.plot(alfa3, cd3, color='k', label= 'Experimental data', linewidth=3)
plt.plot(alfa5, cd5, color='b', label='Theoretical data', linewidth=3)
plt.grid()
plt.legend(fontsize=fs, loc=1)

plt.figure(1)
plt.ylabel(r'$C_{l}$', fontsize=fs)#, **font)
plt.xlabel(r'$\alpha$', fontsize=fs)#, **font)
plt.xticks(fontsize=fs)
plt.yticks(fontsize=fs)

# plt.plot(alfa2, cl2, color='r', label= 'used')
plt.plot(alfa4, cl4, color='k', label= 'Experimental data', linewidth=3)
plt.plot(alfa5, cl5, color='b', label= 'Theoretical data', linewidth=3)
plt.grid()
plt.legend(fontsize=fs, loc=1)
# plt.ylabel(r'C_{l}')
# plt.xlabel(r'\alpha')
plt.show()
