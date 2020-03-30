# necessary plotting input
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

""""Powers Curve Template"""

# size the figure frame: the combination [10,7.5] and the fontsizes with the scaling '0.85\linewidth' in latex works best
fig=plt.figure(figsize=(10,7.5))

# initialize single subplot
ax = plt.subplot(111)


# here set what are your x and y coordinates 
x = 0

# plot all the graphs in one figure
plt.plot(x,df_power_diagram['actual_induced_power_updated']/1000, 'r-')
plt.plot(x,df_power_diagram['actual_parasite_power']/1000, 'r-.')
plt.plot(x,df_power_diagram['actual_profile_power']/1000, 'r:')
plt.plot(x,df_power_diagram['total_power']/1000, 'b--')
plt.legend(['Induced Power', 'Parasite Power', 'Profile Power', 'Total Power'], fontsize=18)


# set the y tick format such that thousands are separated by commas for easier reading
ax.set_yticklabels(['{:,}'.format(int(x)) for x in ax.get_yticks().tolist()])

# set the axis labels; use latex format for units 
plt.xlabel(r'Speed ($\frac{m}{s}$)', fontsize=18)
plt.ylabel(r'Induced Power ($W$)',fontsize=18)

# set the font size, style and
# insert grid
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.rc('font', family='serif')
plt.grid(True)
plt.show()
# save picture as eps
plt.savefig('power_curve_embedded.eps')

""""Bar Plot Template"""

# size the figure frame: the combination [10,7.5] and the fontsizes with the scaling '0.85\linewidth' in latex works best
fig=plt.figure(figsize=(13,8))

# initialize single subplot
ax = plt.subplot(111)

# set the values to plt in the bar chart

power_phases=pd.Series({'Hover': P_in/1000,'Climb':climb_power/1000, 'Cruise':cruise_power/1000})

power_phases.plot(kind='bar')

# set the y tick format such that thousands are separated by commas for easier reading
ax.set_yticklabels(['{:,}'.format(int(x)) for x in ax.get_yticks().tolist()])

# set the axis labels; use latex format for units 
plt.xlabel(r'Flight Phases', fontsize=18)
plt.ylabel(r'Power ($kW$)',fontsize=18)

# set the font size, style and insert grid
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.grid(True)
plt.savefig('power_phases.eps')
