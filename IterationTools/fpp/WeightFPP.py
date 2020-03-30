##
import numpy as np
import scipy as sp
import math
import matplotlib.pyplot as plt
from IterationTools import Mission

cruise_altitude = 4000 #m
aircraft = Mission.TiltRotor(cruise_altitude, 350, 270, 220, 4000., 9.80665, 8, 45, 6.17, 0.9, 1, 1.15, 0.75)
aircraft.performFlight(1.225)
data = aircraft.data
timestamps=aircraft.timestamps
powert=data['powert']
energy=data['energy']
time=data['time']
powert=np.around(powert, decimals=1)
time=np.around(time, decimals=1)

##
# PREPARE DATA FOR MISSION
Pmin= 0
TOe = np.around(timestamps['cruise speed'], decimals=1)
#tdesign = np.around(timestamps['Min Height'], decimals=1)
tendcruise = np.around(timestamps['Min Height'], decimals=1)
#TOe= 680
Toend=np.where(time==TOe)[0]
# EnergyTO=trapz(time(1:Toend),powert(1:Toend)-PfcTO)/(60*60)   #kwh
Pmax=np.amax(powert)
ratioenergies=(59.9+2831.3-265)/(time[-1]-3031.4)
#tdesign=3897.8
#tendcruise=time[-1]#3790
##
time=time-0.1
#tc=tendcruise-tdesign
dc=10
#inindex=np.where(time==tdesign)[0]
end2=np.where(time==tendcruise)[0]
teto=np.array(time[::10])
Pmfc=np.linspace(Pmin,Pmax,len(teto))
Pmfc=np.array([Pmfc])
energyto=np.array([energy[::dc]])
powertTO=np.array([powert[::dc]])
##
#START REAL PROGRAM
#tempoTO=time[0]:((time[-1]-time[0])/(size(teto,2)-1)):time[-1];
teto2, Pmfc2=np.meshgrid(teto,Pmfc)
delta, Pmfc3=np.meshgrid(energyto,Pmfc)
RHS=(powertTO-Pmfc2)*(teto[2]-teto[1])
Toend=np.around(Toend/10,decimals=0)
Toend=Toend.astype(int)[0]
end2=np.around(end2/10,decimals=0)
end2=end2.astype(int)[0]
RHS[:, Toend:end2]=RHS[:, Toend:end2]+0.3*Pmfc.T
##
#CHECK REQUIRMENTS
RHS[RHS==0]=1e-2
RHS[RHS<0]=0
RHS2=RHS
RHS=np.cumsum(RHS,axis=1)
RHS=RHS[:,-1]
res=RHS
##
activitytime=np.empty(len(res))
for i in range(len(res)-1):
    activitytime[i]=teto[[j for j, e in enumerate(RHS2[i,:]) if e != 0][-1]]
##
activitytime[-1]=0
activitytime=np.array([activitytime])
teto=np.array([teto])
res=np.array([res])

##
#plot res
fig = plt.figure()
plt.plot(teto.T,powertTO.T,'r')
plt.scatter(activitytime,Pmfc)
# ax.view_init(elev=30., azim=240)
# ax.set_xlabel('$X(m)$', fontsize=10, rotation=0,linespacing=8.2)
# ax.set_ylabel('$Z(m)$',fontsize=10)
# ax.yaxis.labelpad=10
# ax.xaxis.labelpad=15
# ax.set_title('Aileron deflection');
plt.show()
##
fig = plt.figure()
plt.plot(teto.T,energyto.T)
plt.plot(teto.T,np.fliplr(res).T)
plt.show()
##
#Main propulsion system
eff=0.4477
SEh2=141.8  #MJ/kg
mjtokwh=0.277778
SWmp=1.73+1  #kg/kg
SEmp=0.5682 #kg/kw
SPcc=0.195  #kg/kw
SEsc= 4E-3  #kwh/kg
SPsc= 10    #kw/kg 2.5
PSF=1.3

##
#battery data
SEb= 0.460
SPb= 3
##
#WEIGHT BATTERIES
PSF=1.3
wbe1=(RHS)*mjtokwh*1E-3/SEb
wbe2=(Pmax-Pmfc)/SPb;
wbe1=np.array([wbe1])
#wbe2=np.array([wbe2])
wbe22=wbe2*PSF;
wsc1=(RHS)*mjtokwh*1E-3/SEsc
wsc2=(Pmax-Pmfc)/SPsc
wbefinal=np.maximum(wbe1, wbe2)
wbefinal2=np.maximum(wbe1, wbe22)

##
plt.plot(Pmfc.T, wbe1.T)
plt.plot(Pmfc.T, wbe2.T)
plt.show()

##
# Weight of main propulsion system
energymp=Pmfc*(TOe+time[-1]-tendcruise+0.7*(tendcruise-TOe))
wemp=energymp/eff/(SEh2*1E3)*(SWmp)
wpmp=Pmfc*SEmp
wpmp2=wpmp*PSF
wpcc=Pmax*SPcc*np.ones(Pmfc.size)
wpcc2=wpcc*PSF
Wtot=wpmp+wbefinal+wemp+wpcc
Wtot2=wpmp2+wbefinal2+wemp+wpcc2

##
plt.plot(Pmfc.T,wpmp.T,linewidth=5)
plt.plot(Pmfc.T,wemp.T,linewidth=5)
plt.plot(Pmfc.T,wbefinal.T,linewidth=5)
plt.plot(Pmfc.T,wpcc.T,linewidth=5)
plt.plot(Pmfc.T,Wtot.T,linewidth=5)
plt.show()
##
plotforme=np.concatenate((361.6*np.ones(len(time[1:Toend*10:dc])), 0.7*361.6*np.ones(len(time[Toend*10:end2*10:dc])),361.6*np.ones(len(time[end2*10::dc]))), axis=None)
plotforme=np.array([plotforme])
plt.plot(time,powert)
plt.plot(teto.T,plotforme.T)
plt.ylim((0,1300))
plt.show()
# xlabel("Time [s]")
# ylabel("Power required[kw]")
# lgd=legend("Power Required","Fuel Cell output power","location","North");
# lgd.FontSize=30;
# set(gca,'FontSize',30)
# ylim([0 1300])

##
#final results

M=np.amin(Wtot)
I=np.where(Wtot==np.amin(Wtot))[1][0]
pesopsf=Wtot2[0,I]
POWERRRRR=Pmfc[0,I]
energybatterykwh=RHS[I]*mjtokwh*1E-3
energybatterykj=RHS[I]*1E-3
wbefinal[0,I]
volumedensity=620
volumeb=(RHS[I])*mjtokwh*1E-3/volumedensity
fuelweight=1/2.73*wemp[0,I]




##

