##
import numpy as np
##
#fuel cell
Vfccell=0.65
Ifccell=450
nfccell=215#335
nseries=2
nparallel=4 #2
Vfc= nfccell*Vfccell

Vfcstack=Vfc*nseries
Ifcstack=Ifccell*nparallel
print('Vfcstack',Vfcstack)
print('Ifcstack',Ifcstack)
##
#batteries
space=2
Vbatt=3.5
pmax=717.4E3#851E3
rho=1000
Ispec=1400
minseries=np.ceil(Vfcstack/Vbatt)
Ireq=pmax/(minseries*Vbatt)

print('Ireq', Ireq)
diametro=21
storage= 77.26E3
l=70.5
V=(diametro/2)**2*np.pi*l
Vl=V*10**(-6)
ec=Vl*rho
spe=0.448
rows=np.ceil(storage/(ec*minseries))
width= ((rows+1)/2)*diametro+space*(rows/2+1)
lenght= (minseries/2)*diametro+space*(minseries/2+1)
print('series',minseries)
print('parallel', rows)
print('width', width)
print('lenght',lenght)
crow=1400*(ec/spe/1000)
currenttotb=crow*rows
mcell=Ireq/Ispec*1/rows
print('mass cell', mcell)
print('vbatt', minseries*Vbatt)
print('current batt for motor', Ireq)
print('current batt', currenttotb )
##


