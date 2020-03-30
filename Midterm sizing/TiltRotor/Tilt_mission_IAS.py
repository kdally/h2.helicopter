import math
import matplotlib.pyplot as plt
import scipy.io
import numpy as np
import pickle
from matplotlib import rc
import pandas as pd
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
# from matplotlib import style
# style.use('classic')

def density(height):
    temp = 288.15-0.0065*height
    pressure = 101325*math.pow(temp/288.15,(-9.80665/(-0.0065*287.05)))
    return pressure/(287.05*temp)

def exportMatlab(x,namex,y,namey,filename):
    scipy.io.savemat(filename, dict(namex=x,namey=y))
    return 'Data Exported'

class TiltRotor():

    def __init__(self,cruise_altitude,cruise_speed,cruise_range,climb_speed,mass,g,roc,stall,A,e,M,dl,eta_prop):
        self.cruise_altitude =  cruise_altitude #m
        self.cruise_speed = cruise_speed/3.6 #km/h
        self.climb_speed = climb_speed/3.6 #km/h speed forward in climb, also assumed in descent
        self.cruise_range = cruise_range #km    
        self.mass = mass #kg
        self.g = g #m/s
        self.weight = mass*g
        self.roc = roc #m/s rate of climb
        self.height = 0 #m
        self.energy = 0   #kJ
        self.velocityTAS = 0 #m/s  True Airspeed
        self.velocityIAS = 0 #m/s Indicated Airspeed
        self.acceleration = 0 #m/s^2 Forward Acceleration
        self.xposition = 0 #m
        self.powert = 0 #kW Total Power to shaft
        self.poweru = 0 #kW Power for lift due to tilt rotor thrust
        self.powerfor = 0 #kW #Power to overcome parasite drag
        self.poweracc = 0 #kW #Power to accelerate forward
        self.powerc = 0 #Power to climb
        self.angle = math.pi/2. #current angle of tilt rotor
        self.cl = 1 #Cl of aircraft intialised at 1
        self.cd0 = 0.017 #cd0 of aircraft assumed
        self.A = A #aspect ratio
        self.e = e #oswald efficiency
        self.S = mass*g/(self.cl*0.5*1.225*stall**2) #surface area of aircraft
        self.M = M #figure of merit
        self.N = 2 #number of rotors
        self.B = 3 #number of blades
        self.r = 0.5*math.sqrt((self.mass/(60*2))*4/math.pi) #radius of blades
        self.omega = 85 #RPS rad/s
        self.c = 0.2358
        self.dl = dl #download factor
        self.eta_prop = eta_prop #propeller efficieny
        self.lift = 0 #aircraft lift
        self.drag = 0 #aircraft drag
        self.eta_engine = 0.4 #efficiency for energy transfer from fuel to shaft
        self.eta_trans = 1#0.96
        self.eta_inst = 0.98
        self.SPEC = 141.8e6 #specific energy consumption (46.6 for jetfuel)
        self.vinduced = 0
        self.data = {'time':[], 'height':[], 'density':[], 'xposition':[], 'velocityTAS':[], 'velocityIAS':[], 'acceleration':[], 'powert':[], 'poweru': [], 'powerfor': [], 'poweracc': [], 'powerc':[], 'energy': [], 'angle': [], 'lift/drag': [], 'weight':[], 'vinduced':[]} #dictionary of lists of data
        self.timestamps = {}
    
    def findForce(self,coefficient,rho): 
        # find aerodynamic force
        return 0.5*rho*self.S*coefficient*self.velocityTAS**2
    
    def c_l_blade(self, L_rotor, rho):
        return 6*L_rotor / (self.N*self.B*rho*self.omega**2*self.c*self.r**3)

    def c_d_blade_1(self,x):
        return 20441 * x**6 -1505.9 * x**5 + 442.64 *x**4 - 63.716 * x**3 + 4.7524 * x**2 - 0.0961 * x + 0.0058

    def p_profile(self, c_d_blade, rho):
        return np.where(c_d_blade != 0.0058, self.N * self.B * c_d_blade * rho * self.omega**3 * self.r**4  * self.c /8, 0)
     
    def inducedPower(self,thrust,rho): 
        #find induced power
        induced = (thrust*self.dl/self.M)*math.sqrt((-0.5*self.velocityTAS**2)+math.sqrt(((0.5*self.velocityTAS**2)**2)+self.vHover(thrust,rho)**4))*(1/(self.eta_trans*self.eta_inst))/1000.
        self.vinduced = math.sqrt((-0.5*self.velocityTAS**2)+math.sqrt(((0.5*self.velocityTAS**2)**2)+self.vHover(thrust,rho)**4))*(1/(self.eta_trans*self.eta_inst))
        C_l_blade = self.c_l_blade(thrust,rho)
        aoaob = C_l_blade/(2*math.pi)
        C_d_blade = self.c_d_blade_1(aoaob)
        P_profile = self.p_profile(C_d_blade,rho)
        return induced+P_profile/1000.
        # return induced
   
    def testinducedPower(self,thrust,v,rho): 
        #ignore this
        return (thrust*self.dl/self.M)*math.sqrt((-0.5*v**2)+math.sqrt(((0.5*v**2)**2)+self.vHover(thrust,rho)**4))/1000.


    def vHover(self,thrust,rho): 
        #velocity induced in hover
        return math.sqrt((self.dl*thrust)/(2*self.N*rho*math.pi*self.r**2))     
    
    def velPower(self,thrust): 
        #thrust*velocity
        return thrust*self.velocityTAS/1000

    def dragCoefficient(self,cl): 
        #find drag coefficient based 
        return self.cd0 + (cl**2)/(math.pi*self.A*self.e)

    def convertToIAS(self,velocityTAS,rho,rho0):
        #convert from true airspeed to indicated airspeed
        return velocityTAS*math.sqrt(rho/rho0)

    def convertToTAS(self,velocityIAS,rho,rho0):
        #convert from indicated airspeed to true airspeed
        return velocityIAS*math.sqrt(rho0/rho)

    def findAngle(self,T_perp,T_par):
        #find angle based on trhust perpendicular and parallel
        return math.atan(T_perp/T_par)
    
    def appendData(self,time,density):
        #update lists
        self.data['time'].append(time)
        self.data['height'].append(self.height)
        self.data['density'].append(density)
        self.data['xposition'].append(self.xposition)
        self.data['acceleration'].append(self.acceleration)
        self.data['velocityTAS'].append(self.velocityTAS)
        self.data['velocityIAS'].append(self.velocityIAS)
        self.data['powert'].append(self.powert)
        self.data['poweru'].append(self.poweru)
        self.data['powerfor'].append(self.powerfor)
        self.data['poweracc'].append(self.poweracc)
        self.data['powerc'].append(self.powerc)
        self.data['energy'].append(self.energy)
        self.data['angle'].append(self.angle)
        if self.drag!= 0: self.data['lift/drag'].append(self.lift/self.drag) 
        else: self.data['lift/drag'].append(0)
        self.data['weight'].append(self.weight/self.g)
        self.data['vinduced'].append(self.vinduced)
        return

    def updateWeight(self,dt):
        #decreases weight by fuel used
        self.weight -= (self.powert*self.g*dt/(self.eta_engine*(self.SPEC/1000)))
        return

    
rho0 = 1.225 #sea level
IASswitch = True #constant to check if climb speed reached

#mission profile
hover_time = 60 #s 
cruise_altitude = 4000 #m

aircraft = TiltRotor(cruise_altitude,350,300,220,4000.,9.80665,8,60,6.17,0.9,1,1.15,0.75)


climb_acceleration = aircraft.climb_speed/80 #time to reach climbspeed (divide by seconds required.)

time = 0 #s #intitalise time
dt = 0.1 #s #time step
test_lift = True

#climb
while aircraft.height < aircraft.cruise_altitude:
    
    rho = density(aircraft.height) #find density

    #hover phase
    while time < hover_time:
        aircraft.poweru = aircraft.inducedPower(aircraft.weight,rho)
        aircraft.powert = aircraft.poweru
        aircraft.energy += aircraft.powert*dt
        time += dt
        aircraft.appendData(time,rho)
        aircraft.updateWeight(dt)
        aircraft.timestamps['Begin Climb'] = time
    
    aircraft.height += aircraft.roc*dt

    aircraft.velocityIAS = aircraft.convertToIAS(aircraft.velocityTAS,rho,rho0)

    #check if climb velocity is reached
    if (aircraft.velocityIAS < aircraft.climb_speed) and IASswitch:
        aircraft.acceleration = climb_acceleration
        T_acc = aircraft.acceleration*aircraft.mass #thrust required for acceleration
        aircraft.timestamps['IAS Climb Velocity'] = time

    else:
        IASswitch = False #set that climb speed has been reached
        temp_target_TAS = aircraft.convertToTAS(aircraft.climb_speed,rho,rho0) #set target speed in TAS
        aircraft.acceleration = (temp_target_TAS-aircraft.velocityTAS)/dt
        T_acc = aircraft.acceleration*aircraft.mass


    aircraft.velocityTAS += aircraft.acceleration*dt #update velocity

    aircraft.lift = aircraft.findForce(aircraft.cl,rho) #find lift 
    cd = aircraft.dragCoefficient(aircraft.cl)
    aircraft.drag = aircraft.findForce(cd,rho)

    T_perp = max(aircraft.weight - aircraft.lift,0.) #either set required perpendicular thrust to counteract weight or zero if lift covers it
    
    if test_lift and T_perp == 0.:
        test_lift = False
        aircraft.timestamps['Full Wing Flight'] = time

    aircraft.poweru = aircraft.inducedPower(T_perp,rho) 
    aircraft.powerfor = aircraft.velPower(aircraft.drag)/aircraft.eta_prop

    if T_perp == 0.: aircraft.cl = (aircraft.weight)/(aircraft.S*0.5*rho*aircraft.velocityTAS**2)

    
    aircraft.poweracc = aircraft.velPower(T_acc)
    aircraft.powerc = aircraft.roc*(aircraft.lift + T_perp/2)/1000  #power to climb is a component of both lift due to wing and perpendicular thrust
    aircraft.powert = aircraft.poweru+aircraft.powerfor+aircraft.poweracc+aircraft.powerc
    aircraft.energy += aircraft.powert*dt

    aircraft.xposition += aircraft.velocityTAS*dt
    aircraft.angle = aircraft.findAngle(T_perp,aircraft.drag+T_acc)

    time += dt
    aircraft.appendData(time,rho)
    aircraft.updateWeight(dt)


#set unused variables for next phase to zero
aircraft.poweru = 0
aircraft.powerc = 0

aircraft.timestamps['Cruise Altitude'] = time

#cruise
time_cruise = time
cruise_acceleration = (aircraft.cruise_speed-aircraft.convertToTAS(aircraft.climb_speed,rho,rho0))/120 #acceleration to reach cruise speed in given seconds (divisor)

while aircraft.xposition<aircraft.cruise_range*1000:

    while aircraft.velocityTAS < aircraft.cruise_speed:
        #accelerate to crusie speed
        aircraft.cl = (aircraft.weight)/(aircraft.S*0.5*rho*aircraft.velocityTAS**2)
        aircraft.lift = aircraft.findForce(aircraft.cl,rho)
        cd = aircraft.dragCoefficient(aircraft.cl)
        aircraft.drag = aircraft.findForce(cd,rho)
        aircraft.acceleration = min(cruise_acceleration,(aircraft.cruise_speed-aircraft.velocityTAS)/dt) #take minimum of cruise acceleration 
        T_acc = aircraft.acceleration*aircraft.mass
        aircraft.velocityTAS += aircraft.acceleration*dt
        aircraft.velocityIAS = aircraft.convertToIAS(aircraft.velocityTAS,rho,rho0)

        aircraft.xposition += aircraft.velocityTAS*dt
        aircraft.powerfor = aircraft.velPower(aircraft.drag)/aircraft.eta_prop
        aircraft.poweracc = aircraft.velPower(T_acc)
        aircraft.powert = aircraft.powerfor+aircraft.poweracc
        aircraft.energy += aircraft.powert*dt

        time += dt
        aircraft.appendData(time,rho)
        aircraft.updateWeight(dt)
        aircraft.timestamps['cruise speed'] = time
    
    aircraft.cl = (aircraft.weight)/(aircraft.S*0.5*rho*aircraft.velocityTAS**2)
    aircraft.lift = aircraft.findForce(aircraft.cl,rho)
    cd = aircraft.dragCoefficient(aircraft.cl)
    aircraft.drag = aircraft.findForce(cd,rho)
    aircraft.acceleration = 0
    aircraft.xposition += aircraft.velocityTAS*dt
    aircraft.poweracc = 0
    aircraft.powerfor = aircraft.velPower(aircraft.drag)/aircraft.eta_prop
    aircraft.powert = aircraft.powerfor
    aircraft.energy += aircraft.powert*dt

    time += dt
    aircraft.appendData(time,rho)
    aircraft.updateWeight(dt)

test_lift = True

#descent
min_height = 500. #m minimum height before decelrating from climb speed to 0 in forward speed
descent = aircraft.roc #m/s
final_deceleration = (0 - aircraft.convertToTAS(aircraft.climb_speed,density(min_height),rho0))/(min_height/descent) #deceleration from minheight to ground
aircraft.timestamps['descent'] = time
T_perp = 0
while aircraft.height>0:
    #while not at the ground yet

    rho = density(aircraft.height)
    aircraft.velocityIAS = aircraft.convertToIAS(aircraft.velocityTAS,rho,rho0)
    aircraft.height -= descent*dt
    
    if aircraft.velocityIAS > aircraft.climb_speed+0.1:
        #if not reached climb (descent) speed yet
        aircraft.acceleration = -climb_acceleration
        aircraft.timestamps['IAS Descent'] = time
       
    
    elif aircraft.height > min_height:
        #if reached climb speed but not at min_height
        temp_target_TAS = aircraft.convertToTAS(aircraft.climb_speed,rho,rho0)
        aircraft.acceleration = (temp_target_TAS-aircraft.velocityTAS)/dt
        aircraft.timestamps['Min Height'] = time

    else:
        #final part of descent
        aircraft.acceleration = final_deceleration


    aircraft.velocityTAS += aircraft.acceleration*dt
        
    aircraft.lift = aircraft.findForce(aircraft.cl,rho)
    T_perp = max(aircraft.weight - aircraft.lift, 0)

    # if test_lift and T_perp != 0.:
    #    test_lift = False
    #    aircraft.timestamps['lift lost'] = time


    aircraft.poweru = aircraft.inducedPower(T_perp,rho)

    if T_perp == 0: aircraft.cl = (aircraft.weight)/(aircraft.S*0.5*rho*aircraft.velocityTAS**2)

    
    cd = aircraft.dragCoefficient(aircraft.cl)
    aircraft.drag = aircraft.findForce(cd,rho)

    aircraft.powerfor = aircraft.velPower(aircraft.drag)/aircraft.eta_prop
    aircraft.powert = aircraft.poweru+aircraft.powerfor

    aircraft.energy += aircraft.powert*dt

    aircraft.xposition += aircraft.velocityTAS*dt
    aircraft.angle = aircraft.findAngle(T_perp,aircraft.drag+T_acc)

    time += dt
    aircraft.appendData(time,rho)
    aircraft.updateWeight(dt)

#second hover maneuver
count = 0 #temporary counter to time last hover period
aircraft.timestamps['hover 2'] = time

while count < hover_time:
        aircraft.poweru = aircraft.inducedPower(aircraft.weight,rho)
        aircraft.powert = aircraft.poweru
        aircraft.energy += aircraft.powert*dt
        count += dt
        time += dt
        aircraft.appendData(time,rho)
        aircraft.updateWeight(dt)





data = aircraft.data #point to new dictionary to avoid long coding
print(aircraft.energy)
print(max(data['powert']))
print(aircraft.r)
b = math.sqrt(aircraft.S*aircraft.A)
c = b/aircraft.A
print('b',b)
print('c',c)


# size the figure frame: the combination [10,7.5] and the fontsizes with the scaling '0.85\linewidth' in latex works best
# fig=plt.figure(figsize=(10,7.5)) #curve
fig=plt.figure(figsize=(10,7.5)) #phases

# initialize single subplot
ax = plt.subplot(111)

power_phases=pd.Series({'Hover': data['powert'][0],'Climb': max(data['powert']), 'Cruise':data['powert'][10000],'Descent':data['powert'][34000]})
print(power_phases)
# power_phases.plot(kind='bar')


x = 'time'
y = 'powert'    
end = 5700
time_plot = data[x][0:end]
y_plot = data[y][0:end]
plt.plot(time_plot,y_plot, 'r-')
plt.xlabel(r'Time $\left[s\right]$', fontsize=28)
ax.set_yticklabels(['{:,}'.format(int(x)) for x in ax.get_yticks().tolist()])
plt.ylabel(r'Power $\left[kW\right]$',fontsize=28)
# plt.xlabel(r'Flight Phase', fontsize=18)
countl = 0
for n in aircraft.timestamps.keys():
    i = aircraft.timestamps[n]
    if i < max(time_plot):
        margin = 7    
        countl += 1
        plt.axvline(x = i, color = 'b', linestyle = '--', linewidth = '1')
        index = data[x].index(i)
        if n == 'Cruise Altitude':
            margin = -90
        
        elif n == 'Begin Climb':
            margin = -70
        plt.text(i+margin,data[y][index]+20,n,fontsize = 18)
# plt.legend()
plt.xticks(fontsize=24)
plt.yticks(fontsize=24)
plt.rc('font', family='serif')
plt.grid(True)
plt.show()


# needed = ['time','weight','powert','energy']
# export = {}
# for i in range(len(needed)):
#     export[needed[i]] = data[needed[i]]
# scipy.io.savemat('bendata.mat', export)

# with open('tiltrotor.pickle', 'wb') as f:
#     pickle.dump(aircraft, f)

# pickle_in = open('tiltrotor.pickle', 'rb')
# clf = pickle.load(pickle_in)






