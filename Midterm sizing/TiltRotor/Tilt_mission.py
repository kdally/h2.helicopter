
import math
import matplotlib.pyplot as plt

def density(height):
    temp = 288.15-0.0065*height
    pressure = 101325*math.pow(temp/288.15,(-9.80665/(-0.0065*287.05)))
    return pressure/(287.05*temp)

class TiltRotor():

    def __init__(self,cruise_power,hover_power,max_power,cruise_altitude,cruise_speed,climb_speed,mass,g,roc,S,A,e,M,r):
        self.cruise_power = cruise_power #kW
        self.hover_power = hover_power #kW
        self.max_power =  max_power#kW
        self.cruise_altitude =  cruise_altitude #m
        self.cruise_speed = cruise_speed/3.6 #km/h
        self.climb_speed = climb_speed/3.6
        self.mass = mass #kg
        self.g = g #m/s
        self.roc = roc #m/s
        self.height = 0
        self.energy = 0        
        self.velocity = 0
        self.acceleration = 0
        self.xposition = 0
        self.powert = 0
        self.poweru = 0
        self.powerfor = 0
        self.poweracc = 0
        self.powerc = 0
        self.angle = math.pi/2.
        self.cl = 1
        self.A = A
        self.e = e
        self.S = S
        self.M = M
        self.r = r
        self.data = {'time':[], 'height':[], 'density':[], 'xposition':[], 'velocity':[], 'acceleration':[], 'powert':[], 'poweru': [], 'powerfor': [], 'poweracc': [], 'powerc':[], 'energy': [], 'angle': []}
    
    def findForce(self,coefficient,rho,velocity,S):
        return 0.5*rho*S*coefficient*velocity**2
    
    def appendData(self,time,density):
        self.data['time'].append(time)
        self.data['height'].append(self.height)
        self.data['density'].append(density)
        self.data['xposition'].append(self.xposition)
        self.data['acceleration'].append(self.acceleration)
        self.data['velocity'].append(self.velocity)
        self.data['powert'].append(self.powert)
        self.data['poweru'].append(self.poweru)
        self.data['powerfor'].append(self.powerfor)
        self.data['poweracc'].append(self.poweracc)
        self.data['powerc'].append(self.powerc)
        self.data['energy'].append(self.energy)
        self.data['angle'].append(self.angle)
        return

    


#mission profile
hover_time = 60 #s
cruise_altitude = 2000 #m

aircraft = TiltRotor(379.88,1275.977148,1364.27,cruise_altitude,350,220,4000.,9.80665,7,24.38,10,0.9,0.65,2.72)


time = 0 #s
dt = 1 #s
rho = 1.225

while aircraft.height < aircraft.cruise_altitude:
    #hover phase
    while time < hover_time:
        time += dt
        aircraft.poweru = (aircraft.mass*aircraft.g/aircraft.M)*math.sqrt(aircraft.mass*aircraft.g/(4*rho*math.pi*(aircraft.r**2)))/1000
        aircraft.powert = aircraft.poweru
        aircraft.energy += aircraft.powert*dt
        aircraft.appendData(time,rho)
    
    time += dt
    aircraft.height += aircraft.roc*dt
    rho = density(aircraft.height)

    lift = aircraft.findForce(aircraft.cl,rho,aircraft.velocity,aircraft.S)
    cd = (aircraft.cl**2)/(math.pi*aircraft.A*aircraft.e) + 0.002
    drag = aircraft.findForce(cd,rho,aircraft.velocity,aircraft.S)

    if aircraft.velocity < aircraft.climb_speed:
        T_acc = aircraft.mass*aircraft.roc*aircraft.cruise_speed/aircraft.cruise_altitude
    else:
        T_acc = 0
    aircraft.acceleration = T_acc/aircraft.mass

    T_perp = max(aircraft.mass*aircraft.g - lift,0)
    aircraft.poweru = (T_perp/aircraft.M)*math.sqrt(T_perp/(4*rho*math.pi*(aircraft.r**2)))/1000
    aircraft.powerfor = drag*aircraft.velocity/1000
    
    aircraft.poweracc = T_acc*aircraft.velocity/1000
    aircraft.powerc = aircraft.roc*(lift + T_perp/2)/1000
    aircraft.powert = aircraft.poweru+aircraft.powerfor+aircraft.poweracc+aircraft.powerc
    aircraft.energy += aircraft.powert*dt

    aircraft.velocity += aircraft.acceleration*dt
    aircraft.xposition += aircraft.velocity*dt

    aircraft.appendData(time,rho)

data = aircraft.data
plt.plot(data['time'],data['powert'])
plt.show()

    





