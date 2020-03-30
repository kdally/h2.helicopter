import math
import numpy as np

MTOM = 4000
MTOW = MTOM*9.81
dl = 1.15
M = 1
velocityTAS = 0
N = 2
B = 3
r = 3.15
c = 0.2358
omega = 85




def density(height):
    temp = 288.15-0.0065*height
    pressure = 101325*math.pow(temp/288.15,(-9.80665/(-0.0065*287.05)))
    return pressure/(287.05*temp)

class TiltRotor():

    def __init__(self,cruise_altitude,cruise_speed,cruise_range,climb_speed,mass,g,roc,stall,A,e,M,dl,eta_prop):
        #aircraft characteristics
        self.mass = mass #kg
        self.cl = 1.2 #Cl of self intialised at 1.2 (cl max())
        self.c = 0.2358
        self.S = mass*g/(self.cl*0.5*1.225*stall**2) #surface area of self
        self.M = M #figure of merit
        self.N = 2 #number of rotors
        self.B = 3 #number of blades
        self.diskloading = 60
        self.r = 0.5*math.sqrt((self.mass/(self.diskloading*2))*4/math.pi) #radius of blades
        self.eta_prop = eta_prop #propeller efficieny

        #for performFlight()
        self.cruise_altitude =  cruise_altitude #m
        self.cruise_speed = cruise_speed/3.6 #km/h
        self.climb_speed = climb_speed/3.6 #km/h speed forward in climb, also assumed in descent
        self.cruise_range = cruise_range #km
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
        self.cd0 = 0.017 #cd0 of self assumed
        self.A = A #aspect ratio
        self.e = e #oswald efficiency
        self.omega = 85 #RPS rad/s
        self.dl = dl #download factor
        self.lift = 0 #self lift
        self.drag = 0 #self drag
        self.eta_engine = 0.4 #efficiency for energy transfer from fuel to shaft
        self.SPEC = 141.8e6 #specific energy consumption (46.6 for jetfuel)
        self.vinduced = 0
        self.data = {'time':[], 'height':[], 'density':[], 'xposition':[], 'velocityTAS':[], 'velocityIAS':[], 'acceleration':[], 'powert':[], 'poweru': [], 'powerfor': [], 'poweracc': [], 'powerc':[], 'energy': [], 'angle': [], 'lift/drag': [], 'weight':[], 'vinduced':[]} #dictionary of lists of data
        self.timestamps = {}

    def c_l_blade(self, L_rotor, rho):
        return 6*L_rotor / (self.N*self.B*rho*self.omega**2*self.c*self.r**3)

    def c_d_blade_1(self,x):
        return 20441 * x**6 -1505.9 * x**5 + 442.64 *x**4 - 63.716 * x**3 + 4.7524 * x**2 - 0.0961 * x + 0.0058

    def p_profile(self, c_d_blade, rho):
        return np.where(c_d_blade != 0.0058, self.N * self.B * c_d_blade * rho * self.omega**3 * self.r**4  * self.c /8, 0)

    def inducedPower(self,thrust,rho):
        #find induced power
        self.eta_trans = 1#0.96
        self.eta_inst = 0.98
        induced = (thrust*self.dl/self.M)*math.sqrt((-0.5*self.velocityTAS**2)+math.sqrt(((0.5*self.velocityTAS**2)**2)+self.vHover(thrust,rho)**4))*(1/(self.eta_trans*self.eta_inst))/1000.
        self.vinduced = math.sqrt((-0.5*self.velocityTAS**2)+math.sqrt(((0.5*self.velocityTAS**2)**2)+self.vHover(thrust,rho)**4))*(1/(self.eta_trans*self.eta_inst))
        C_l_blade = self.c_l_blade(thrust,rho)
        aoaob = C_l_blade/(2*math.pi)
        C_d_blade = self.c_d_blade_1(aoaob)
        P_profile = self.p_profile(C_d_blade,rho)

        return induced+P_profile/1000.
        # return induced

    def vHover(self,thrust,rho):
        #velocity induced in hover
        return math.sqrt((self.dl*thrust)/(2*self.N*rho*math.pi*self.r**2))


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

    def performFlight(self,rho0):

        IASswitch = True #constant to check if climb speed reached

        hover_time = 60 #s

        climb_acceleration = self.climb_speed/80 #time to reach climbspeed (divide by seconds required.)

        time = 0 #s #intitalise time
        dt = 0.1 #s #time step
        test_lift = True

        #climb
        while self.height < self.cruise_altitude:

            rho = density(self.height) #find density

            #hover phase
            while time < hover_time:
                self.poweru = self.inducedPower(self.weight,rho)
                self.powert = self.poweru
                self.energy += self.powert*dt
                time += dt
                self.appendData(time,rho)
                self.updateWeight(dt)
                self.timestamps['Begin Climb'] = time

            self.height += self.roc*dt

        return

def vHover(thrust,rho):
    #velocity induced in hover
    return math.sqrt((dl*thrust)/(2*N*rho*math.pi*r**2))

def c_l_blade(L_rotor, rho):
    return 6*L_rotor / (N*B*rho*omega**2*c*r**3)

def c_d_blade_1(x):
    return 20441 * x**6 -1505.9 * x**5 + 442.64 *x**4 - 63.716 * x**3 + 4.7524 * x**2 - 0.0961 * x + 0.0058

def p_profile(c_d_blade, rho):
    return np.where(c_d_blade != 0.0058, N * B * c_d_blade * rho * omega**3 * r**4  * c /8, 0)


def inducedPower(thrust,rho):
    #find induced power
    eta_trans = 1#0.96
    eta_inst = 0.98
    induced = (thrust*dl/M)*math.sqrt((-0.5*velocityTAS**2)+math.sqrt(((0.5*velocityTAS**2)**2)+vHover(thrust,rho)**4))*(1/(eta_trans*eta_inst))/1000.
    C_l_blade = c_l_blade(thrust,rho)
    aoaob = C_l_blade/(2*math.pi)
    C_d_blade = c_d_blade_1(aoaob)
    P_profile = p_profile(C_d_blade,rho)

    return induced+P_profile/1000.




cruise_altitude = 2000 #m


print(inducedPower(MTOW, 1.225))

aircraft = TiltRotor(cruise_altitude,350,270,220,4000.,9.80665,8,45,6.17,0.9,1,1.15,0.75)
print(aircraft.inducedPower(aircraft.weight, 1.225))


aircraft.performFlight(1.225)
print(aircraft.data['powert'][0])
print(aircraft.S, aircraft.A)

