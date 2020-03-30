import math
import matplotlib.pyplot as plt
import numpy as np
# from matplotlib import style
# style.use('classic')

def density(height):
    temp = 288.15-0.0065*height
    pressure = 101325*math.pow(temp/288.15,(-9.80665/(-0.0065*287.05)))
    return pressure/(287.05*temp)

class WingEmbd():

    def __init__(self,cruise_altitude,cruise_speed,cruise_range,climb_speed,mass,g,roc,stall,A,e,M,r,dl,eta_prop):
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
        self.cl = 1.2 #Cl of aircraft intialised at 1.2
        self.cd0 = 0.017 #cd0 of aircraft assumed
        self.A = A #aspect ratio
        self.e = e #oswald efficiency
        self.S = mass*g/(self.cl*0.5*1.225*stall**2) #surface area of aircraft
        self.M = M #figure of merit
        self.r = r #radius of blades
        self.dl = dl #download factor
        self.eta_prop = eta_prop #propeller efficieny
        self.eta_duct = eta_duct #duct efficiency
        self.lift = 0 #aircraft lift
        self.drag = 0 #aircraft drag
        self.eta_engine = 0.32 #efficiency for energy transfer from fuel to shaft
        self.SPEC = 46.4e6 #specific energy consumption (46.6 for jetfuel)
        self.data = {'time':[], 'height':[], 'density':[], 'xposition':[], 'velocityTAS':[], 'velocityIAS':[], 'acceleration':[], 'powert':[], 'poweru': [], 'powerfor': [], 'poweracc': [], 'powerc':[], 'energy': [], 'angle': [], 'lift/drag': [], 'weight':[]} #dictionary of lists of data

    def findForce(self,coefficient,rho):
        # find aerodynamic force
        return 0.5*rho*self.S*coefficient*self.velocityTAS**2

    def inducedPowerhover(self,thrust,rho):
        #find induced power
        return (thrust*self.dl/self.M)*math.sqrt((-0.5*self.velocityTAS**2)+math.sqrt(((0.5*self.velocityTAS**2)**2)+self.vHover(rho)**4))/1000.

    def testinducedPower(self,thrust,v,rho):
        #ignore this
        return (thrust*self.dl/self.M)*math.sqrt((-0.5*v**2)+math.sqrt(((0.5*v**2)**2)+self.vHover(rho)**4))/1000.

    def vHover(self,rho):
        #velocity induced in hover
        return math.sqrt((self.dl*self.mass*self.g)/(4*rho*math.pi*self.r**2))

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



    def lift(self, MTOW, rho, V, S, S_o, C_lmax):
        L_mc = 0.5 * rho * V**2 * S * C_lmax             # max lift produced by closed wing
        L_mo = 0.5 * rho * V**2 * S_o * C_lmax           # max lift produced by open wing
        # L_w = np.where(L_mc > MTOW, MTOW, L_mo)
        # C_l = np.where(L_mc > MTOW, 2 * MTOW / (rho * V**2 * S), C_lmax)
        # S_w   = np.where(L_mc > MTOW, S, S_o)
        #
        # return L_w, C_l, S_w
        return L_mo, C_lmax, S_o

    def c_d(self, C_l, C_d0, A, e):
        return C_d0 + C_l**2 / (pi * A * e)

    def drag(self, rho, V, S, C_d, ):
        return  0.5 * rho * V**2 * S * C_d


    def c_l_blade(self, L_rotor, N, B, rho, omega, c, R):
        return 6*L_rotor / (N*B*rho*omega**2*c*R**3)

    def collective(self, L_rotor, N, B, rho, omega, c, R, theta_t, v_i):
        num_1 = L_rotor / (N * B * rho * omega**2 * pi * c)
        num_2 = R**2 * (theta_t * R - v_i / omega) / 2
        col   = (num_1 - num_2) * R**3 / 3
        return np.where(v_i != 0, col, 0)

    def v_induced(self, L_rotor, rho, V, c, R, C_dr, N):
        A_p = N * c * R
        A = N * pi * R**2
        Thrust = (1 + (C_dr * A_p)/A) * L_rotor
        v_induced = np.sqrt(-((V**2)/2) + np.sqrt((V**4)/4 + (Thrust/(2 * rho * A))**2))
        return np.where(Thrust != 0, v_induced, 0), Thrust

    def aoa_cuad(self, theta_t, R, v_i, omega, alfa_c):
        return theta_t * R - v_i / omega + alfa_c

    def aoa_lin(self, rho, R_o, R_f, omega, c, v_i, N, B, L_rotor, theta_t, theta_r):
        return (L_rotor / (N * B * rho * omega**2 * c * pi) - theta_r * (R_f**3/3 - R_o**3/3) - theta_t * (R_f**3/4 -R_o**3/4) + v_i*(R_f**2/(2*omega)- R_o**2/(2*omega))) / (R_f**3/3 - R_o**3/3)

    def omegas_cuad(self, rho, R, c, v_i, N, B, L_rotor, theta_t):
        A = 0.5 * rho * R**3 * c * pi * theta_t * N * B
        B = -0.5 * rho * R**2 * c * pi * v_i * N * B
        C = -L_rotor
        omega_1 = (-B - np.sqrt(B**2 - 4 * A * C)) / (2 * A)
        omega_2 = (-B + np.sqrt(B**2 - 4 * A * C)) / (2 * A)
        return omega_1, omega_2

    def omegas_lin(self, rho, R, c, v_i, N, B, L_rotor, theta_t, theta_r):
        A = N * B * rho * c * pi * (R**3 * theta_r / 3 + R**3 * theta_t / 4)
        B = N * B * rho * c * pi * (R**2 * v_i/ 2)
        C = -L_rotor
        omega_1 = (-B - np.sqrt(B**2 - 4 * A * C)) / (2 * A)
        omega_2 = (-B + np.sqrt(B**2 - 4 * A * C)) / (2 * A)
        return omega_1, omega_2

    def c_d_blade(self, alfa_c, omega, R, v_i, theta_t):
        x = theta_t * R - v_i / omega + alfa_c
        c_d = 20441 * x**6 -1505.9 * x**5 + 442.64 *x**4 - 63.716 * x**3 + 4.7524 * x**2 - 0.0961 * x + 0.0058
        return np.where(v_i != 0, c_d, 0)

    def c_d_blade_1(self, x):
        c_d = 20441 * x**6 -1505.9 * x**5 + 442.64 *x**4 - 63.716 * x**3 + 4.7524 * x**2 - 0.0961 * x + 0.0058
        return c_d

    def p_parasite(self, D, V, eta_prop):
        p_p = D * V/eta_prop
        return p_p

    def p_induced(self, Thrust, v_induced, FM, eta_tr, eta_inst, eta_duct):
        return Thrust*v_induced/FM * (1/(eta_tr*eta_inst))*eta_duct

    def p_profile(self, N, B, c_d_blade, rho, omega, R, c):
        return np.where(c_d_blade != 0.0058, N * B * c_d_blade * rho * omega**3 * R**4  * c /8, 0)

    def p_climb(self, MTOW, RoC, P_r):
        return MTOW*RoC

    def l_check(self, rho, omega, R, c, N, B, theta_t, v_i):
        return 0.5 * rho * omega**2 * R**2 * c * pi * N * B * (theta_t * R - v_i/omega)


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
        self.data['weight'].append(self.weight)
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
cruise_speed = 350
cruise_range = 300
climb_speed = 220
mass = 3600
g = 9.80665
roc = 8
stall = 45
A = 5
e = 0.8
M = 0.65
r = 2.2
dl = 1
eta_prop = 0.75
eta_duct = 0.6





aircraft = WingEmbd(cruise_altitude,cruise_speed,cruise_range,climb_speed,mass,g,roc,stall,A,e,M,r,dl,eta_prop, eta_duct)

climb_acceleration = aircraft.climb_speed/80 #time to reach climbspeed (divide by seconds required.)

time = 0 #s #intitalise time
dt = 0.1 #s #time step

#climb
while aircraft.height < aircraft.cruise_altitude:

    rho = density(aircraft.height) #find density


    #hover phase
    while time < hover_time:
        aircraft.poweru = aircraft.inducedPower(aircraft.mass*aircraft.g,rho)
        aircraft.powert = aircraft.poweru
        aircraft.energy += aircraft.powert*dt
        time += dt
        aircraft.appendData(time,rho)
    aircraft.updateWeight(dt)






    aircraft.height += aircraft.roc*dt

    aircraft.velocityIAS = aircraft.convertToIAS(aircraft.velocityTAS,rho,rho0)

    #check if climb velocity is reached
    if (aircraft.velocityIAS < aircraft.climb_speed) and IASswitch:
        aircraft.acceleration = climb_acceleration
        T_acc = aircraft.acceleration*aircraft.mass #thrust required for acceleration

    else:
        IASswitch = False #set thtat climb speed has been reached
        temp_target_TAS = aircraft.convertToTAS(aircraft.climb_speed,rho,rho0) #set target speed in TAS
        aircraft.acceleration = (temp_target_TAS-aircraft.velocityTAS)/dt
        T_acc = aircraft.acceleration*aircraft.mass


    aircraft.velocityTAS += aircraft.acceleration*dt #update velocity

    aircraft.lift = aircraft.findForce(aircraft.cl,rho) #find lift
    cd = aircraft.dragCoefficient(aircraft.cl)
    aircraft.drag = aircraft.findForce(cd,rho)

    T_perp = max(aircraft.mass*aircraft.g - aircraft.lift,0) #either set required perpendicular thrust to counteract weight or zero if lift covers it
    aircraft.poweru = aircraft.inducedPower(T_perp,rho)
    aircraft.powerfor = aircraft.velPower(aircraft.drag)/aircraft.eta_prop

    if T_perp == 0: aircraft.cl = (aircraft.weight)/(aircraft.S*0.5*rho*aircraft.velocityTAS**2)


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


#cruise
time_cruise = time
cruise_acceleration = (aircraft.cruise_speed-aircraft.convertToTAS(aircraft.climb_speed,rho,rho0))/120 #acceleration to reach cruise speed in given seconds (divisor)

while aircraft.xposition<aircraft.cruise_range*1000:

    while aircraft.velocityTAS < aircraft.cruise_speed:
        #accelerate to crusie speed
        aircraft.cl = (aircraft.mass*aircraft.g)/(aircraft.S*0.5*rho*aircraft.velocityTAS**2)
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



#descent
min_height = 500. #m minimum height before decelrating from climb speed to 0 in forward speed
final_deceleration = (0 - aircraft.convertToTAS(aircraft.climb_speed,density(min_height),rho0))/(min_height/aircraft.roc) #deceleration from minheight to ground

while aircraft.height>0:
    #while not at the ground yet

    rho = density(aircraft.height)
    aircraft.velocityIAS = aircraft.convertToIAS(aircraft.velocityTAS,rho,rho0)
    aircraft.height -= aircraft.roc*dt

    if aircraft.velocityIAS > aircraft.climb_speed+0.1:
        #if not reached climb (descent) speed yet
        aircraft.acceleration = -climb_acceleration


    elif aircraft.height > min_height:
        #if reached climb speed but not at min_height
        temp_target_TAS = aircraft.convertToTAS(aircraft.climb_speed,rho,rho0)
        aircraft.acceleration = (temp_target_TAS-aircraft.velocityTAS)/dt

    else:
        #final part of descent
        aircraft.acceleration = final_deceleration


    aircraft.velocityTAS += aircraft.acceleration*dt


    T_perp = max(aircraft.mass*aircraft.g - aircraft.lift,0)
    aircraft.poweru = aircraft.inducedPower(T_perp,rho)

    if T_perp == 0: aircraft.cl = (aircraft.weight)/(aircraft.S*0.5*rho*aircraft.velocityTAS**2)

    aircraft.lift = aircraft.findForce(aircraft.cl,rho)
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

while count < hover_time:
        aircraft.poweru = aircraft.inducedPower(aircraft.mass*aircraft.g,rho)
        aircraft.powert = aircraft.poweru
        aircraft.energy += aircraft.powert*dt
        count += dt
        time += dt
        aircraft.appendData(time,rho)
        aircraft.updateWeight(dt)


print(aircraft.energy)
data = aircraft.data #point to new list to avoid long coding
x = 'time'
y = 'powert'
plt.plot(data[x],data[y])
plt.xlabel(x)
plt.ylabel(y)
# plt.axvline(x = 1000, color ='r')
# plt.legend()
plt.show()

