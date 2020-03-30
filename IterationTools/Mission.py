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

def temperature(height):
    return 288.15-0.0065*height
class TiltRotor():

    def __init__(self,rotor,cruise_altitude,cruise_speed,cruise_range,climb_speed,mass,g,roc,stall,A,S,e,diskloading,dl,eta_prop):
        #aircraft characteristics
        self.props = rotor
        self.mass = mass #kg
        self.cl = 1.2 #Cl of self intialised at 1.2 (cl max())
        self.c = 0.2358
        self.S = S#mass*g/(self.cl*0.5*1.225*stall**2) #surface area of self
        self.cwing = math.sqrt(A*S)/A
        self.M = 1 #figure of merit
        self.N = 2 #number of rotors
        self.B = 3 #number of blades
        self.diskloading = diskloading
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
        self.data = {'time':[], 'height':[], 'density':[], 'xposition':[], 'velocityTAS':[], 'velocityIAS':[], 'acceleration':[], 'powert':[], 'poweru': [], 'powerfor': [], 'poweracc': [], 'powerc':[], 'energy': [], 'angle': [], 'lift/drag': [], 'weight':[], 'vinduced':[], 'temp':[], 'drag': [], 'vperp': [], 'thrusta': [], 'powera': [], 'check': []} #dictionary of lists of data
        self.timestamps = {}
        self.previous = 1000
        self.vperp = 0
        self.thrusta = 0
        self.powera = 0
        self.check = 0
        self.alexoutput = []
    
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

    def inducedAlex(self,thrust,roc):
        self.check = 0
        vperp = self.velocityTAS*math.cos(self.angle)+roc*math.cos((math.pi/2)-self.angle)
        self.vperp = vperp
        self.thrusta = thrust
        # print('int(vperp)', int(vperp))
        # print('thrust', int(round(thrust/2 / 200.0)) * 200)
        power = self.props[int(vperp)][int(round(thrust/2 / 200)) * 200]
        power = power*(1+(self.dl-1)*math.sin(self.angle))
        if power > 10**8:
            power = self.velPower(thrust)*500/self.eta_prop
            self.check = 1
        elif np.isnan(power):
            self.check = 2
            power = self.previous
        self.powera = 2*power/1000
        
        return 2*power/1000

    def roundtime(self):
        self.data['time'] = list(np.round(self.data['time'],1))
        return

    def valuesForAlex(self):
        times = [0.1,self.timestamps['Begin Climb'],self.timestamps['Begin Climb']+0.2, self.timestamps['IAS Climb Velocity']+1, self.timestamps['Cruise Altitude'], self.timestamps['Cruise Altitude']+0.5, self.timestamps['cruise speed'], self.timestamps['cruise speed']+0.3, self.timestamps['descent'], self.timestamps['IAS Descent']+0.5, self.timestamps['Min Height'], self.timestamps['hover 2'], self.data['time'][-1]]
        index = []
        rl = []
        check = 0
        for i in times:
            index.append(self.data['time'].index(round(i,1)))
        for i in index:
            temp = []
            temp.append(self.data['thrusta'][i])
            temp.append(self.data['vperp'][i])
            temp.append(self.data['height'][i])
            temp.append(self.data['angle'][i])
            temp.append(self.data['time'][i])
            rl.append(temp)
            if self.data['thrusta'][i] > check: check = self.data['thrusta'][i]
        # print('max thrust a actual', max(self.data['thrusta']))
        # print('max thrust a alex adata', check)
        return rl
        
    def testInducedPower(self,thrust,rho): 
        self.eta_trans = 1#0.96
        self.eta_inst = 0.98
        #ignore this
        return (thrust*self.dl/self.M)*math.sqrt((-0.5*self.velocityTAS**2)+math.sqrt(((0.5*self.velocityTAS**2)**2)+self.vHover(thrust,rho)**4))*(1/(self.eta_trans*self.eta_inst))/1000.

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
        self.data['temp'].append(temperature(self.height))
        self.data['drag'].append(self.drag)
        self.data['vperp'].append(self.vperp)
        self.data['thrusta'].append(self.thrusta)
        self.data['powera'].append(self.powera)
        self.data['check'].append(self.check)
        self.previous = self.powert
        return

    def updateWeight(self,dt):
        #decreases weight by fuel used
        self.weight -= 0#(self.powert*self.g*dt/(self.eta_engine*(self.SPEC/1000)))
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
                inducedpp = self.testInducedPower(self.weight,rho)
                #===================================================
                self.poweru = self.inducedAlex(self.weight, 0)
                #===================================================
                self.powert = self.poweru
                self.energy += self.powert*dt
                time += dt
                self.appendData(time,rho)
                self.updateWeight(dt)
                self.timestamps['Begin Climb'] = time
            
            self.height += self.roc*dt

            self.velocityIAS = self.convertToIAS(self.velocityTAS,rho,rho0)

            #check if climb velocity is reached
            if (self.velocityIAS < self.climb_speed) and IASswitch:
                self.acceleration = climb_acceleration
                T_acc = self.acceleration*self.mass #thrust required for acceleration
                self.timestamps['IAS Climb Velocity'] = time

            else:
                IASswitch = False #set that climb speed has been reached
                temp_target_TAS = self.convertToTAS(self.climb_speed,rho,rho0) #set target speed in TAS
                self.acceleration = (temp_target_TAS-self.velocityTAS)/dt
                T_acc = self.acceleration*self.mass


            self.velocityTAS += self.acceleration*dt #update velocity

            self.lift = self.findForce(self.cl,rho) #find lift 
            cd = self.dragCoefficient(self.cl)
            self.drag = self.findForce(cd,rho)

            T_perp = max(self.weight - self.lift,0.) #either set required perpendicular thrust to counteract weight or zero if lift covers it
            self.T_perp = T_perp
            if test_lift and T_perp == 0.:
                test_lift = False
                self.timestamps['Full Wing Flight'] = time

            #===================================================
            # self.poweru = self.inducedPower(T_perp,rho) 
            # self.powerfor = self.velPower(self.drag)/self.eta_prop
            thrust = math.sqrt(math.pow(T_perp,2)+math.pow(self.drag,2))
            self.poweru = self.inducedAlex(thrust+T_acc, self.roc)
            self.powerfor = 0
            #===================================================

            if T_perp == 0.: self.cl = (self.weight)/(self.S*0.5*rho*self.velocityTAS**2)

            
            self.poweracc = self.velPower(T_acc)
            self.powerc = self.roc*(self.lift + T_perp/2)/1000  #power to climb is a component of both lift due to wing and perpendicular thrust
            self.powert = self.poweru+self.powerfor+self.poweracc+self.powerc
            self.energy += self.powert*dt

            self.xposition += self.velocityTAS*dt
            self.angle = self.findAngle(T_perp,self.drag+T_acc)

            time += dt
            self.appendData(time,rho)
            self.updateWeight(dt)


        #set unused variables for next phase to zero
        self.poweru = 0
        self.powerc = 0
        self.roc = 0

        self.timestamps['Cruise Altitude'] = time

        #cruise
        time_cruise = time
        cruise_acceleration = (self.cruise_speed-self.convertToTAS(self.climb_speed,rho,rho0))/120 #acceleration to reach cruise speed in given seconds (divisor)

        while self.xposition<self.cruise_range*1000:

            while self.velocityTAS < self.cruise_speed-0.1:
                #accelerate to crusie speed
                self.cl = (self.weight)/(self.S*0.5*rho*self.velocityTAS**2)
                self.lift = self.findForce(self.cl,rho)
                cd = self.dragCoefficient(self.cl)
                self.drag = self.findForce(cd,rho)
                self.acceleration = min(cruise_acceleration,(self.cruise_speed-self.velocityTAS)/dt) #take minimum of cruise acceleration 
                T_acc = self.acceleration*self.mass
                self.velocityTAS += self.acceleration*dt
                self.velocityIAS = self.convertToIAS(self.velocityTAS,rho,rho0)

                self.xposition += self.velocityTAS*dt
                # self.powerfor = self.velPower(self.drag)/self.eta_prop
                self.powerfor = self.inducedAlex(self.drag+T_acc, self.roc)
                self.poweracc = self.velPower(T_acc)
                self.powert = self.powerfor+self.poweracc
                self.energy += self.powert*dt

                time += dt
                self.appendData(time,rho)
                self.updateWeight(dt)
                self.timestamps['cruise speed'] = time
            
            self.cl = (self.weight)/(self.S*0.5*rho*self.velocityTAS**2)
            self.lift = self.findForce(self.cl,rho)
            cd = self.dragCoefficient(self.cl)
            self.drag = self.findForce(cd,rho)
            self.acceleration = 0
            self.xposition += self.velocityTAS*dt
            self.poweracc = 0
            # self.powerfor = self.velPower(self.drag)/self.eta_prop
            self.powerfor = self.inducedAlex(self.drag,self.roc)

            self.powert = self.powerfor
            self.energy += self.powert*dt

            time += dt
            self.appendData(time,rho)
            self.updateWeight(dt)

        test_lift = True

        #descent
        min_height = 500. #m minimum height before decelrating from climb speed to 0 in forward speed
        descent = 8 #m/s #deceleration from minheight to ground
        self.timestamps['descent'] = time
        T_perp = 0
        while self.height>0:
            #while not at the ground yet

            rho = density(self.height)
            self.velocityIAS = self.convertToIAS(self.velocityTAS,rho,rho0)
            self.height -= descent*dt
            
            if self.velocityIAS > self.climb_speed+0.1:
                #if not reached climb (descent) speed yet
                self.acceleration = -climb_acceleration
                self.timestamps['IAS Descent'] = time
            
            
            elif self.height > min_height:
                #if reached climb speed but not at min_height
                temp_target_TAS = self.convertToTAS(self.climb_speed,rho,rho0)
                self.acceleration = (temp_target_TAS-self.velocityTAS)/dt
                self.timestamps['Min Height'] = time

            else:
                #final part of descent
                descent = 1.67
                final_deceleration = (0 - self.convertToTAS(self.climb_speed,density(min_height),rho0))/(min_height/descent)
                self.acceleration = final_deceleration
                


            self.velocityTAS += self.acceleration*dt
                
            self.lift = self.findForce(self.cl,rho)
            T_perp = max(self.weight - self.lift, 0)
            self.powerc = -descent*(T_perp/2)/1000

            # if test_lift and T_perp != 0.:
            #    test_lift = False
            #    self.timestamps['lift lost'] = time

            if T_perp == 0: self.cl = (self.weight)/(self.S*0.5*rho*self.velocityTAS**2)

            
            cd = self.dragCoefficient(self.cl)
            self.drag = self.findForce(cd,rho)


            # self.powerfor = self.velPower(self.drag)/self.eta_prop
            thrust = math.sqrt(math.pow(0*T_perp,2)+math.pow(self.drag,2))
            self.poweru = self.inducedAlex(thrust, 0)
            self.powerfor = 0
            self.powert = self.poweru+self.powerfor+self.powerc

            self.energy += self.powert*dt

            self.xposition += self.velocityTAS*dt
            self.angle = self.findAngle(T_perp,self.drag+T_acc)

            time += dt
            self.appendData(time,rho)
            self.updateWeight(dt)

        #second hover maneuver
        count = 0 #temporary counter to time last hover period
        self.timestamps['hover 2'] = time

        while count < hover_time:
                self.poweru = self.inducedAlex(self.weight,0)
                self.powert = self.poweru
                self.energy += self.powert*dt
                count += dt
                time += dt
                self.appendData(time,rho)
                self.updateWeight(dt)

        return

    def plotFlight(self,stringx,stringy):
        print('energy',self.energy)
        print('max power', max(self.data['powert']))
        b = math.sqrt(self.S*self.A)
        c = b/self.A
        print('r',self.r)
        print('b',b)
        print('c',c)
        print('S', self.S)
        # size the figure frame: the combination [10,7.5] and the fontsizes with the scaling '0.85\linewidth' in latex works best
        # fig=plt.figure(figsize=(10,7.5)) #curve
        fig=plt.figure(figsize=(10,7.5)) #phases

        # initialize single subplot
        ax = plt.subplot(111)

        power_phases=pd.Series({'Hover': self.data['powert'][0],'Climb': max(self.data['powert']), 'Cruise':self.data['powert'][10000],'Descent':self.data['powert'][34000]})
        print(power_phases)
        # power_phases.plot(kind='bar')
    
        end = len(self.data[stringx])#5700
        time_plot = self.data[stringx][0:end]
        y_plot = self.data[stringy][0:end]
        plt.plot(time_plot,y_plot, 'r-')
        plt.xlabel(r'Time $\left[s\right]$', fontsize=28)
        # ax.set_yticklabels(['{:,}'.format(int(x)) for x in ax.get_yticks().tolist()])
        plt.ylabel(r'Power $\left[kW\right]$',fontsize=28)
        # plt.xlabel(r'Flight Phase', fontsize=18)
        countl = 0
        if stringx == 'time':
            for n in self.timestamps.keys():
                i = self.timestamps[n]
                if i < max(time_plot):
                    margin = 7    
                    countl += 1
                    plt.axvline(x = i, color = 'b')#, linestyle = '--', linewidth = '1')
                    index = self.data[stringx].index(i)
                    if n == 'Cruise Altitude':
                        margin = -90
                    
                    elif n == 'Begin Climb':
                        margin = -70
                    plt.text(i+margin,self.data[stringy][index]+20,n,fontsize = 18)
        # plt.legend()
        plt.xticks(fontsize=24)
        plt.yticks(fontsize=24)
        plt.rc('font', family='serif')
        plt.grid(True)
        plt.show()


def main():
    pass

if __name__ == "__main__":
    cruise_altitude = 2000 #m


    aircraft = TiltRotor(cruise_altitude,350,270,220,4000.,9.80665,8,45,6.17,0.9,60,1.15,0.75) #create tiltRotor object
    aircraft.performFlight(1.225) #simulate the mission
    aircraft.plotFlight('time','powert') #plot mission profile for x and y value
    print(aircraft.timestamps)



# needed = ['time','weight','powert','energy']
# export = {}
# for i in range(len(needed)):
#     export[needed[i]] = data[needed[i]]
# scipy.io.savemat('bendata.mat', export)

# with open('tiltrotor.pickle', 'wb') as f:
#     pickle.dump(aircraft, f)

# pickle_in = open('tiltrotor.pickle', 'rb')
# clf = pickle.load(pickle_in)
