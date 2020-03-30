##
import numpy as np

def libfuelcell():
    Data = np.array([49, 63, 81, 98, 125])
    weight = np.array([21, 25, 29, 34, 43])
    name = np.array(["Power Cell S3 167", "Power Cell S3 215", "Power Cell S3 275", "Power Cell S3 335", "Power Cell S3 455"])
    return Data, weight, name

def libbattery():
    Data = np.array([49, 63, 81, 98, 125])
    weight = np.array([21, 25, 29, 34, 43])
    name = np.array(["Power Cell S3 167", "Power Cell S3 215", "Power Cell S3 275", "Power Cell S3 335", "Power Cell S3 455"])
    return Data, weight, name

def libelectricmotor():
    Data = np.array([160, 200, 60, 82, 87, 150, 150, 0.2, 200, 200, 260])
    weight = np.array([28.2, 37, 13.6, 33.5, 43, 33.5, 43, 0.3, 19.9, 16, 50])
    name = np.array(["YASA P400 RS", "YASA 750 R", "DHX FALCON", "HVH250 Standard", "HVH250 HT High Torque", "HVH250 High Flow Cooling", "HVH250 HT High Flow Cooling",
                     "Maxon EC-4pole", "EMRAX 268", "Magnax AXF225", "Siemens SP260D"])
    return Data, weight, name

def libmotorcontroller():
    Data = np.array([100, 220])
    weight = np.array([3.5, 9.1])
    name = np.array(["Fraunhofer inverter", "ROHM"])
    return Data, weight, name

def libdcdcconverter():
    Data = np.array([200])
    weight = np.array([3.2])
    name = np.array(["Fraunhofer IISB"])
    return Data, weight, name

def libbattery():
    Data1 = np.array([0.12, 0.141, 0.142, 0.448])
    Data2 = np.array([3, 2.84, 30.13, 4.48])
    name = np.array(["LIP", "Anto 1", "Anto 2", "NAC//Si/C"])
    return Data1, Data2, name
def libcompressor():
    Data = np.array([0.1, 0.016, 0.025, 0.005, 0.055, 0.017, 0.024, 0.22, 0.39, 0.63, 0.84])
    weight = np.array([4.2, 0.6, 0.6, 0.11, 0.7, 1.5, 1.5, 2.9, 5.1, 6.0, 6.4])
    name = np.array(["Aeristech","Celeroton CT17-700","Celeroton CT17-1000","Celeroton CT-15-150","Celeroton CT-14-1000","Celeroton CT 17.700GB","Celeroton CT-17.1000GB","Rotrex C15-60","Rotrex C30-54","Rotrex C38-91/92","Rotrex C38R-112"])
    return Data, weight, name

##
# from IterationTools.fpp import Componentspicker
# from IterationTools import Mission
#
# cruise_altitude = 4000 #m
# aircraft = Mission.TiltRotor(cruise_altitude, 350, 270, 220, 4000., 9.80665, 8, 45, 6.17, 0.9, 1, 1.15, 0.75)
# aircraft.performFlight(1.225)
# data = aircraft.data
# timestamps=aircraft.timestamps
# powert=data['powert']
# energy=data['energy']
# time=data['time']
# powert=np.around(powert, decimals=1)
# time=np.around(time, decimals=1)
# from IterationTools.Mission import TiltRotor
# ##
# Req=np.amax(powert)/2
# Data, weight, name=libelectricmotor()
# name, number, mass, rating = Componentspicker.choosecomponents(Req, Data, weight, name)
# print("name=",name)
# print("number=",number)
# print("mass=",mass)
# print("rating=",rating)
##

