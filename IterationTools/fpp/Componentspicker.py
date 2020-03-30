##
#from IterationTools.fpp import Componentslibrary
import numpy as np

def choosecomponents(Req, Data, weight, name):
    nFC=np.ceil(Req/Data)
    mass=weight*(nFC+1)
    minm=np.amin(mass)
    i=np.where(mass==minm)[0]
    return name[i], nFC[i], mass[i], nFC[i]* Data[i]
##
# Req=361
# Data=np.array([49, 63, 81, 98, 125])
# weight=np.array([21, 25, 29, 34, 43])
# name=np.array(["A167", "B215", "C275", "D335", "E455"])
##
#name, number, mass, rating = choosecomponents(Req, Data, weight, name)
##
# print("name=",name)
# print("number=",number)
# print("mass=",mass)
# print("rating=",rating)
##
# Data, weight, name=Componentslibrary.libfuelcell()
# Req=np.array([200, 300, 500, 50])
def choosecomponentsarray(Req, Data, weight, name):
    #name2=np.chararray(len(Req), unicode=True)
    name2=[]
    number = np.empty(len(Req))
    mass = np.empty(len(Req))
    rating = np.empty(len(Req))

    for i in range(len(Req)):
        name3, number3, mass3, rating3 = choosecomponents(Req[i], Data, weight, name)
        name2.append(name3[0])
        number[i], mass[i], rating[i]= number3[0], mass3[0], rating3[0]
    return name2, number, mass, rating

# name2, number, mass, rating= choosecomponentsarray(Req, Data, weight, name)
##
# def choosecomponents(Req1, Req2, Data1, Data2, name):
#
#     mass=weight*(nFC+1)
#     minm=np.amin(mass)
#     i=np.where(mass==minm)[0]
#     return name[i], nFC[i], mass[i], nFC[i]* Data[i]
# def choosecomponentsarraybattery(Req1, Req2, Data1, Data2, name):
#     #name2=np.chararray(len(Req), unicode=True)
#     name2=[]
#     number = np.empty(len(Req))
#     mass = np.empty(len(Req))
#     rating = np.empty(len(Req))
#
#     for i in range(len(Req)):
#         name3, number3, mass3, rating3 = choosecomponents(Req[i], Data, weight, name)
#         name2.append(name3[0])
#         number[i], mass[i], rating[i]= number3[0], mass3[0], rating3[0]
#     return name2, number, mass, rating