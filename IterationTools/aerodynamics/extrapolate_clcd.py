import numpy as np
import pandas as pd
import os

names = ['0009', '0011', '0012', '0018', '23015', '64210']
columns = ['alpha', 'cl', 'cd', 'cd_p', 'cm', 'top_x', 'bot_x']

dfs = []
names = ['0012']
for name in names:
    lower = []
    lowerlower = []
    upper = []
    upperupper = []
    # print(name)
    info = pd.read_excel('0012clcd.xlsx').values#, sheet_name=name).values
    # print(info)

    clmax = np.max(info[:, 1])
    idx_stall = np.where(info[:, 1] == clmax)
    alfa_stall = info[idx_stall, 0]
    cd_stall = info[idx_stall, 2]

    cdmax = 1.11+0.018*10
    A1 = cdmax/2
    B1 = cdmax
    # print(np.max(info[:, 1]))
    A2 = (clmax -cdmax*np.sin(np.radians(alfa_stall))*np.cos(np.radians(alfa_stall)))*np.sin(np.radians(alfa_stall))/(np.cos(np.radians(alfa_stall))**2)
    B2 = (cd_stall-cdmax*np.sin(np.radians(alfa_stall))**2)/np.cos(np.radians(alfa_stall))

    for alfa in np.arange(info[-1, 0]+0.25, 90, 0.25):
        cl = A1 * np.sin(np.radians(2*alfa)) + A2*np.cos(np.radians(alfa))**2/np.sin(np.radians(alfa))
        cd = B1*np.sin(np.radians(alfa))**2 + B2*np.cos(np.radians(alfa))
        new_entry = [alfa, cl, cd, 0, 0, 0, 0]
        upper.append(new_entry)

    for alfa in np.arange(-90, info[0, 0]-0.2, 0.25):
        cl = A1 * np.sin(np.radians(2*alfa)) + A2*np.cos(np.radians(alfa))**2/np.sin(np.radians(alfa))
        cd = B1*np.sin(np.radians(alfa))**2 + B2*np.cos(np.radians(alfa))
        new_entry = [alfa, cl, cd, 0, 0, 0, 0]
        lower.append(new_entry)


    info = np.append(info, upper, axis=0)
    info = np.append(lower, info, axis=0)
    idx_mid = np.where(info[:,0] == 0)[0]

    print(info)
    print(name)
    for i, alfa in enumerate(np.arange(-180, -90.1, 0.25)):
        # print('start')
        # print(info[idx_mid, 0])
        # print(idx_mid)
        # print(i)
        # print(idx_mid+i)
        # print(len(info))
        cl = info[idx_mid+i,1]
        cd = info[idx_mid+i,2]
        new_entry = [alfa, cl, cd, 0, 0, 0, 0]
        lowerlower.append(new_entry)

    for i, alfa in enumerate(np.arange(90, 180.1, 0.25)):
        cl = info[i,1]
        # print(cl)
        cd = info[i,2]
        new_entry = [alfa, cl, cd, 0, 0, 0, 0]
        upperupper.append(new_entry)

    info = np.append(lowerlower, info, axis=0)
    info = np.append(info, upperupper, axis=0)

    info = pd.DataFrame(info, columns=columns)
    dfs.append(info)


pre = os.path.dirname(os.path.realpath(__file__))
fname = 'airfoils_try.xlsx'
path = os.path.join(pre, fname)

writer = pd.ExcelWriter(path, engine='xlsxwriter')
for i, df in enumerate(dfs):
    df.to_excel(writer, sheet_name=names[i], index=False)

writer.save()
