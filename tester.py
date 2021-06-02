#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 16 14:17:52 2021

@author: michaelpokornik
"""
from joblib import Parallel, delayed
import multiprocessing

inputs = range(10) 
def processInput(i):
    return i * i

num_cores = multiprocessing.cpu_count()

results = Parallel(n_jobs=num_cores)(delayed(LocatePh)(ex,ey,phx,phy,i) for i in inputs)
print(results)
def LocatePh(ex,ey,phx,phy,i):
    locx = np.array(np.isin(ex,phx[i]))
    # locy = np.array(np.isin(ey,phy[i]))
    # loct = np.array(locx * locy)
    # maxv = np.max(loct)
    maxv = np.max(locx)
    if maxv == 1:
        return (i)
    else:
        print("not " ,i)
    
    
ex = ED.x[ED.eleInd > 0]
ex = ex[:,np.newaxis]
ey = ED.y[ED.eleInd > 0]
ey = ey[:,np.newaxis]
epos = np.hstack((ex,ey))
epos = epos[epos[:,0].argsort()]
phx = Ph.x
phx = phx[:,np.newaxis]
phy = Ph.y
phy = phy[:,np.newaxis]
ppos = np.hstack((phx,phy))
ppos = ppos[ppos[:,0].argsort()]
m = np.size(ppos,0)
n = np.size(epos,0)
xi = 0
exi = 0
ind = np.zeros((m,))
   
while exi < n:
    if xi == m:
        break
    if epos[exi][0] < ppos[xi][0]:
        exi +=1
    elif epos[exi][0] == ppos[xi][0]:
        tempx = xi
        while epos[exi][0] == ppos[tempx][0]:
            if epos[exi][1] != ppos[tempx][1]:
                tempx +=1
            else:
                ED.eleInd[exi] = 1
                exi+=1
                break
    else:
        xi+=1
                        
    
ex = ED.x[ED.eleInd > 0]
ex = ex[:,np.newaxis]
ey = ED.y[ED.eleInd > 0]
ey = ey[:,np.newaxis]
epos = np.hstack((ex,ey))
phx = Ph.x
phx = phx[:,np.newaxis]
phy = Ph.y
phy = phy[:,np.newaxis]
ppos = np.hstack((phx,phy))
m = np.size(ppos,0)
n = np.size(epos,0)
xi = 0
exi = 0
ind = np.zeros((m,))
t0 = time.time()
while xi < m:
    if exi == n:
        break
    if ppos[xi][0] < epos[exi][0]:
        xi +=1
    elif ppos[xi][0] == epos[exi][0]:
        print('look')
        if ppos[xi][1] == epos[exi][1]:
            ind[xi] = 1
        else:
            exi+=1
    else:
        exi+=1
    #print(xi) 
t01 = time.time()
print(t01-t0)




# def processInput(ED,Ph,i):
#     locx = np.array(np.isin(ED.x[ED.eleInd > 0],Ph.x[i]))
#     locy = np.array(np.isin(ED.y[ED.eleInd > 0],Ph.y[i]))
#     loct = np.array(locx * locy)
#     maxv = np.max(loct)
#     if maxv == 1:
#         return (Ph.x[idx],Ph.y[idx])
# inputs = range(Ph.x.size)
# (phx,phy) = Parallel(n_jobs=3)(delayed(processInput)(ED,Ph,i) for i in inputs) 
    
# print(phx,phy)