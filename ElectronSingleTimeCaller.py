#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  9 13:35:20 2021

@author: michaelpokornik
"""
"""
Created on Mon Mar  8 10:38:18 2021

@author: michaelpokornik
"""
import matplotlib.pyplot as plt
import numpy as np
import ElectronDataClass as EDC



target = 'uniform' #uniform or channel
cut = '30deg' # or normal
#ED = EDC.ElectronDataClass({'Path':'/Volumes/4TBWD/PhotonProject/track_carbon/uniform/30deg/'})
ED = EDC.ElectronDataClass({'Path':'/Volumes/4TBWD/PhotonProject/track_carbon/' + f'{target}' + '/' + f'{cut}' + '/'})

#%%
t = ED.LoadElectronsAtTime('id',40,{'XLim':(-5,5,1e6),'YLim':(-4,4,1e6)})
rad_angle = np.arctan2(ED.py[ED.eleInd > 0],ED.px[ED.eleInd > 0])
deg_angle = np.rad2deg(rad_angle)
phi_bins = np.linspace(-180, 180,181)
digitized = np.digitize(deg_angle,phi_bins)
normCounts = ED.w[ED.eleInd > 0] / np.mean(np.diff(phi_bins))
phiMidPoints = (phi_bins[1:] + phi_bins[:-1]) / 2
LineOut = np.bincount(digitized-1, weights=normCounts)
LineOutMissing = np.zeros((phiMidPoints.shape[0] - LineOut.shape[0]))
LineOut = np.append(LineOut,LineOutMissing)
fig = plt.figure()
ax = fig.gca()
ax.plot(phiMidPoints,LineOut)
#plt.rcParams.update({'font.size': 14})
ax.set_title(r'uniform ' + f'{cut}'  + f' {t * 1e15:.2f} fs',fontsize=22)
ax.set_xlabel(r'$\Phi [\circ]$',fontsize=16)
ax.set_ylabel(r'$\frac{dN}{d\Phi}$',fontsize=16)
ax.set_xticks(np.linspace(-180,180,4))
ax.set_ylim([0, 1e13])
ax.grid()
