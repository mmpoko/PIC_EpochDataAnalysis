#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import sdf_helper as sdf_helper
import DensityClass as DC
import numpy as np
from scipy import stats
from matplotlib.colors import LinearSegmentedColormap
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.colors as col


def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx


MeV = (1.6022e-19)*(1.0e6) 
c_speed = 299792458 # m/s
mc = 9.11e-31 * c_speed #kg m /s
J2MeV = 1 * (6.242e18) / (1.0e6)

target = 'targets_3d' #uniform or channel
EfileName = 'dens'
species = 'Derived_Average_Particle_Energy_averaged_Electron_t'
run = 'run1'
duration = 120
n=5

EnergyUnits = 'MeV'
EnC = DC.DensityClass({'Path':'/mnt/e/' + f'{target}' + '/' + f'{run}' + '/'})
EnC.LoadEkbarAtTime(EfileName,n,species,(1/J2MeV))
cbarV = [0,1.75,3.25]
nlevels = 20
#slices of the 3d
ind = {}
xvalues = np.linspace(-3,3,num=7)
xvalues = [-3,-2,0,2,3]
yvalues = np.linspace(-3,3,num=7)
yvalues = [-3,-2,0,2,3]
zval = [1e-6,3e-6,5e-6,5.2e-6]
subt = ['z=1um','z=3um','z=5um','z=5.2um']
lblsz=36
ticksz=30
cbarft=30
cbart=40
for i, val in enumerate(zval):
   ind[i] = find_nearest(EnC.z,val)



fig,ax = plt.subplots(1,4,figsize=(20,10))
clevels = np.linspace(cbarV[0],cbarV[-1],num=nlevels)
xg,yg = np.meshgrid(EnC.x*1e6,EnC.y*1e6)
for i in range(4):
    pltM = ax[i].contourf(xg,yg,np.transpose(EnC.Ek[:,:,ind[i]]),cmap='jet',vmin=cbarV[0],vmax=cbarV[-1],extend='both',alpha=.75,levels=clevels)
    ax[i].set_title(subt[i],fontsize=22)
    ax[1].set_xlabel('X[um]',fontsize=lblsz)
    ax[i].set_xticks(xvalues[1:-1])
    ax[0].set_ylabel('Y[um]',fontsize=lblsz)
    ax[i].set_yticks(yvalues[1:-1])
    ax[i].tick_params(axis='x',labelsize=ticksz)
    ax[i].tick_params(axis='y',labelsize=ticksz)
    ax[i].grid(color='grey', linestyle='-', linewidth=.5)

#ax[0].set_xticklabels([])
ax[1].set_yticklabels([])
ax[2].set_yticklabels([])
ax[3].set_yticklabels([])

fig.subplots_adjust(top=0.7)
fig.subplots_adjust(hspace=0, wspace=.05)
cbar_ax = fig.add_axes([0.41,.8,.2,.1])
cbar=fig.colorbar(pltM, cax=cbar_ax,orientation='horizontal')
cbar.set_ticks(cbarV)
cbar.ax.tick_params(labelsize=cbarft)
cbar.ax.set_title(r'$\overline{E_{k}}$ [MeV]',fontsize=cbart)
# plt.colorbar(m, boundaries=np.linspace(0, 2, 6))