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

def CrossSections(CS,EnC,ind):
    if CS['sliceD'] == 'x':
        dat = EnC.Ek[ind,:,:]
    elif CS['sliceD'] == 'y':
        dat = EnC.Ek[:,ind,:]
    else:
        dat = EnC.EK[:,:,ind]
    return dat

MeV = (1.6022e-19)*(1.0e6) 
c_speed = 299792458 # m/s
mc = 9.11e-31 * c_speed #kg m /s
J2MeV = 1 * (6.242e18) / (1.0e6)

target = 'targets_3d' #uniform or channel
EfileName = 'dens'
species = 'Derived_Average_Particle_Energy_averaged_Electron_t'
run = 'run1'
duration = 120
n=32

EnergyUnits = 'MeV'
EnC = DC.DensityClass({'Path':'/Volumes/Elements/' + f'{target}' + '/' + f'{run}' + '/'})
EnC.LoadEkbarAtTime(EfileName,n,species,(1/J2MeV))
cbarV = [0,15,30]
nlevels = 20
clevels = np.linspace(cbarV[0],cbarV[-1],num=nlevels)
#slices of the 3d
CrossSection = {'pltx':EnC.z,'plty':EnC.x,'sliceD':'y'}
xvalues = np.linspace(-3,3,num=7)
xvalues = [-2.5,-2,0,2,2.5]
zvalues = np.linspace(-3,6,num=10)
zvalues = [-2.5,0,2.5,5,7.5]
yval = [2e-6,1e-6,0]
subt = ['y=2um','y=1um','y=0um']
ln=len(subt)
xg,yg = np.meshgrid(EnC.z*1e6,EnC.x*1e6)

lblsz=36
ticksz=30
cbarft=30
cbart=40
ind = {}


for i, val in enumerate(yval):
   ind[i] = find_nearest(EnC.y,val)
   
#fig,ax = plt.subplots(1,n,figsize=(20,10))
fig,ax = plt.subplots(ln,1,figsize=(10,20))


for i in range(ln):
    dat = CrossSections(CrossSection, EnC, ind[i])
    pltM = ax[i].contourf(xg,yg,dat,cmap='jet',vmin=cbarV[0],vmax=cbarV[-1],extend='both',alpha=.75,levels=clevels)
    ax[i].set_title(subt[i],color='white',fontsize=40,loc='left',x=0, y=0.8)
    ax[i].set_xticks(zvalues[1:-1])
    ax[i].set_yticks(xvalues[1:-1])
    ax[i].set_xlim(zvalues[0],zvalues[-1])
    ax[i].set_ylim(xvalues[0],xvalues[-1])
    ax[i].tick_params(axis='x',labelsize=ticksz)
    ax[i].tick_params(axis='y',labelsize=ticksz)
    ax[i].grid(color='grey', linestyle='-', linewidth=.5)

#ax[0].set_xticklabels([])
ax[-1].set_xlabel('Z[um]',fontsize=lblsz)
ax[1].set_ylabel('X[um]',fontsize=lblsz)
ax[0].set_xticklabels([])
ax[1].set_xticklabels([])
#ax[3].set_yticklabels([])

fig.subplots_adjust(bottom=.2)
fig.subplots_adjust(hspace=0.02, wspace=0)
cbar_ax = fig.add_axes([0.3,.95,.4,.05])
cbar=fig.colorbar(pltM, cax=cbar_ax,orientation='horizontal')
cbar.set_ticks(cbarV)
cbar.ax.tick_params(labelsize=cbarft)
cbar.ax.set_title(r'$\overline{E_{k}}$ [MeV]',fontsize=cbart)
# plt.colorbar(m, boundaries=np.linspace(0, 2, 6))