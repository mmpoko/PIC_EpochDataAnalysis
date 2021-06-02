#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 12 21:51:37 2021

@author: michaelpokornik
"""
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
import FieldClass as FC
import matplotlib.animation as ani
import DensityClass as DC
from matplotlib.colors import LinearSegmentedColormap
from mpl_toolkits.axes_grid1 import make_axes_locatable

#%%
fig,axs = plt.subplots(1,1,figsize=(10,6))
#%%
n=80
SMALL_SIZE = 30
MEDIUM_SIZE = 40
BIGGER_SIZE = 40
normFac = 3.103e14 #E0
xticks=[-10,0,10,20,30,40,50]
yticks=[-10,-5,0,5,10]
cdticks = [0,15,30]
cfticks = [0,.5,1]
#ncrit cal
omega_l = 2*np.pi*((3e8)/(.8e-6)) #rad/s
e = 1.602e-19 #c
m_e = 9.109e-31 #kg
# for gaussian ncrit = m_e * omega_l * omega_l /4/np.pi/e/e
eps0 = 8.85e-12
ncrit = m_e * eps0 * omega_l * omega_l /e/e
subSampInt = 1
cbarOn = False
#efields
#Fields = FC.FieldClass({'Path':'/Volumes/4TBWD/PhotonProject/short_pulse/uniform/30deg/r1/'})

Fields = FC.FieldClass({'Path':'/Volumes/Elements/targets_eta/uniform/30deg/run4/'})
#Density
species = 'Electron'
Dens = DC.DensityClass({'Path':Fields.Path})
#%% plot fields
t = Fields.LoadFieldAtTime('e_fields',n,'Electric_Field_Ey')
#t = Fields.LoadFieldAtTime('b_fields',n,'Magnetic_Field_Bz_')
t2=Dens.LoadDensityAtTime('dens',n,species,ncrit)
dcolor = 'jet'
dalphacolor = 'jet_alpha'

#absey = np.abs(Fields.bz[::subSampInt,::subSampInt] /1e6)
absey = np.abs(Fields.ey[::subSampInt,::subSampInt] / normFac)
#absey[absey<.5]=np.nan
den = Dens.dens[::subSampInt,::subSampInt]
xp,yp = np.meshgrid(Dens.x[::subSampInt] * 1e6,Dens.y[::subSampInt] * 1e6)
# get blue red alpha colormap#################
ncolors = 256
color_array = plt.get_cmap('Greys')(range(ncolors))
    # change alpha values
color_array[:,-1] = np.linspace(0.0,1.0,ncolors)
#color_array[64:,-1] = 1
    # create a colormap object
map_object = LinearSegmentedColormap.from_list(name='greys_alpha',colors=color_array)
plt.register_cmap(cmap=map_object)
m = plt.cm.ScalarMappable(cmap=plt.get_cmap('greys_alpha'))
m.set_array(np.transpose(absey))
m.set_clim([cfticks[0],cfticks[-1]])
# get density colormap#########################
ncolors = 256
color_array = plt.get_cmap(dcolor)(range(ncolors))
# change alpha values
color_array[:,-1] = np.linspace(0.0,1.0,ncolors)
# create a colormap object
map_object = LinearSegmentedColormap.from_list(name=dalphacolor,colors=color_array)
# register this new colormap with matplotlib
plt.register_cmap(cmap=map_object)
m1 = plt.cm.ScalarMappable(cmap=plt.get_cmap(dalphacolor))
m1.set_array(np.transpose(den))
m1.set_clim([cdticks[0],cdticks[-1]])

#plot pcolor
imd = axs.pcolormesh(xp,yp,np.transpose(den),vmin=cdticks[0],vmax=cdticks[-1],cmap=plt.get_cmap(dalphacolor))
xd,yd = np.meshgrid(Fields.x[::subSampInt]*1e6,Fields.y[::subSampInt]*1e6)
imf = axs.contour(xd,yd,np.transpose(absey),vmin=cfticks[0],vmax=cfticks[-1],cmap=plt.get_cmap("greys_alpha"))

#set up colorbars
if cbarOn == True:
    cd = fig.colorbar(m1,ax=axs,orientation='horizontal',shrink=.2,pad=.3)
    cd.set_label(r'$n/n_{cr}$',fontsize=MEDIUM_SIZE)
    cd.set_ticks(cdticks)
    cd.ax.tick_params(labelsize=MEDIUM_SIZE)    
    cf = fig.colorbar(m,ax=axs,shrink=.4)
    cf.set_label(r'$Bz [MT]$',fontsize=MEDIUM_SIZE)
    #cf.set_label(r'$|Ey/E0|$',fontsize=MEDIUM_SIZE)
    cf.set_ticks(cfticks)
    cf.ax.tick_params(labelsize=MEDIUM_SIZE)

axs.set_xlabel('x [$\mu m$]')
axs.set_xticks(xticks)
axs.set_yticks(yticks)
axs.set_xlim([xticks[0],xticks[-1]])
axs.set_ylim([yticks[0],yticks[-1]])
axs.grid(color='grey', linestyle='-', linewidth=.5)
axs.set_ylabel('y [$\mu m$]')

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=MEDIUM_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
#plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

fig.tight_layout()


