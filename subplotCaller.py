#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 22 15:19:10 2021

@author: michaelpokornik
"""
import matplotlib.pyplot as plt
import numpy as np
import ElectronDataClass as EDC
import PhotonDataClass as PDC
import FieldClass as FC
import matplotlib.animation as ani
import makesubplotFrame as makesub
import makesubplotFrameBWR as makesubbwr
import DensityClass as DC
import matplotlib.gridspec as gridspec
#%%create figure
#fig,ax = plt.subplots(2,3,figsize=(16,8),gridspec_kw={'width_ratios': [2,3,3]})
fig,ax = plt.subplots(1,1,figsize=(8,4))
#%%run
MeV = (1.6022e-19)*(1.0e6) 
c_speed = 299792458 # m/s
mc = 9.11e-31 * c_speed #kg m /s
J2MeV = 1 * (6.242e18) / (1.0e6)
omega_l = 2*np.pi*((3e8)/(.8e-6)) #rad/s
e = 1.602e-19 #c
m_e = 9.109e-31 #kg
# for gaussian ncrit = m_e * omega_l * omega_l /4/np.pi/e/e
eps0 = 8.85e-12
ncrit = m_e * eps0 * omega_l * omega_l /e/e
ncrit_c = ncrit/6.0
m_c = m_e * 1836 * 12
wpi2 = ncrit_c * (6*e)*(6*e)/m_c/eps0
wpi=np.sqrt(wpi2)
ts = 1/wpi
print(ts)
target = 'uniform_im' #uniform or channel
cut = '30deg' # or normal
run = 'run1'
EnergyUnits = 'MeV'
species = 'Electron'
Restricte = {'E0EtaLim':(1,J2MeV)}
Restrictp = {}
paramse = {'phi':(-180,180,181), 'E0':(0,900,301,J2MeV)}
paramsp = {'phi':(-180,180,181), 'E0':(1,400,201,J2MeV)}

x_particle = [-60,-30,0,30,60]
y_elec = [0,300,600,750]
y_ph =[1,150,300,400]
n = 70
n2=n
ebarlims=[8,14]
phbarlims=[10,20]

normFac = 3.103e14 #E0
xticks=[0,20,40]
yticks=[-15,-10,-5,0,5,10,15]
cbareabslims = [0,1]
cbarelims = [-1,1]
#cbardlims = [0,130]
cbardlims = [0,15]
#cbardlims for uniform is [0 30]


#create the class objects
ED = EDC.ElectronDataClass({'Path':'/Volumes/Elements/targets_eta/' + f'{target}' + '/' + f'{cut}' + '/' + f'{run}' + '/'})
#ED = EDC.ElectronDataClass({'Path':'/Volumes/4TBWD/PhotonProject/short_pulse/' + f'{target}' + '/' + f'{cut}' + '/' + f'{run}' + '/'})
Ph = PDC.PhotonDataClass({'Path':ED.Path})
Fields = FC.FieldClass({'Path':ED.Path})
Dens = DC.DensityClass({'Path':ED.Path})


#'Ylbl':(r'$E_{electron}$ [MeV]')
#r'$E_{photon}$ [MeV]'
#r'$Y[\mu m]$'
plotN = 0
elecInputs={'fileName':('e_loc'),'xlims':(x_particle),'ylims':(y_elec),'cbarOn':(True),'XAxLbl':(False),'YAxLbl':(False),'Xlbl':(r'$\theta[^{\circ}]$'),
            'Ylbl':(''),'params':(paramse),'LastYTick':(False),
            'cbarlims':(ebarlims),'cbarLbl':(r'$\frac{d^2 logN_{e}}{d\theta dE_{e}}$'),'Restrict':(Restricte),'time':(n)}

phInputs={'fileName':('PhotonDetail'),'ClassObj':(Ph),'xlims':(x_particle),'ylims':(y_ph),'cbarOn':(True),'XAxLbl':(False),'YAxLbl':(True),
          'Xlbl':(r'$\theta[^{\circ}]$'),'Ylbl':(''),'params':(paramsp),'LastYTick':(False),
            'cbarlims':(phbarlims),'cbarLbl':(r'$\frac{d^2 logN_{ph}}{d\theta dE_{ph}}$'),'Restrict':({}),'time':(n2)}

fieldInputs={'fileName':('e_fields'),'field':('Electric_Field_Ey'),'XAxLbl':(True),'YAxLbl':(True),'Xlbl':(r''),'Ylbl':(r''),
                 'time':(n),'norm':(normFac),'xlims':(xticks),'ylims':(yticks),'cbarOn':(False),'cbarLbl':(r'$E_{y}/E_{0}$'),'LastYTick':(False),'FirstYTick':(False),
                 'FirstXTick':(True),'LastXTick':(True),'cbarlims':(cbarelims)}

DensInputs={'fileName':('dens'),'species':('Electron'),'time':(n),'norm':(ncrit),'xlims':(xticks),'ylims':(yticks),'logscale':(False),'Xlbl':(r''),
            'LastYTick':(True),'FirstYTick':(True),'Ylbl':(''),'cbarOn':(False),'XAxLbl':(True),'YAxLbl':(False), 'cbarLbl':(r'$n_{e}/n_{cr}$'),'cbarlims':(cbardlims)}

if DensInputs['logscale']==True:
    DensInputs['cbarlims']=[.1,250]
    DensInputs['cbarLbl'] = r'$\log{|n_{C}/n_{cr}|}$'
    DensInputs['lgcbarlims']=[.1,1,10,100]
    
    
TICK_SIZE = 24
SMALL_SIZE = 28
MEDIUM_SIZE = 30
BIGGER_SIZE = 40
vshrink=.5
hshrink=.5
vasp=8
hasp=12
plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=MEDIUM_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
#plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title
fontArgs = [TICK_SIZE,SMALL_SIZE,MEDIUM_SIZE,BIGGER_SIZE,vshrink,hshrink,vasp,hasp]
#fig=makesub.makesubplotFrame(fig,plotN,'electron',ED,Ph,Fields,Dens,elecInputs,phInputs,fontArgs)
#fig=makesub.makesubplotFrame(fig,plotN,'field',ED,Ph,Fields,Dens,fieldInputs,DensInputs,fontArgs)

fig=makesubbwr.makesubplotFrameBWR(fig,plotN,'field',ED,Ph,Fields,Dens,fieldInputs,DensInputs,fontArgs)
ax=fig.get_axes()
ax[plotN].tick_params(axis='both',labelsize=MEDIUM_SIZE)
ax[plotN].spines["bottom"].set_linewidth(2)
ax[plotN].spines["left"].set_linewidth(2)
ax[plotN].spines["right"].set_linewidth(2)
ax[plotN].spines["top"].set_linewidth(2)
fig.subplots_adjust(wspace=0.0, hspace=0.04)
#%% set up gridspec template
# gs1 = gridspec.GridSpec(1,5) #template allows for 3 to 2 ratio
# gs1.update(left=0.1,right=.65,wspace=.05)
# ax1=plt.subplot(gs1[0,:-2])
# ax2=plt.subplot(gs1[0,3:])
# #next one
# gs2= gridspec.GridSpec(1,1)
# gs2.update(left=.75,right=.9)
# ax3=plt.subplot(gs2[0,0])

#move subplots
# ax[0].set_position(gs1[0,:-2].get_position(fig))
# ax[1].set_position(gs1[0,3:].get_position(fig))
# ax[2].set_position(gs2[0,0].get_position(fig))






