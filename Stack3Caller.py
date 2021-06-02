#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  8 10:38:18 2021

@author: michaelpokornik
"""
import matplotlib.pyplot as plt
import numpy as np
import ElectronDataClass as EDC
import PhotonDataClass as PDC
import FieldClass as FC
import matplotlib.animation as ani
import make_3stack_frame as make3
import DensityClass as DC
#%%electrons and photons
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
target = 'uniform' #uniform or channel
cut = '30deg' # or normal
run = 'run5'
EnergyUnits = 'MeV'
duration = 31

noLabel=True #if true, then don't include the x,y labels 
#Restricte = {'XLim':(-5,50,1e6),'YLim':(-15,15,1e6),'E0Lim':(20,1000,J2MeV)}
Restricte = {'E0EtaLim':(1,J2MeV)}
#Restrictp = {'E0Lim':(1,300,J2MeV)}
Restrictp = {}
paramse = {'phi':(-180,180,181), 'E0':(0,900,301,J2MeV)}
paramsp = {'phi':(-180,180,181), 'E0':(1,400,201,J2MeV)}

ED = EDC.ElectronDataClass({'Path':'/Volumes/Elements/targets_eta/' + f'{target}' + '/' + f'{cut}' + '/' + f'{run}' + '/'})
Ph = PDC.PhotonDataClass({'Path':'/Volumes/Elements/targets_eta/' + f'{target}' + '/' + f'{cut}' + '/' + f'{run}' + '/'})

param_stre = '_e_'
param_strp = '_p_'

for key in Restricte:
    if key != 'E0EtaLim':
        param_stre = param_stre + f'{key}_{Restricte[key][0]}_{Restricte[key][1]}_'
    else:
        param_stre = param_stre + f'{key}_{Restricte[key][0]}_{EnergyUnits}_'      
for key in Restrictp:
    param_strp = param_strp + f'{key}_{Restrictp[key][0]}_{Restrictp[key][1]}_'
    
param_str = param_stre + param_strp

pltName =  f'elec_ph_ang_2D{target}_{cut}_{param_str}_{run}_moreticks.gif'
#fig,ax = plt.subplots(2,1,figsize=(15,15))
fig,ax = plt.subplots(3,1,figsize=(8,12))
#%%Fields
n = 80
SMALL_SIZE = 24
MEDIUM_SIZE = 34
BIGGER_SIZE = 40
normFac = 3.103e14 #E0
xticks=[0,20,40,60]
yticks=[-10,0,10]
Fields = FC.FieldClass({'Path':ED.Path})
#%%Density
species = 'Electron'
Dens = DC.DensityClass({'Path':ED.Path})
#%%run      
#with density
fig=make3.make_3stack_frame(ED,Ph,Fields,noLabel,fig,f'ang_{target}_{cut}',22,'e_loc',Restricte,
        paramse,'PhotonDetail',Restrictp,paramsp,False,[-60,-30,0,30,60],[0,300,500,800],[10,14],[-60,-30,0,30,60],[0,200,400],
        [8,19],[0,60],[-15,15],[0,1],3.013e14,'e_fields','|Ey/E0|',n,
        'Electric_Field_Ey',[-10,0,10,20,30,40,50],[-10,-5,0,5,10],24,30,30,Dens,species,ncrit,[0,15,30],False)



