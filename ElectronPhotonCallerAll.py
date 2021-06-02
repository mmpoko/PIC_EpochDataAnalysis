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
import matplotlib.animation as ani


MeV = (1.6022e-19)*(1.0e6) 
c_speed = 299792458 # m/s
mc = 9.11e-31 * c_speed #kg m /s
J2MeV = 1 * (6.242e18) / (1.0e6)

target = 'channel' #uniform or channel
cut = '30deg' # or normal
run = 'run1'
EnergyUnits = 'MeV'
duration = 31
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
fig,ax = plt.subplots(3,1,figsize=(15,15))
#%%

def make_frame(ED,Ph,*args):
    
    fig = args[0]
    axs = args[1]
    axs = fig.get_axes()
    #axs = fig.gca()
    plotName = args[2]
    n = args[3]
    ename = args[4]
    Restricte = args[5]
    paramse = args[6]
    pname = args[7]
    Restrictp = args[8]
    paramsp = args[9]
    
    t = ED.LoadElectronsAtTime(f'{ename}',n,Restricte)
    t2 = Ph.LoadPhotonsAtTime(f'{pname}',2*n,Restrictp)
    Ph.KeepOnlyEmittingParticles(ED)
    

    fig.suptitle(fr'{plotName} ' +  fr'{t * 1e15:.2f} fs',fontsize=22)
    
    #%% if pseudocolor
    HistV, Histx, Histy  = ED.Pseudo2D(paramse)
    X,Y = np.meshgrid(Histx,Histy)
    logHist = np.log(HistV)    
    im = axs[0].pcolormesh(X,Y,np.transpose(logHist),vmin=np.nanmin(logHist),vmax=np.nanmax(logHist), cmap=plt.get_cmap('jet'))
    im.set_clim([9,13])
    xtickmarks = np.array([-180,-150,-120,-90,-60,-30,0,30,60,90,120,150,180])
    axs[0].grid(color='grey', linestyle='-', linewidth=.5)
    ce = fig.colorbar(im,ax=axs[0])  
    #ce.set_clim([9,13])
    axs[0].set_xlabel(r'$\Phi [\circ]$',fontsize=20)
    axs[0].set_ylabel(r'$E [MeV]$',fontsize=20)
    ce.set_label(r'electron $\frac{d^2 logN}{d\phi dE}$',fontsize=20)
    axs[0].set_xticks(xtickmarks)
    axs[0].tick_params(labelsize=16)
    axs[0].tick_params(axis='x',rotation=-45)
    axs[0].axvline(color="grey",linewidth=2)
    axs[0].set_xlim((-90,90))
    axs[0].set_ylim((0,900))
    if hasattr(Ph,'peid'):
        if Ph.peid.size !=0:
            HistVp, Histxp, Histyp = Ph.Pseudo2D(paramsp)
            Xp,Yp = np.meshgrid(Histxp,Histyp)
            logHistp = np.log(HistVp)
            imp = axs[1].pcolormesh(Xp,Yp,np.transpose(logHistp),vmin=np.nanmin(logHistp),vmax=np.nanmax(logHistp), cmap=plt.get_cmap('jet'))
            cp = fig.colorbar(imp,ax=axs[1]) 
            imp.set_clim([9,16])
            #cp.set_clim([10,16])
            cp.set_label(r'photon $\frac{d^2 logN}{d\phi dE}$',fontsize=20)
    
    axs[1].grid(color='grey', linestyle='-', linewidth=.5) 
    axs[1].set_xlabel(r'$\Phi [\circ]$',fontsize=20)
    axs[1].set_ylabel(r'$E [MeV]$',fontsize=20)
    axs[1].set_xticks(xtickmarks)
    axs[1].tick_params(labelsize=16)
    axs[1].tick_params(axis='x',rotation=-45)
    axs[1].axvline(color="grey",linewidth=2)
    axs[1].set_xlim((-90,90))
    axs[1].set_ylim((0,300))


    #%% if lineout
    # #electron
    # LineOute,phiMidPoints = ED.LineOut('phi','E0',J2MeV)
    # axs[0].plot(phiMidPoints,LineOute)
    # axs[0].set_xlabel(r'$\Phi [\circ]$',fontsize=16)
    # #ax0.set_ylabel(r'$E [MeV]$',fontsize=16)
    # axs[0].set_ylabel(r'$dEd\Phi$',fontsize=16)
    # axs[0].set_xticks(np.linspace(-180,180,11))
    # axs[0].set_ylim([0, .01])
    # axs[0].grid(color='k', linestyle='-', linewidth=.5)
    # #photon
    # if hasattr(Ph,'peid'):
    #     if Ph.peid.size !=0:
    #         LineOutp,phiMidPoints = ED.LineOut('phi','E0',J2MeV)
    #         axs[1].plot(phiMidPoints,LineOutp)
            
    # axs[1].grid(color='k', linestyle='-', linewidth=.5) 
    # axs[1].set_xlabel(r'$\Phi [\circ]$',fontsize=16)
    # axs[1].set_ylabel(r'$dEd\Phi$',fontsize=16)
    # axs[1].set_xticks(np.linspace(-180,180,11))
    # axs[1].set_ylim([0, .01])
    
#%%


moviewriter = ani.PillowWriter(fps=10)
with moviewriter.saving(fig,pltName,duration):
    for j in range(duration):
        make_frame(ED,Ph,fig,ax,f'ang_{target}_{cut}',j,'e_loc',Restricte,paramse,'PhotonDetail',Restrictp,paramsp)
        #plt.figure(num=1, figsize=(10,20),dpi=100) 
        fig.tight_layout(pad=3.0)
        moviewriter.grab_frame()
        fig.clf()
        fig.add_subplot(121)
        fig.add_subplot(122)

    moviewriter.finish()

