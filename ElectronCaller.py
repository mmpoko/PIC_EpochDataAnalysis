#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  8 10:38:18 2021

@author: michaelpokornik
"""
import matplotlib.pyplot as plt
import numpy as np
import ElectronDataClass as EDC
import matplotlib.animation as ani


MeV = (1.6022e-19)*(1.0e6) 
c_speed = 299792458 # m/s
mc = 9.11e-31 * c_speed #kg m /s
J2MeV = 1 * (6.242e18) / (1.0e6)

target = 'uniform' #uniform or channel
pltType = 'Pseudo2D_allemitting'
cut = '30deg' # or normal
run = ''
duration = 120
Restrict = {'E0EtaLim':(1,J2MeV)} #eliminate all electrons whos E [MeV] * eta < 1 MeV photon
#Restrict = {'XLim':(-10,10,1e6),'YLim':(-10,10,1e6)} 
EnergyUnits = 'MeV'
param_str = ''
for key in Restrict:
    if key != 'E0EtaLim':
        param_str = param_str + f'{key}_{Restrict[key][0]}_{Restrict[key][1]:.3f}_'
    else:
        param_str = param_str + f'{key}_{Restrict[key][0]}_{EnergyUnits}_'


params = {'phi':(-180,180,181), 'E0':(0,900,651,J2MeV)}

pltName = f'{target}_{cut}_{pltType}_{param_str}_{run}'+'.gif'

#Restrict = {'IDLim':(np.zeros(1,100))}
#/Volumes/Elements/targets_eta/
ED = EDC.ElectronDataClass({'Path':'/Volumes/4TBWD/PhotonProject/eta/' + f'{target}' + '/' + f'{cut}' + '/' + f'{run}' + '/'})
fig = plt.figure(figsize=(20,10))

#%%

def make_frame(ED,*args):
    fig = args[0]
    fname = args[1]
    plotName = args[2]
    n = args[3]
    Restrict = args[4]
    params = args[5]
    t = ED.LoadElectronsAtTime(f'{fname}',n,Restrict)
    # HistV, Histx, Histy  = ED.Pseudo2D(params)
    # X,Y = np.meshgrid(Histx,Histy)
    ax0 = fig.gca()
    HistV, Histx, Histy  = ED.Pseudo2D(params)
    X,Y = np.meshgrid(Histx,Histy)
    logHist = np.log(HistV)    
    im = ax0.pcolormesh(X,Y,np.transpose(logHist),vmin=np.nanmin(logHist),vmax=np.nanmax(logHist), cmap=plt.get_cmap('jet'))
    im.set_clim([9,15])
    xtickmarks = np.array([-180,-150,-120,-90,-60,-30,0,30,60,90,120,150,180])
    ax0.grid(color='grey', linestyle='-', linewidth=.5)
    ce = fig.colorbar(im,ax=ax0)    
    ax0.set_xlabel(r'$\Phi [\circ]$',fontsize=16)
    ax0.set_ylabel(r'$E [MeV]$',fontsize=16)
    ce.set_label(r'electron $\frac{d^2 logN}{d\phi dE}$',fontsize=16)
    ax0.set_xticks(xtickmarks)
    ax0.axvline(color="black",linewidth=3)
    # logHist = np.log(HistV)
    # im = ax0.pcolormesh(X,Y,np.transpose(logHist),vmin=np.nanmin(logHist),vmax=np.nanmax(logHist), cmap=plt.get_cmap('jet'))
    # im.set_clim([9,13])
    #ax0.grid(color='k', linestyle='-', linewidth=.5)
    #LineOut,phiMidPoints = ED.LineOut('phi','E0',J2MeV)
    #ax0.plot(phiMidPoints,LineOut)
    ax0.set_title(fr'{plotName} ' +  fr'{t * 1e15:.2f} fs',fontsize=22)
    ax0.tick_params(labelsize=16)
    ax0.tick_params(axis='x',rotation=-45)

    #ax0.set_xlabel(r'$\Phi [\circ]$',fontsize=16)
    #ax0.set_ylabel(r'$dEd\Phi$',fontsize=16)
    #ax0.set_xticks(np.linspace(-180,180,21))
    #ax0.set_ylim([0, 1])
    #ax.grid()#


moviewriter = ani.PillowWriter(fps=20)

with moviewriter.saving(fig,pltName,duration):
    for j in range(duration):
        make_frame(ED,fig,'e_loc',f'ang_elec_{target}_{cut}',j,Restrict,params)
        moviewriter.grab_frame()
        fig.clf()
    moviewriter.finish()

