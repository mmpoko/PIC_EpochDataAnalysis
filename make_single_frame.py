#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  5 13:20:40 2021

@author: michaelpokornik
"""
import matplotlib.pyplot as plt
import numpy as np

def make_single_frame(ED,Ph,*args):
    
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
    cbarOn = args[10]
    x1lims = args[11]
    y1lims = args[12]
    x2lims = args[13]
    y2lims = args[14]
    cbar1lims = args[15]
    cbar2lims = args[16]
    SMALL_SIZE = args[17]
    MEDIUM_SIZE = args[18]
    BIGGER_SIZE = args[19]
    t = ED.LoadElectronsAtTime(f'{ename}',n,Restricte)
    t2 = Ph.LoadPhotonsAtTime(f'{pname}',2*n,Restrictp)
    Ph.KeepOnlyEmittingParticles(ED)
    

    #fig.suptitle(fr'{plotName} ' +  fr'{t * 1e15:.2f} fs',fontsize=22)
    
    #%% if pseudocolor
    HistV, Histx, Histy  = ED.Pseudo2D(paramse)
    X,Y = np.meshgrid(Histx,Histy)
    logHist = np.log(HistV)    
    im = axs[0].pcolormesh(X,Y,np.transpose(logHist),vmin=np.nanmin(logHist),vmax=np.nanmax(logHist), cmap=plt.get_cmap('jet'))
    im.set_clim(cbar1lims)
    axs[0].grid(color='grey', linestyle='-', linewidth=.5)
    if cbarOn == True:
        ce = fig.colorbar(im,ax=axs[0])    
        ce.set_label(r'electron $\frac{d^2 logN}{d\phi dE}$',fontsize=MEDIUM_SIZE)
        ce.ax.tick_params(labelsize=SMALL_SIZE) 
        #ce.set_clim(cbar1lims)
    axs[0].set_xlabel(r'$\Phi [\circ]$',fontsize=MEDIUM_SIZE)
    axs[0].set_ylabel(r'$E [MeV]$',fontsize=MEDIUM_SIZE)
    axs[0].axvline(color="grey",linewidth=2)
    xtickmarks = np.array([-180,-150,-120,-90,-60,-30,0,30,60,90,120,150,180])
    axs[0].set_xticks(xtickmarks)
    axs[0].tick_params(labelsize=SMALL_SIZE)
    axs[0].tick_params(axis='x',rotation=-45)
    axs[0].set_xlim(x1lims)
    axs[0].set_ylim(y1lims)
    if hasattr(Ph,'peid'):
        if Ph.peid.size !=0:
            HistVp, Histxp, Histyp = Ph.Pseudo2D(paramsp)
            Xp,Yp = np.meshgrid(Histxp,Histyp)
            logHistp = np.log(HistVp)
            imp = axs[1].pcolormesh(Xp,Yp,np.transpose(logHistp),vmin=np.nanmin(logHistp),vmax=np.nanmax(logHistp), cmap=plt.get_cmap('jet'))
            if cbarOn == True:
                cp = fig.colorbar(imp,ax=axs[1]) 
                cp.set_label(r'photon $\frac{d^2 logN}{d\phi dE}$',fontsize=MEDIUM_SIZE)
                cp.ax.tick_params(labelsize=SMALL_SIZE) 
                #cp.set_clim(cbar2lims)
            imp.set_clim(cbar2lims)
            
    
    axs[1].grid(color='grey', linestyle='-', linewidth=.5) 
    axs[1].axvline(color="grey",linewidth=2)
    axs[1].set_xlabel(r'$\Phi [\circ]$',fontsize=MEDIUM_SIZE)
    axs[1].set_ylabel(r'$E [MeV]$',fontsize=MEDIUM_SIZE)
    axs[1].set_xticks(xtickmarks)
    axs[1].tick_params(labelsize=SMALL_SIZE)
    axs[1].tick_params(axis='x',rotation=-45)
    axs[1].set_xlim(x2lims)
    axs[1].set_ylim(y2lims)
    fig.tight_layout(pad=1.0)
    return fig
    #fig.savefig("name.svg")


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
    