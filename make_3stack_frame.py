#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 14 10:07:22 2021

@author: michaelpokornik
"""
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  5 13:20:40 2021

@author: michaelpokornik
"""
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap
from mpl_toolkits.axes_grid1 import make_axes_locatable


def make_3stack_frame(ED,Ph,FC,noLbl,*args):
    
    fig = args[0]
    axs = fig.get_axes()
    #axs = fig.gca()
    plotName = args[1]
    n = args[2]
    ename = args[3]
    Restricte = args[4]
    paramse = args[5]
    pname = args[6]
    Restrictp = args[7]
    paramsp = args[8]
    cbarOn = args[9]
    x1lims = args[10]
    y1lims = args[11]
    cbar1lims = args[12]
    x2lims = args[13]
    y2lims = args[14]
    cbar2lims = args[15]
    x3lims = args[16]
    y3lims = args[17]
    cbar3lims = args[18]
    fnorm = args[19]
    fname = args[20]
    cbar3name = args[21]
    n2 = args[22]
    fvar = args[23]
    x3ticks = args[24]
    y3ticks = args[25]
    SMALL_SIZE = args[26]
    MEDIUM_SIZE = args[27]
    BIGGER_SIZE = args[28]
    t = ED.LoadElectronsAtTime(f'{ename}',n,Restricte)
    t2 = Ph.LoadPhotonsAtTime(f'{pname}',n,Restrictp)
    t3 = FC.LoadFieldAtTime(f'{fname}',n2,fvar)
    #normalize field
    FC.ey=FC.ey/fnorm
    Ph.KeepOnlyEmittingParticles(ED)
    subSampIntd=1
    subSampIntf = 1



    #fig.suptitle(fr'{plotName} ' +  fr'{t * 1e15:.2f} fs',fontsize=22)
    
    #%% if pseudocolor
    #electron
    HistV, Histx, Histy  = ED.Pseudo2D(paramse)
    X,Y = np.meshgrid(Histx,Histy)
    logHist = np.log(HistV)    
    del HistV
    im = axs[0].pcolormesh(X,Y,np.transpose(logHist),vmin=cbar1lims[0],vmax=cbar1lims[-1], cmap=plt.get_cmap('jet'))
    del X,Y,logHist
    axs[0].grid(color='grey', linestyle='-', linewidth=.5)
    if cbarOn == True:
        ce = fig.colorbar(im,ax=axs[0],shrink=.3)    
        ce.set_label(r'electron $\frac{d^2 logN}{d\theta dE}$',fontsize=MEDIUM_SIZE)
        ce.ax.tick_params(labelsize=SMALL_SIZE) 
        ce.set_ticks(cbar1lims)
        #ce.set_clim(cbar1lims)
    if noLbl == False:
        axs[0].set_xlabel(r'$\theta [^{\circ}]$',fontsize=MEDIUM_SIZE)
        axs[0].set_ylabel(r'Electron $E [MeV]$',fontsize=MEDIUM_SIZE)
    axs[0].axvline(color="grey",linewidth=2)
    axs[0].set_xticks(x1lims)
    axs[0].set_yticks(y1lims)
    axs[0].tick_params(labelsize=SMALL_SIZE)
    #axs[0].tick_params(axis='x',rotation=-45)
    axs[0].set_xlim([x1lims[0],x1lims[-1]])
    axs[0].set_ylim([y1lims[0],y1lims[-1]])
    #photon
    if hasattr(Ph,'peid'):
        if Ph.peid.size !=0:
            HistVp, Histxp, Histyp = Ph.Pseudo2D(paramsp)
            Xp,Yp = np.meshgrid(Histxp,Histyp)
            logHistp = np.log(HistVp)
            del HistVp
            imp = axs[1].pcolormesh(Xp,Yp,np.transpose(logHistp),vmin=cbar2lims[0],vmax=cbar2lims[-1], cmap=plt.get_cmap('jet'))
            del Xp,Yp,logHistp
            if cbarOn == True:
                cp = fig.colorbar(imp,ax=axs[1],shrink=.3) 
                cp.set_label(r'photon $\frac{d^2 logN}{d\theta dE}$',fontsize=MEDIUM_SIZE)
                cp.ax.tick_params(labelsize=SMALL_SIZE) 
                cp.set_ticks(cbar2lims)
            
            
    
    axs[1].grid(color='grey', linestyle='-', linewidth=.5) 
    axs[1].axvline(color="grey",linewidth=2)
    if noLbl == False:
        axs[1].set_xlabel(r'$\theta [^{\circ}]$',fontsize=MEDIUM_SIZE)
        axs[1].set_ylabel(r'Photon $E [MeV]$',fontsize=MEDIUM_SIZE)
    axs[1].set_xticks(x2lims)
    axs[1].set_yticks(y2lims)
    axs[1].tick_params(labelsize=SMALL_SIZE)
    #axs[1].tick_params(axis='x',rotation=-45)
    axs[1].set_xlim([x2lims[0],x2lims[-1]])
    axs[1].set_ylim([y2lims[0],y2lims[-1]])
    #field
    if noLbl == False:
        axs[2].set_xlabel(r'$x [\mu m]$',fontsize=MEDIUM_SIZE)
        axs[2].set_ylabel(r'$y [\mu m]$',fontsize=MEDIUM_SIZE)



    
    if len(args) > 29:
        Dens = args[29]
        species = args[30]
        normF = args[31]
        cbardenslims = args[32]
        cbarOn2 = args[33]
        t4=Dens.LoadDensityAtTime('dens',n2,species,normF)
        xd,yd = np.meshgrid(Dens.x[::subSampIntd]*1e6,Dens.y[::subSampIntd]*1e6)
        # get colormap
        ncolors = 256
        color_array = plt.get_cmap('jet')(range(ncolors))
    
        # change alpha values
        color_array[:,-1] = np.linspace(0.0,1.0,ncolors)
        # create a colormap object
        map_object = LinearSegmentedColormap.from_list(name='jet_alpha',colors=color_array)
    
        # register this new colormap with matplotlib
        plt.register_cmap(cmap=map_object)
        m1 = plt.cm.ScalarMappable(cmap=plt.get_cmap('jet_alpha'))
        den = Dens.dens[::subSampIntd,::subSampIntd]
        m1.set_array(np.transpose(den))
        m1.set_clim([cbardenslims[0],cbardenslims[-1]])
        imd = axs[2].pcolormesh(xd,yd,np.transpose(den),vmin=cbardenslims[0],vmax=cbardenslims[-1],cmap=plt.get_cmap('jet_alpha'))
        
        del den
        if cbarOn2==True:
            cd = fig.colorbar(m1,ax=axs[2],orientation='horizontal',pad=.5,shrink=.3)
            cd.set_label(r'$n/n_{cr}$',fontsize=MEDIUM_SIZE)
            cd.set_ticks(cbardenslims)
            #cd.set_clim([0,100])
            cd.ax.tick_params(labelsize=SMALL_SIZE)
    else:
        xd,yd = np.meshgrid(FC.x[::subSampIntf] * 1e6,FC.y[::subSampIntf] * 1e6)
    
    # get colormap
    ncolors = 256
    color_array = plt.get_cmap('Greys')(range(ncolors))
    color_array[:,-1] = np.linspace(0.0,1.0,ncolors)
    color_array[64:,-1] = 1
    #color_array[:,-1] = np.append(np.linspace(1.0,0.0,(ncolors//2)+1)[:-1],np.linspace(0.0,1.0,(ncolors//2)+1)[1:],axis=None)
        # change alpha values
    #color_array[64:,-1] = 1  
        # create a colormap object
    map_object = LinearSegmentedColormap.from_list(name='greys_alpha',colors=color_array)
    plt.register_cmap(cmap=map_object)
   
    
    #absey = abs(FC.ey[::subSampIntf,::subSampIntf])
    absey = FC.ey[::subSampIntf,::subSampIntf]
    
    xd,yd = np.meshgrid(FC.x[::subSampIntf] * 1e6,FC.y[::subSampIntf] * 1e6)
    imf = axs[2].contour(xd,yd,np.transpose(absey),vmin=cbar3lims[0],vmax=cbar3lims[-1], cmap=plt.get_cmap('greys_alpha'))
    del xd,yd
    imf.set_clim(cbar3lims)
    m = plt.cm.ScalarMappable(cmap=plt.get_cmap('greys_alpha'))
    m.set_array(np.transpose(absey))
    del absey
    m.set_clim(cbar3lims)
    axs[2].tick_params(labelsize=SMALL_SIZE)
    if cbarOn == True:
        cf = fig.colorbar(m,ax=axs[2],shrink=.3)
        cf.set_label(cbar3name,fontsize=SMALL_SIZE)
        cf.ax.tick_params(labelsize=SMALL_SIZE)
        cf.set_ticks(cbar3lims)
        #cf.set_clim([-1,1])
    axs[2].grid(color='grey', linestyle='-', linewidth=.5)
    axs[2].set_xticks(x3ticks)
    axs[2].set_yticks(y3ticks)
    axs[2].set_xlim([x3ticks[0],x3ticks[-1]])
    axs[2].set_ylim([y3ticks[0],y3ticks[-1]])
    fig.tight_layout()
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
    