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


def makesubplotFrame(fig,axn,pltType,ED,Ph,FC,Dens,args1,args2,txtargs):
    axs=fig.get_axes()
    cax=axs[axn]
    if pltType=='photon' or pltType =='electron':
        t = ED.LoadElectronsAtTime(args1['fileName'],args1['time'],args1['Restrict'])
        t2 = Ph.LoadPhotonsAtTime(args2['fileName'],args2['time'],args2['Restrict'])
        Ph.KeepOnlyEmittingParticles(ED)
        if pltType=='electron':
            HistV, Histx, Histy  = ED.Pseudo2D(args1['params'])
            X,Y = np.meshgrid(Histx,Histy)
            logHist = np.log(HistV)    
            del HistV
            im = cax.pcolormesh(X,Y,np.transpose(logHist),vmin=args1['cbarlims'][0],vmax=args1['cbarlims'][-1], cmap=plt.get_cmap('jet'))
            del X,Y,logHist
            cax.grid(color='grey', linestyle='-', linewidth=.5)
            if args1['cbarOn'] == True:
                ce = fig.colorbar(im,ax=cax,shrink=txtargs[-3],orientation='horizontal',aspect=txtargs[-1])    
                ce.set_label(args1['cbarLbl'],fontsize=txtargs[1])
                ce.set_ticks(args1['cbarlims'])
                ce.ax.tick_params(labelsize=txtargs[0])
            if args1['XAxLbl'] == True:
                cax.set_xlabel(args1['Xlbl'])
            else:
                cax.xaxis.set_ticklabels([])
            if args1['YAxLbl']==True:    
                cax.set_ylabel(args1['Ylbl'])
            else:
                cax.yaxis.set_ticklabels([])
                
            cax.axvline(color="grey",linewidth=2)
            cax.set_xticks(args1['xlims'])
            if args1['LastYTick'] == True:
                cax.set_yticks(args1['ylims'])
            else:
                cax.set_yticks(args1['ylims'][:-1])  
            cax.set_xlim([args1['xlims'][0],args1['xlims'][-1]])
            cax.set_ylim([args1['ylims'][0],args1['ylims'][-1]])
            return fig
        elif pltType=='photon':
            if hasattr(Ph,'peid'):
                if Ph.peid.size !=0:
                    HistV, Histx, Histy  = Ph.Pseudo2D(args2['params'])
                    X,Y = np.meshgrid(Histx,Histy)
                    logHist = np.log(HistV)    
                    del HistV
                    im = cax.pcolormesh(X,Y,np.transpose(logHist),vmin=args2['cbarlims'][0],vmax=args2['cbarlims'][-1], cmap=plt.get_cmap('jet'))
                    del X,Y,logHist
                    cax.grid(color='grey', linestyle='-', linewidth=.5)
                    if args2['cbarOn'] == True:
                        ce = fig.colorbar(im,ax=cax,shrink=txtargs[-3],orientation='horizontal',aspect=txtargs[-1])    
                        ce.set_label(args2['cbarLbl'],fontsize=txtargs[1])
                        ce.set_ticks(args2['cbarlims'])
                        ce.ax.tick_params(labelsize=txtargs[0])
                    if args2['XAxLbl'] == True:
                        cax.set_xlabel(args2['Xlbl'])
                    else:
                        cax.xaxis.set_ticklabels([])
                    if args2['YAxLbl']==True:    
                        cax.set_ylabel(args2['Ylbl'])
                    else:
                        cax.yaxis.set_ticklabels([])
                    cax.axvline(color="grey",linewidth=2)
                    cax.set_xticks(args2['xlims'])
                    if args2['LastYTick'] == True:
                        cax.set_yticks(args2['ylims'])
                    else:
                        cax.set_yticks(args2['ylims'][:-1])  
                    cax.set_xlim([args2['xlims'][0],args2['xlims'][-1]])
                    cax.set_ylim([args2['ylims'][0],args2['ylims'][-1]])
                    return fig
            return fig
    else:
        subSampInt=1
        t=FC.LoadFieldAtTime(args1['fileName'],args1['time'],args1['field'])
        #normalize field
        FC.ey=FC.ey/args1['norm']
        t=Dens.LoadDensityAtTime(args2['fileName'],args2['time'],args2['species'],args2['norm'])
        cax.grid(color='grey', linestyle='-', linewidth=.5) 
        if args2['XAxLbl'] == True:
            cax.set_xlabel(args2['Xlbl'])
        else:
            cax.xaxis.set_ticklabels([])
        if args2['YAxLbl']==True:    
            cax.set_ylabel(args2['Ylbl'])
        else:
            cax.yaxis.set_ticklabels([])
        xd,yd = np.meshgrid(Dens.x[::subSampInt]*1e6,Dens.y[::subSampInt]*1e6)
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
        den = Dens.dens[::subSampInt,::subSampInt]
        m1.set_array(np.transpose(den))
        m1.set_clim([args2['cbarlims'][0],args2['cbarlims'][-1]])
        imd = cax.pcolormesh(xd,yd,np.transpose(den),vmin=args2['cbarlims'][0],vmax=args2['cbarlims'][-1],cmap=plt.get_cmap('jet_alpha'))
        del den,xd,yd
        if args2['cbarOn']==True:
            cd = fig.colorbar(m1,ax=cax,shrink=txtargs[-4],aspect=txtargs[-2])
            cd.set_label(args2['cbarLbl'],fontsize=txtargs[1])
            cd.set_ticks(args2['cbarlims'])
            cd.ax.tick_params(labelsize=txtargs[0])
        # get colormap
        ncolors = 256
        color_array = plt.get_cmap('Greys')(range(ncolors))
        color_array[:,-1] = np.linspace(0.0,1.0,ncolors)
        #color_array[64:,-1] = 1
        map_object = LinearSegmentedColormap.from_list(name='greys_alpha',colors=color_array)
        plt.register_cmap(cmap=map_object)    
        absey = abs(FC.ey[::subSampInt,::subSampInt])
        xd,yd = np.meshgrid(FC.x[::subSampInt] * 1e6,FC.y[::subSampInt] * 1e6)
        imf = cax.contour(xd,yd,np.transpose(absey),vmin=args1['cbarlims'][0],vmax=args1['cbarlims'][-1], cmap=plt.get_cmap('greys_alpha'))
        del xd,yd
        m = plt.cm.ScalarMappable(cmap=plt.get_cmap('greys_alpha'))
        m.set_array(np.transpose(absey))
        del absey
        m.set_clim([args1['cbarlims'][0],args1['cbarlims'][-1]])
        if args1['cbarOn'] == True:
            cf = fig.colorbar(m,ax=cax,shrink=txtargs[-3],orientation='horizontal',aspect=txtargs[-1])
            cf.set_label(args1['cbarLbl'],fontsize=txtargs[1])
            cf.set_ticks(args1['cbarlims'])
            cf.ax.tick_params(labelsize=txtargs[0])
            #cf.ax.xaxis.set_ticks_position('top')
            #cf.ax.xaxis.set_label_position('top')
        cax.grid(color='grey', linestyle='-', linewidth=.5)
        cax.set_xticks(args1['xlims'])
        if args1['LastYTick'] == True:
            if args1['FirstYTick'] == True:
                cax.set_yticks(args1['ylims'])
            else:
                cax.set_yticks(args1['ylims'][1:])
        else:
            if args1['FirstYTick']==True:
                cax.set_yticks(args1['ylims'][:-1])  
            else:
                cax.set_yticks(args1['ylims'][1:-1]) 
                
        cax.set_xlim([args1['xlims'][0],args1['xlims'][-1]])
        cax.set_ylim([args1['ylims'][0],args1['ylims'][-1]])
        return fig


    #         LineOutp,phiMidPoints = ED.LineOut('phi','E0',J2MeV)
    #         axs[1].plot(phiMidPoints,LineOutp)
            
    # axs[1].grid(color='k', linestyle='-', linewidth=.5) 
    # axs[1].set_xlabel(r'$\Phi [\circ]$',fontsize=16)
    # axs[1].set_ylabel(r'$dEd\Phi$',fontsize=16)
    # axs[1].set_xticks(np.linspace(-180,180,11))
    # axs[1].set_ylim([0, .01])
    