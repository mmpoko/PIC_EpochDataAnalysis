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
import time

MeV = (1.6022e-19)*(1.0e6) 
c_speed = 299792458 # m/s
mc = 9.11e-31 * c_speed #kg m /s
J2MeV = 1 * (6.242e18) / (1.0e6)

target = 'uniform' #uniform or channel
cut = 'normal2' # or normal
duration = 4
#Restricte = {'XLim':(-5,50,1e6),'YLim':(-15,15,1e6),'E0Lim':(20,1000,J2MeV)}
Restricte = {'XLim':(-5,50,1e6),'YLim':(-15,15,1e6),'E0Lim':(0,1000,J2MeV)}
Restrictp = {'E0Lim':(1,300,J2MeV)}
paramse = {'phi':(-180,180,181), 'E0':(0,900,301,J2MeV)}
paramsp = {'phi':(-180,180,181), 'E0':(1,300,101,J2MeV)}
#Restrict = {'IDLim':(np.zeros(1,100))}
ED = EDC.ElectronDataClass({'Path':'/Volumes/4TBWD/PhotonProject/track_carbon/' + f'{target}' + '/' + f'{cut}' + '/'})
Ph = PDC.PhotonDataClass({'Path':'/Volumes/4TBWD/PhotonProject/track_carbon/' + f'{target}' + '/' + f'{cut}' + '/'})

fig = plt.figure()
#fig, axs = plt.subplots(1,2)
#axs = fig.gca()
#%%

def make_frame(ED,Ph,*args):
    
    fig = args[0]
    #axs = args[1]
    axs = fig.gca()
    plotName = args[1]
    n = args[2]
    ename = args[3]
    Restricte = args[4]
    paramse = args[5]
    pname = args[6]
    Restrictp = args[7]
    paramsp = args[8]
    
    t = ED.LoadElectronsAtTime(f'{ename}',n,Restricte)
    Ph.KeepOnlyEmittingParticles(ED)
    #t2 = Ph.PositionExactRestrict(ED)
    fig.suptitle(fr'{plotName} ' +  fr'{t * 1e15:.2f} fs',fontsize=22)
    
    #params = {'phi':(-180,180,181), 'E0':(20,600,101,J2MeV)}
    HistV, Histx, Histy  = ED.Pseudo2D(paramse)
    #HistVp, Histxp, Histyp = Ph.Pseudo2D(paramsp)

    X,Y = np.meshgrid(Histx,Histy)
    #Xp,Yp = np.meshgrid(Histxp,Histyp)
    logHist = np.log(HistV)
    #logHistp = np.log(HistVp)
    
    im = axs.pcolormesh(X,Y,np.transpose(logHist),vmin=np.nanmin(logHist),vmax=np.nanmax(logHist), cmap=plt.get_cmap('jet'))
    im.set_clim([9,13])
    axs.grid(color='k', linestyle='-', linewidth=.5)
    c = fig.colorbar(im,ax=axs)    
    axs.set_xlabel(r'$\Phi [\circ]$',fontsize=16)
    axs.set_ylabel(r'$E [MeV]$',fontsize=16)
    
    # imp = axs[1].pcolormesh(Xp,Yp,np.transpose(logHistp),vmin=np.nanmin(logHistp),vmax=np.nanmax(logHistp), cmap=plt.get_cmap('jet'))
    # imp.set_clim([9,13])
    # axs[1].grid(color='k', linestyle='-', linewidth=.5)
    # cp = fig.colorbar(im,ax=axs[0])    
    # axs[1].set_xlabel(r'$\Phi [\circ]$',fontsize=16)
    # axs[1].set_ylabel(r'$E [MeV]$',fontsize=16)
    
    
    #ax.set_xticks(np.linspace(-180,180,4))
    #ax.set_ylim([0, 1e14])
    #ax.grid()#

#%%

metadata = dict(title=f'elec_ang_{target}_{cut}', artist='Matplotlib',comment='Movie support!')
ffmpegWriter = ani.writers['ffmpeg']
moviewriter = ffmpegWriter(fps=15,metadata=metadata)
param_stre = '_e_'
param_strp = '_p_'

for key in Restricte:
    param_stre = param_stre + f'{key}_{Restricte[key][0]}_{Restricte[key][1]}_'
for key in Restrictp:
    param_strp = param_strp + f'{key}_{Restrictp[key][0]}_{Restrictp[key][1]}_'
param_str = param_stre + param_strp


#load photons at end tiem to get acuumulated photons
pname = 'PhotonDetail'
n = 400
Restrictp = {}
t2 = Ph.LoadPhotonsAtTime(f'{pname}',n,Restrictp)

#keep only emitting electrons
#t = ED.LoadElectronsAtTime('id',350,Restricte)
#Ph.KeepOnlyEmittingParticles(ED)

with moviewriter.saving(fig, f'elec_ang_2D{target}_{cut}_{param_str}.mp4',duration):
    for j in range(duration):
        make_frame(ED,Ph,fig,f'{target} ' + f'{cut} ',j,'id',Restricte,paramse,'PhotonDetail',Restrictp,paramsp)
        moviewriter.grab_frame()
        fig.clf()
    moviewriter.finish()

#%%
#ax.show()
#eMeV = ED.E0[ED.eleInd > 0] / (1e6) * 6.242e18 
#eMeV = ED.E0[ED.eleInd > 0] * ED.w[ED.eleInd > 0] / (1e6) * 6.242e18