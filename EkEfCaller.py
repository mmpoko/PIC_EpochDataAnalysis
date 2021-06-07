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
import FieldClass as FC

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

def CrossSections(CS,data,dname,ind):
    if CS['sliceD'] == 'x':
        dat = eval(f'data.{dname}[ind,:,:]')
    elif CS['sliceD'] == 'y':
        dat = eval(f'data.{dname}[:,ind,:]')
    else:
        dat = eval(f'data.{dname}[:,:,ind]')
    return dat

MeV = (1.6022e-19)*(1.0e6) 
c_speed = 299792458 # m/s
mc = 9.11e-31 * c_speed #kg m /s
J2MeV = 1 * (6.242e18) / (1.0e6)
eps0 = 8.854e-12 #F/m

target = 'targets_3d' #uniform or channel
EfileName = 'edens'
species = 'Derived_Average_Particle_Energy_averaged_Electron_t'
run = 'run30cpm'
duration = 120
n=25

EnergyUnits = 'MeV'
EnC = DC.DensityClass({'Path':'/Volumes/Elements/' + f'{target}' + '/' + f'{run}' + '/'})
t = EnC.LoadEkbarAtTime(EfileName,n,species,(1/J2MeV))
#fields
I0 = 1e21
E0 = np.sqrt(2*I0/eps0/c_speed)
#FC = FC.FieldClass({'Path':'/Volumes/Elements/' + f'{target}' + '/' + f'{run}' + '/','dim':3})
FC = FC.FieldClass({'Path':EnC.Path,'dim':3})
efield = 'Electric_Field_Ex'
t2 = FC.LoadFieldAtTime('e_fields',n,efield)

#Energy Contours
cbarV = [0,15,30]
nlevels = 20
clevels = np.linspace(cbarV[0],cbarV[-1],num=nlevels)
#slices of the 3d
CrossSection = {'pltx':EnC.z,'plty':EnC.y,'sliceD':'x'}
xvalues = np.linspace(-3,3,num=7)
xvalues = [-2.5,-2,0,2,2.5]
zvalues = np.linspace(-3,6,num=10)
zvalues = [-10,-5,0,5,10]
yval = [2e-6,1e-6,0]
subt = ['y=2um','y=1um','y=0um']
ln=len(subt)
xg,yg = np.meshgrid(EnC.z*1e6,EnC.y*1e6)

#Field Contours
cbarfV = [-3,0,3]
nlevels = 20
clevelsf = np.linspace(cbarfV[0],cbarfV[-1],num=nlevels)
#slices of the 3d
CrossSectionf = {'pltx':FC.z,'plty':FC.x,'sliceD':'y'}
xfg,yfg = np.meshgrid(FC.z*1e6,FC.x*1e6)

#plt formatting
lblsz=36
ticksz=30
cbarft=30
cbart=40

ind = {}
indf = {}

for i, val in enumerate(yval):
   ind[i] = find_nearest(EnC.x,val)
   indf[i] = find_nearest(FC.y,val)
   
#fig,ax = plt.subplots(1,n,figsize=(20,10))
fig,ax = plt.subplots(ln,1,figsize=(10,20))


for i in range(ln):
    dat = CrossSections(CrossSection, EnC,'Ek', ind[i]) #energy 
    datf = CrossSections(CrossSectionf, FC,'ex',  indf[i]) #field
    datf = datf/ E0
    #energy
    pltM = ax[i].contourf(xg,yg,dat,cmap='jet',vmin=cbarV[0],vmax=cbarV[-1],extend='both',alpha=1,levels=clevels)
    #field
    pltF = ax[i].contourf(xfg,yfg,datf,cmap='bwr_r',vmin=cbarfV[0],vmax=cbarfV[-1],extend='both',alpha=0,levels=clevelsf)
    ax[i].set_title(subt[i],color='black',fontsize=40,loc='left',x=0, y=0.8)
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
#cbar=fig.colorbar(pltF, cax=cbar_ax,orientation='horizontal')
cbar.set_ticks(cbarV)
#cbar.set_ticks(cbarfV)
cbar.ax.tick_params(labelsize=cbarft)
cbar.ax.set_title(r'$\overline{E_{k}}$ [MeV]',fontsize=cbart)
#cbar.ax.set_title(r'$E_{x}/E_{0}$',fontsize=cbart)
# plt.colorbar(m, boundaries=np.linspace(0, 2, 6))