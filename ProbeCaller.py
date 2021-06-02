from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import sdf_helper as sdf_helper
import ProbeDataClass as PDC
import numpy as np
from scipy import stats
from matplotlib.colors import LinearSegmentedColormap
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.colors as col

MeV = (1.6022e-19)*(1.0e6) 
c_speed = 299792458 # m/s
mc = 9.11e-31 * c_speed #kg m /s
J2MeV = 1 * (6.242e18) / (1.0e6)

target = 'targets_3d' #uniform or channel
probeFileName = 'probe'
probeName = 'tube_foil_probe'
run = 'run1'
duration = 120
n=28

#Restrict = {'XLim':(0,np.nan,1e6),'YLim':(0,np.nan,1e6)} 
Restrict = {} 
EnergyUnits = 'MeV'





pltName = f'{target}_{run}_{probeName}'+'.gif'

#Restrict = {'IDLim':(np.zeros(1,100))}
#/Volumes/Elements/targets_eta/
PD = PDC.ProbeDataClass({'Path':'/Volumes/Elements/' + f'{target}' + '/' + f'{run}' + '/','probeName':probeName})
PD.LoadProbeAtTime(probeFileName,n,Restrict)
fig = plt.figure(figsize=(20,10))
nphibin = 181
phiBins = np.linspace(-90,90,num=nphibin)
phiMid = (phiBins[1:] + phiBins[0:-1]) /2
# thetaBins = np.linspace(-180,180,num=361)  # for forward and backward ,total sphere
nthetabin = 181
thetaBins = np.linspace(-90,90,num=nthetabin)
thetaMid = (thetaBins[1:] + thetaBins[0:-1]) /2
#calc after max gamma
ngambin = 51
gamBins = np.linspace(1,np.max(PD.gam),num=ngambin)
gamMid = (gamBins[1:] + gamBins[0:-1]) /2
phiInds = np.digitize(PD.phi,phiBins,right=True)-1
thetaInds = np.digitize(PD.theta,thetaBins,right=True)-1
gamInds = np.digitize(PD.gam,gamBins,right=True)-1


w_res = np.zeros(PD.w.shape)
w_res[PD.eleInd>0]=PD.w[PD.eleInd>0]
#ret = stats.binned_statistic_dd([gamInds,thetaInds,phiInds],PD.w, bins=[gamBins,thetaBins,phiBins],statistic='sum')
ret = stats.binned_statistic_dd([gamInds,thetaInds,phiInds],w_res, bins=[ngambin-1,nthetabin-1,nphibin-1],statistic='sum')
bincounts = ret.statistic
pdf = bincounts / np.sum(bincounts) / np.max(np.diff(gamBins)) / np.max(np.diff(phiBins)) / np.max(np.diff(thetaBins))
pdf = np.transpose(pdf)

#%%
fig = plt.figure()
ax1 = fig.add_subplot(111, projection='3d')
#ax1 = fig.add_subplot(111)

#Use one less than bin edges to give rough bin location
X, Y = np.meshgrid(thetaMid,phiMid)
ct = []
pdfmask = np.ma.array(pdf,mask=pdf==0.0)
for i in range(ngambin-1):
    if np.nansum(pdf[:,:,i]) > 0:
        ct.append(i)
ctsub= ct[::np.int(np.ceil(len(ct)/5))]
#ctsub = 

#Loop over range of slice locations (default histogram uses 10 bins)
##
ncolors = 256
color_array = plt.get_cmap('jet')(range(ncolors))
color_array[:,-1] = np.linspace(0,.5,num=256)
map_object = LinearSegmentedColormap.from_list(name='jet_alpha',colors=color_array)
plt.register_cmap(cmap=map_object)    
m = plt.cm.ScalarMappable(cmap=plt.get_cmap('jet_alpha'))
m.set_clim(0,np.max(pdf))

for ii in ctsub:
     cs = ax1.contourf(X,Y,pdf[:,:,ii], 
                       zdir='z', 
                       offset=gamMid[ii], 
                       level=100, 
                       cmap=plt.get_cmap('jet_alpha')
                       )

cf = fig.colorbar(m,ax=[ax1],shrink=.2,location='right')
cf.set_label(r'$\frac{d^{3}N}{d\theta d\phi d\gamma}$',fontsize=20)
     
ax1.set_xlabel('theta')
ax1.set_ylabel('phi')
ax1.set_zlabel(r'$\gamma$')
ax1.set_xlim(-90, 90)
ax1.set_ylim(-90, 90)
ax1.set_zlim(gamMid[0],gamMid[-1])
ax1.view_init(20,30)
#plt.colorbar(cs)
plt.show()
#ax1.view_init(30,30)
#fig