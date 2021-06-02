#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
import sdf_helper as sdfhelper
import os
import glob

class ElectronDataClass:
    def __init__(self, *args,**kwargs):
        self.MeV = (1.6022e-19)*(1.0e6) 
        self.c_speed = 299792458 # m/s
        self.mc = 9.11e-31*self.c_speed #kg m /s
        self.pxread = False
        self.pyread = False
        for dictionary in args:
            for key in dictionary:
                setattr(self,key,dictionary[key])
        for key in kwargs:
            setattr(self,key,kwargs[key])
    
    def LoadElectronsAtTime(self,elFileName,fileNum, *restrict):
        elFilePath = f"{self.Path}" + f'{elFileName}' + f'{fileNum:04d}' + '.sdf'
        curFile = sdfhelper.getdata(fr'{elFilePath}')
        print(f"time = {curFile.Header['time'] * 1e15} fs")
        #load in electron ids first
        for field in dir(curFile):
            if 'Particles_ID' in field:
                self.eIds = eval(f'curFile.{field}.data')
                self.eleInd = np.ones(self.eIds.shape,dtype=int)
                break
        #load in rest of the data
        for field in dir(curFile):
            if 'Particles_ID' in field:
                continue
            elif 'Grid_Particles' in field:
                self.x = eval(f'curFile.{field}.data[0]')
                self.y = eval(f'curFile.{field}.data[1]')
                continue
            elif 'Particles_Px' in field:
                self.px = eval(f'curFile.{field}.data')
                self.pxread = True
                continue
            elif 'Particles_Py' in field:
                self.py = eval(f'curFile.{field}.data')
                self.pyread = True
                continue
            elif 'Particles_Weight' in field:
                self.w = eval(f'curFile.{field}.data')
                continue
            elif '_eta_' in field:
                self.eta = eval(f'curFile.{field}.data')
            elif 'Particles_Gamma' in field:
                self.gam = eval(f'curFile.{field}.data')
                self.E0 = self.gam * self.mc * self.c_speed
                continue
            elif self.pxread == True and self.pyread == True:
                p2 = (self.px **2) + (self.py ** 2)
                self.E0 = np.sqrt(1+ (p2/(self.mc ** 2)))* self.mc * self.c_speed
            else:
                continue
                
        for dictionary in restrict:
            for key in dictionary:
                if 'XLim' in key:
                    self.XLimRestrict(dictionary[key][0:-1],dictionary[key][-1])
                elif 'YLim' in key:
                    self.YLimRestrict(dictionary[key][0:-1],dictionary[key][-1])
                elif 'PxLim' in key:
                    self.PxLimRestrict(dictionary[key][0:-1],dictionary[key][-1])
                elif 'PyLim' in key:
                    self.PyLimRestrict(dictionary[key][0:-1],dictionary[key][-1])  
                elif 'GammaLim' in key:
                    self.GammaLimRestrict(dictionary[key][0:-1],dictionary[key][-1])
                elif 'E0Lim' in key:
                    self.E0LimRestrict(dictionary[key][0:-1],dictionary[key][-1])
                elif 'IDLim' in key:
                    self.IDRestrict(dictionary[key])
                elif 'eIDLim' in key:
                    self.IDexactRestrict(dictionary[key][0],dictionary[key][-1])
                elif 'E0EtaLim' in key:
                    self.EtaEnergyRestrict(dictionary[key][0],dictionary[key][-1])
                else:
                    continue
                
        return curFile.Header['time']
    
                
   # restrictions  
    def XLimRestrict(self,bounds,units):
        self.eleInd[(self.x * units < bounds[0]) & (np.isnan(bounds[0]) == False)]= 0
        self.eleInd[(self.x * units > bounds[1]) & (np.isnan(bounds[1]) == False)]= 0
    def YLimRestrict(self,bounds,units):
        self.eleInd[(self.y * units < bounds[0]) & (np.isnan(bounds[0]) == False)]= 0
        self.eleInd[(self.y * units > bounds[1]) & (np.isnan(bounds[1]) == False)]= 0 
    def PxLimRestrict(self,bounds,units):
        self.eleInd[(self.px * units < bounds[0]) & (np.isnan(bounds[0]) == False)]= 0
        self.eleInd[(self.px * units > bounds[1]) & (np.isnan(bounds[1]) == False)]= 0
    def PyLimRestrict(self,bounds,units):
        self.eleInd[(self.py * units < bounds[0]) & (np.isnan(bounds[0]) == False)]= 0
        self.eleInd[(self.py * units > bounds[1]) & (np.isnan(bounds[1]) == False)]= 0
    def GammaLimRestrict(self,bounds,units):
        self.eleInd[(self.gam * units < bounds[0]) & (np.isnan(bounds[0]) == False)]= 0
        self.eleInd[(self.gam * units > bounds[1]) & (np.isnan(bounds[1]) == False)]= 0     
    def E0LimRestrict(self,bounds,units):
        self.eleInd[(self.E0 * units < bounds[0]) & (np.isnan(bounds[0]) == False)]= 0
        self.eleInd[(self.E0 * units > bounds[1]) & (np.isnan(bounds[1]) == False)]= 0
    def IDRestrict(self,ids):
        self.eleInd[ids == 0] = 0
    def EtaEnergyRestrict(self,PhCutOff,units):
        self.eleInd[self.E0 * units * .44 * self.eta < PhCutOff]=0
    def IDExactRestrict(self,ids):
        self.eleInd[~np.isin(self.eIds,ids)] = 0
    def GetPositions(self):
        return self.x[self.eleInd > 0] , self.y[self.eleInd > 0]
    
        
    
    def LineOut(self,LineOutx,LineOuty,*units):
        normVal = 1e6 #take into account 3rd dim = 1 meter
        if LineOutx == "phi":
            rad_angle = np.arctan2(self.py[self.eleInd > 0],self.px[self.eleInd > 0])
            deg_angle = np.rad2deg(rad_angle)
            phi_bins = np.linspace(-180, 180,181)
            digitized = np.digitize(deg_angle,phi_bins)
        else:
            print('no other lineout type has been set up yet')
            
        if LineOuty == "E0":
            normCounts = self.w[self.eleInd > 0] * self.E0[self.eleInd > 0] * units / np.mean(np.diff(phi_bins)) / normVal
        else:
            normCounts = self.w[self.eleInd > 0] / np.mean(np.diff(phi_bins)) / normVal
            
        phiMidPoints = (phi_bins[1:] + phi_bins[:-1]) / 2
        LineOut = np.bincount(digitized-1, weights=normCounts)
        LineOutMissing = np.zeros((phiMidPoints.shape[0] - LineOut.shape[0]))
        LineOut = np.append(LineOut,LineOutMissing)
        return LineOut, phiMidPoints
    
    
    def Pseudo2D(self,*PseudoType):
        
        bins_ = {}
        digitized = {}
        normVal = 1e6 # will divide by 1 million to account for third dim = 1meter
        indC = -1
        for dictionary in PseudoType:
            for key in dictionary:
                indC +=1
                if 'phi' in key:
                    rad_angle = np.arctan2(self.py[self.eleInd > 0],self.px[self.eleInd > 0])
                    deg_angle = np.rad2deg(rad_angle)
                    bins_[indC] = np.linspace(dictionary[key][0],dictionary[key][1],dictionary[key][2])
                    digitized[indC] = np.digitize(deg_angle,bins_[indC])
                    normVal *= np.mean(np.diff(bins_[indC]))
                elif "E0" in key:
                    e_ = self.E0[self.eleInd > 0] * dictionary[key][3] #change energy to the right units
                    bins_[indC] = np.linspace(dictionary[key][0],dictionary[key][1],dictionary[key][2])
                    digitized[indC] = np.digitize(e_,bins_[indC])
                    normVal *= np.mean(np.diff(bins_[indC]))
        MatDimSum = np.int((bins_[0].size -1) * (bins_[1].size -1))
        normCounts = self.w[self.eleInd > 0]
        #shift down due to digitiize indexing
        digitized[0]-=1
        digitized[1]-=1
        #if values are outside bins, just assign to closest bin
        ax1 = digitized[0]
        ax2 = digitized[1]
        digitized[0][ax1 > bins_[0].size-2] = bins_[0].size - 2
        digitized[1][ax2 > bins_[1].size-2] = bins_[1].size - 2
        digitized[0][ax1 < 0] = 0 
        digitized[1][ax2 < 0] = 0
        rolledInd = np.ravel_multi_index(np.array([[digitized[0]],[digitized[1]]],dtype='int64'),(bins_[0].size -1,bins_[1].size -1))
        rolledSum = np.bincount(rolledInd[0],weights=normCounts/normVal)
        rolledSumMissing = np.zeros((MatDimSum - rolledSum.size))
        rolledSum = np.append(rolledSum,rolledSumMissing)
        rolledSum = rolledSum.reshape((bins_[0].size - 1, bins_[1].size - 1))
        return rolledSum , bins_[0], bins_[1]
        
        
        
        
        

        
        
    def LoadElectrons(self,elFileName):
         elFilePath = f"{self.Path}" + f'{elFileName}' + '*sdf'
         self.elList = glob.glob(f'{elFilePath}')
         self.elList = sorted(self.elList)
         #get photon matrix
         lastFile  = sdfhelper.getdata(f'{self.elList[-1:]}')
             # try:
                 
             # except:
             #     break
         