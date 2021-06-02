#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
import sdf_helper as sdfhelper
import os
import glob

class ProbeDataClass:
    def __init__(self, *args,**kwargs):
        self.MeV = (1.6022e-19)*(1.0e6) 
        self.c_speed = 299792458 # m/s
        self.mc = 9.11e-31*self.c_speed #kg m /s
        for dictionary in args:
            for key in dictionary:
                setattr(self,key,dictionary[key])
        for key in kwargs:
            setattr(self,key,kwargs[key])
    
    def LoadProbeAtTime(self,probFileName,fileNum, *restrict):
        probeFilePath = f"{self.Path}" + f'{probFileName}' + f'{fileNum:04d}' + '.sdf'
        curFile = sdfhelper.getdata(fr'{probeFilePath}')
        print(f"time = {curFile.Header['time'] * 1e15} fs")
        #load in electron ids first
        spatial = 'Grid_Probe_' + self.probeName

        self.x = eval(f'curFile.{spatial}.data[0]')
        self.y = eval(f'curFile.{spatial}.data[1]')
        self.z = eval(f'curFile.{spatial}.data[2]')  
        
        pxName = self.probeName + '_Px'
        pyName = self.probeName + '_Py'
        pzName = self.probeName + '_Pz'
        self.px = eval(f'curFile.{pxName}.data')
        self.py = eval(f'curFile.{pyName}.data')  
        self.pz = eval(f'curFile.{pzName}.data')
        self.p2 = self.px * self.px + self.py * self.py + self.pz * self.pz
        
        self.phi = np.arctan2(self.py,self.pz) * 180 / np.pi
        self.theta = np.arctan2(self.px,self.pz) * 180 / np.pi
        self.gam = np.sqrt(1 + (self.p2 / self.mc / self.mc) )
        self.E0 = self.gam * self.mc * self.c_speed
        weight = self.probeName + '_weight'
        self.w = eval(f'curFile.{weight}.data')
        self.eleInd = np.ones(self.w.shape,dtype=int)
        
        for dictionary in restrict:
            for key in dictionary:
                if 'XLim' in key:
                    self.XLimRestrict(dictionary[key][0:-1],dictionary[key][-1])
                elif 'YLim' in key:
                    self.YLimRestrict(dictionary[key][0:-1],dictionary[key][-1])
                elif 'ZLim' in key:
                    self.ZLimRestrict(dictionary[key][0:-1],dictionary[key][-1])
                elif 'PxLim' in key:
                    self.PxLimRestrict(dictionary[key][0:-1],dictionary[key][-1])
                elif 'PyLim' in key:
                    self.PyLimRestrict(dictionary[key][0:-1],dictionary[key][-1])  
                elif 'PzLim' in key:
                    self.PzLimRestrict(dictionary[key][0:-1],dictionary[key][-1])
                elif 'GammaLim' in key:
                    self.GammaLimRestrict(dictionary[key][0:-1],dictionary[key][-1])
                elif 'E0Lim' in key:
                    self.E0LimRestrict(dictionary[key][0:-1],dictionary[key][-1])
                elif 'IDLim' in key:
                    self.IDRestrict(dictionary[key])
                elif 'eIDLim' in key:
                    self.IDexactRestrict(dictionary[key][0],dictionary[key][-1])
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
    def ZLimRestrict(self,bounds,units):
        self.eleInd[(self.z * units < bounds[0]) & (np.isnan(bounds[0]) == False)]= 0
        self.eleInd[(self.z * units > bounds[1]) & (np.isnan(bounds[1]) == False)]= 0    
    def PxLimRestrict(self,bounds,units):
        self.eleInd[(self.px * units < bounds[0]) & (np.isnan(bounds[0]) == False)]= 0
        self.eleInd[(self.px * units > bounds[1]) & (np.isnan(bounds[1]) == False)]= 0
    def PyLimRestrict(self,bounds,units):
        self.eleInd[(self.py * units < bounds[0]) & (np.isnan(bounds[0]) == False)]= 0
        self.eleInd[(self.py * units > bounds[1]) & (np.isnan(bounds[1]) == False)]= 0
    def PzLimRestrict(self,bounds,units):
        self.eleInd[(self.pz * units < bounds[0]) & (np.isnan(bounds[0]) == False)]= 0
        self.eleInd[(self.pz * units > bounds[1]) & (np.isnan(bounds[1]) == False)]= 0    
    def GammaLimRestrict(self,bounds,units):
        self.eleInd[(self.gam * units < bounds[0]) & (np.isnan(bounds[0]) == False)]= 0
        self.eleInd[(self.gam * units > bounds[1]) & (np.isnan(bounds[1]) == False)]= 0     
    def E0LimRestrict(self,bounds,units):
        self.eleInd[(self.E0 * units < bounds[0]) & (np.isnan(bounds[0]) == False)]= 0
        self.eleInd[(self.E0 * units > bounds[1]) & (np.isnan(bounds[1]) == False)]= 0
    def IDRestrict(self,ids):
        self.eleInd[ids == 0] = 0
    def IDExactRestrict(self,ids):
        self.eleInd[~np.isin(self.eIds,ids)] = 0
    def GetPositions(self):
        return self.x[self.eleInd > 0] , self.y[self.eleInd > 0]
    
        
    

    
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
        
        
        
        
        

        
        
         