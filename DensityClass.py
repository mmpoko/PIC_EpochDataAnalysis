#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 14 15:02:11 2021

@author: michaelpokornik
"""
import numpy as np
import sdf_helper as sdfhelper
import os
import glob

class DensityClass:
    def __init__(self, *args,**kwargs):
        self.t = 0
        self.ndim = 3
        for dictionary in args:
            for key in dictionary:
                setattr(self,key,dictionary[key])
        for key in kwargs:
            setattr(self,key,kwargs[key])
    
    def LoadDensityAtTime(self,DensFileName,fileNum,species,normF):
        DensFilePath = f"{self.Path}" + f'{DensFileName}' + f'{fileNum:04d}' + '.sdf'
        curFile = sdfhelper.getdata(fr'{DensFilePath}')
        print(f"time = {curFile.Header['time'] * 1e15} fs")
        #load in electric field quantities
        for field in dir(curFile):
            if f'Derived_Number_Density_{species}' in field:
                self.dens = eval(f'curFile.{field}.data')
                self.dens = self.dens / normF
                self.x = eval(f'curFile.{field}.grid_mid.data[0]')
                self.y = eval(f'curFile.{field}.grid_mid.data[1]')   
                break
   
    def LoadEkbarAtTime(self,EkFileName,fileNum,species,normF):
        EkFilePath = f"{self.Path}" + f'{EkFileName}' + f'{fileNum:04d}' + '.sdf'
        curFile = sdfhelper.getdata(fr'{EkFilePath}')
        print(f"time = {curFile.Header['time'] * 1e15} fs")
        #load in 
        self.x = eval(f'curFile.{species}.grid_mid.data[0]')
        self.y = eval(f'curFile.{species}.grid_mid.data[1]')
        if self.ndim == 3:
            self.z = eval(f'curFile.{species}.grid_mid.data[2]')
        self.Ek = eval(f'curFile.{species}.data') 
        self.Ek = self.Ek / normF
        self.t = curFile.Header['time']
        return curFile.Header['time']