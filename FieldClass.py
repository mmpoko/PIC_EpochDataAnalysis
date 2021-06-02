#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
import sdf_helper as sdfhelper
import os
import glob

class FieldClass:
    def __init__(self, *args,**kwargs):
        self.t = 0
        for dictionary in args:
            for key in dictionary:
                setattr(self,key,dictionary[key])
        for key in kwargs:
            setattr(self,key,kwargs[key])
    
    def LoadFieldAtTime(self,FieldFileName,fileNum, *args):
        FieldVar = args[0]
        FieldFilePath = f"{self.Path}" + f'{FieldFileName}' + f'{fileNum:04d}' + '.sdf'
        curFile = sdfhelper.getdata(fr'{FieldFilePath}')
        print(f"time = {curFile.Header['time'] * 1e15} fs")
        #load in electric field quantities
        for field in dir(curFile):
            if 'Electric_Field_Ex' in field and FieldVar in 'Electric_Field_Ex':
                self.ex = eval(f'curFile.{field}.data')
                self.x = eval(f'curFile.{field}.grid.data[0][1:]')
                self.y = eval(f'curFile.{field}.grid.data[1][0:-1]')   
                break
            if 'Electric_Field_Ey' in field and FieldVar in 'Electric_Field_Ey':
                self.ey = eval(f'curFile.{field}.data')
                self.x = eval(f'curFile.{field}.grid.data[0][0:-1]')
                self.y = eval(f'curFile.{field}.grid.data[1][1:]')   
                break
            if 'Electric_Field_Ez' in field and FieldVar in 'Electric_Field_Ez':
                self.ez = eval(f'curFile.{field}.data')
                self.x = eval(f'curFile.{field}.grid.data[0][0:-1]')
                self.y = eval(f'curFile.{field}.grid.data[1][0:-1]')   
                break
            if 'Magnetic_Field_Bx' in field and FieldVar in 'Magnetic_Field_Bx':
                self.bx = eval(f'curFile.{field}.data')
                self.x = eval(f'curFile.{field}.grid.data[0][0:-1]')
                self.y = eval(f'curFile.{field}.grid.data[1][1:]')   
                break
            if 'Magnetic_Field_By' in field and FieldVar in 'Magnetic_Field_By':
                self.by = eval(f'curFile.{field}.data')
                self.x = eval(f'curFile.{field}.grid.data[0][1:]')
                self.y = eval(f'curFile.{field}.grid.data[1][0:-1]')   
                break
            if 'Magnetic_Field_Bz' in field and FieldVar in 'Magnetic_Field_Bz':
                self.bz = eval(f'curFile.{field}.data')
                self.x = eval(f'curFile.{field}.grid.data[0][1:]')
                self.y = eval(f'curFile.{field}.grid.data[1][1:]')   
                break    
        self.t = curFile.Header['time']
        return curFile.Header['time']



        
    