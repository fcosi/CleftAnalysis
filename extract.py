#!/usr/bin/env python3

# class for the extraction of values from simulation output to pandas DataFrame

import pandas as pd
import os as os

class extract:
    """
    In the constructor give folder from which to extract 
    other arguments are:
    - long cont folder (standard value False)
    - which file to extract (standard is ionicmodel.txt)    
    """
    
    def __init__(self, folder_num, long_cont = False, values = "ionic"):
        self.folder = folder_num
        self.setFolder = False        
        
        self.long_cont = long_cont
        self.determineIfLong()
        
        self.values = values
        self.determineValues()
        
        self.folderDefinition()
        self.printInfo()

    def read(self):
        var = pd.read_csv(self.folder + self.values, delim_whitespace = True)
        return var

    def readMass(self):
        self.values = "massconversionlog.txt"
        return self.read()

    def readIonic(self):
        self.values = "ionicmodel.txt"
        return self.read()

    def printInfo(self):
        print("extraction folder: ", self.folder,
              "\nis it a continuation from long simulations? ", self.long_cont)

    def folderDefinition(self):
        if (not self.setFolder and not self.long_cont):
            self.folder = "../" + str(self.folder) + "/"
            self.setFolder = True
        elif (not self.setFolder and self.long_cont):
            self.folder = "../long-cont/" + str(self.folder) + "/"
            self.setFolder = True
        else:
            print("folder already set or not valid")
    
    def getParameters(self):
        final_param = pd.DataFrame()
        param = pd.read_csv(self.folder + "parameters.txt", sep="=", index_col=0).T
        for val in param:
            try:
                final_param[val] = [float(param[val][0])]
            except ValueError:
                #print(val, " has no numerical value")
                pass
            except KeyError:
                #print(val, " raised an key error")
                pass
        return final_param

    def determineIfLong(self):
        if (self.long_cont == 0):
            self.long_cont = False
        elif (self.long_cont == 1):
            self.long_cont = True
        else:
            self.long_cont = self.long_cont        

    def determineValues(self):
        if ((self.values == "ionic") or (self.values == "ionicmodel") or
            (self.values == "ionicmodel.txt") or (self.values == 1) or (self.values == "ion")):
            self.values = "ionicmodel.txt"
        elif ((self.values == "mass") or (self.values == "massconversion") or
              (self.values == "massconversionlog.txt") or (self.values == 0) or
              (self.values == "massconversionlog")):
            self.values = "massconversionlog.txt"
        else:
            print("no considered/implemented value, skip")        
