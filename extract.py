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
    
    import with:
    from CleftAnalysis import extract
    
    ------------------------------
    to initiate the class call e.g.
    f = extract.extract(folder_num, subfolder_name, values)
    where:
    folder_num : number of the generated cleftdyn folder
    subfolder_name : string of the possible subfolder. If nothing specified, then no subfolder set
                     numerical values possible: 0 and 1 are special cases defined in determineFolder
    values : name or name snippet of the values to be extracted
   
    """
    
    def __init__(self, folder_num, subfolder = False, values = "ionic"):
        self.folder = folder_num
        # setFolder var ensures a folder is not initiated twice
        self.setFolder = False        
        
        self.subfolder = subfolder
        self.determineFolder()
        
        self.values = values
        self.determineValues()
        
        self.folderDefinition()
        self.printInfo()
        self.param = self.getParam()

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
              "\nis it a continuation from long simulations? If yes, which? ", self.subfolder)

    def folderDefinition(self, new_folder = ""):
        if (not self.setFolder and not self.subfolder):
            self.folder = "../" + str(self.folder) + "/"
            self.setFolder = True
        elif (not self.setFolder and (self.subfolder == True)):
            self.folder = "../long-cont/" + str(self.folder) + "/"
            self.setFolder = True
        elif (not self.setFolder and (type(self.subfolder) == str)):
            self.folder = "../" + self.subfolder + "/" + str(self.folder) + "/"
            self.setFolder = True
        else:
            print("folder already set or not valid")

    def determineFolder(self):
        if (self.subfolder == 0):
            self.subfolder = False
        elif (self.subfolder == 1):
            self.subfolder = True
        else:
            self.subfolder = self.subfolder

    def getParam(self):
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

    def getParameters(self):
        return self.param

    def crusInfo(self):
        crunumber = int(self.param["z_discs"][0]*self.param["crus_per_sec_line"][0]**2 - 4)

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
