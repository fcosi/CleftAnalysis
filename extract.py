#!/usr/bin/env python3

# class for the extraction of values from simulation output to pandas DataFrame

import pandas as pd
import os as os
import re as re
import glob as glob

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
        '''
        print info when initiating class
        '''
        print("extraction folder: ", self.folder,
              "\nis it a continuation from long simulations? If yes, which? ", self.subfolder)

    def folderDefinition(self, new_folder = ""):
        '''
        set folder from which to extract data
        Args:
        new_folder: if set can allow changing the actual folder
        '''
        if (not self.setFolder and not self.subfolder):
            self.folder = "../" + str(self.folder) + "/"
            self.setFolder = True
        elif (not self.setFolder and (self.subfolder == True)):
            self.folder = "../long-cont/" + str(self.folder) + "/"
            self.setFolder = True
        elif (not self.setFolder and (type(self.subfolder) == str)):
            self.folder = "../" + self.subfolder + "/" + str(self.folder) + "/"
            self.setFolder = True
        elif (self.setFolder and new_folder):
            print("folder changed to: " + new_folder)
            self.folder = new_folder
        else:
            print("folder already set to: " + self.folder + " or not valid")

    def determineFolder(self):
        ''''
        look for subfolder arg and set it if found
        '''
        if (self.subfolder == 0):
            self.subfolder = False
        elif (self.subfolder == 1):
            self.subfolder = True
        else:
            self.subfolder = self.subfolder

    def getParam(self):
        '''
        get dataframe of parameters from parameter file of cleftdyn (use var names as keys)
        '''
        final_param = pd.DataFrame()
        param = pd.read_csv(self.folder + "parameters.txt", sep="=", index_col=0).T
        for val in param:
            try:
                final_param[val] = [float(param[val][0])]
            except ValueError:
                pass
            except KeyError:
                pass
        return final_param       

    def getParameters(self):
        return self.param

    def crusInfo(self):
        '''
        get total number of RyR and LCCs
        output are lists of the two channel types
        '''
        crunumber = int(self.param["z_discs"][0]*self.param["crus_per_sec_line"][0]**2 - 4)
        ryr_numbers = []
        lcc_numbers = []
        for i in range(crunumber):
            crufile = open(self.folder + "clefts/cleft" + str(i) + ".log")
            cruinfo = crufile.readline()
            crufile.close()
            cruinfo = cruinfo.split(" ")
            # append number of RyR and LCCs to list
            ryr_numbers.append(int(cruinfo[0]))
            lcc_numbers.append(int(cruinfo[2]))
        return ryr_numbers, lcc_numbers

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

    def get_openCRUs(self, year=2018):
        """
        Gets the number of open crus from the outputfile (year is the starting name of the file)
        also gets the according times
        
        Output:
        - list of time steps
        - list of open crus
        """
        out_file = self.folder + str(year) + "*"
        out_file = glob.glob(out_file)[0]
        time = []
        crus = []
        match = " number of open crus: "
        with open(out_file, "r") as fp:
            time_lines = [line for line in fp if match in line]
        for line in time_lines:
            nums = re.findall("\d+", line)
            time.append(float(nums[0] + "." + nums[1]))
            crus.append(int(nums[2]))
        return time, crus
