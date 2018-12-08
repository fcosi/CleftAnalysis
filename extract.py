#!/usr/bin/env python3

# class for the extraction of values from simulation output to pandas DataFrame

import pandas as pd
import os as os
import re as re
import glob as glob
import sys as sys
import numpy as np

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
        self.crunumber =  int(self.param["z_discs"][0]*self.param["crus_per_sec_line"][0]**2 - 4)

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

    def crusInfo(self, getRadius = False, getLocations = False):
        '''get total number of RyR and LCCs. Optional radius and channel locations
        
        output are lists of the two channel types optional if
        variables set to True also output CRU radius and channel
        locations (in lists)
        
        '''
        ryr_numbers = []
        lcc_numbers = []
        radii = []
        ryr_location = []
        lcc_location = []
        for i in range(self.crunumber):
            crufile = open(self.folder + "clefts/cleft" + str(i) + ".log")
            cruinfo = crufile.readline()
            crufile.close()
            cruinfo = cruinfo.strip().replace("  ", " ").split(" ")
            # append number of RyR, LCCs and Radius to list
            ryr_numbers.append(int(cruinfo[cruinfo.index("RyRs") - 1]))
            lcc_numbers.append(int(cruinfo[cruinfo.index("LCCs") - 1]))
            radii.append(float(cruinfo[cruinfo.index("Radius") - 1]))
            # get the locations of the channels only if wished
            if getLocations:
                tmp1, tmp2 = self.__getLocations(ryr_numbers[i], lcc_numbers[i], cruinfo)
                ryr_location.append(tmp1)
                lcc_location.append(tmp2)
        
        # return different lists according to the wished arguments
        if getRadius and not getLocations:
            return ryr_numbers, lcc_numbers, radii
        elif getLocations and not getRadius:
            return ryr_numbers, lcc_numbers, ryr_location, lcc_location
        elif getRadius and getLocations:
            return ryr_numbers, lcc_numbers, radii, ryr_location, lcc_location
        else:
            return ryr_numbers, lcc_numbers

    def __getLocations(self, nryr, nlcc, cru_info_list):
        '''
        To be used only for crusInfo()
        
        gets x y pairs from first line of the cleft logs
        
        returns np 2d array
        '''
        loc_index = cru_info_list.index("Location")
        cru_info_list = cru_info_list[loc_index + 1:]
        ryr_loc = np.zeros([2, nryr]) 
        lcc_loc = np.zeros([2, nlcc])
        
        if ('(' in cru_info_list[0] or '(' in cru_info_list[1]):    
            for ind, pair in enumerate(cru_info_list[:nryr]):
                loc = pair[1:-1].split(";")
                ryr_loc[0,ind] = float(loc[0])
                ryr_loc[1,ind] = float(loc[1])
            for ind, pair in enumerate(cru_info_list[nryr: nryr + nlcc]):
                loc = pair[1:-1].split(";")
                lcc_loc[0,ind] = float(loc[0])
                lcc_loc[1,ind] = float(loc[1])
        else:
            import warnings
            warnings.warn("""\nExtracting the Locations 'old' saving style, be AWARE they probably will cause numerical issues, since the output is not consistent.""")
            
            for ind, num in enumerate(cru_info_list[:nryr-1:2]):
                ryr_loc[0,ind] = float(cru_info_list[ind])
                ryr_loc[1,ind] = float(cru_info_list[ind + 1])
            
            for ind, num in enumerate(cru_info_list[nryr + 1: nryr + nlcc:2]):
                lcc_loc[0,ind] = float(cru_info_list[ind])
                lcc_loc[1,ind] = float(cru_info_list[ind + 1])
        return ryr_loc, lcc_loc

    def saveCleftPlot(self, mergeOutput = False):
        """
        Save a plot of each cleft using the channel Locations
        """
        import os
        outputname = self.folder + "clefts/cleftall.pdf"
        if (os.path.exists(outputname)):
            sys.exit("cleft ouput plots already exist!")
        
        import matplotlib.pyplot as plt
        nR, nL, radii, rloc, lloc = self.crusInfo(getRadius=True, getLocations=True)
        savepath = self.folder + "clefts/cleftplot{}.pdf"
        for i in range(self.crunumber):
            rad = radii[i]
            fig, ax = plt.subplots()
            ax.plot(rloc[i][0,], rloc[i][1,], "go", label="{} RyR".format(nR[i]))
            ax.plot(lloc[i][0,], lloc[i][1,], "bo", label="{} LCC".format(nL[i]))
            circ = plt.Circle((0,0), rad, color='k', fill=False)
            ax.add_artist(circ)
            ax.set_xlabel("x [$\mu m$]")
            ax.set_ylabel("y [$\mu m$]")
            ax.set_xlim((-1.5*rad, 1.5*rad))
            ax.set_ylim((-1.5*rad, 1.5*rad))
            ax.legend()
            fig.savefig(savepath.format(i), bbox_inches = "tight")
        if mergeOutput:            
            savepath = savepath.format("*")
            print(savepath)
            os.system('pdfunite {} {}clefts/cleftall.pdf'.format(savepath, self.folder))


    def getOpenChannels(self):
        '''
        get time, open RyR and open LCC from file in clefts/
        if not existing create it
        
        - Return:
        pandas DataFrame
        '''
        
        outputname = self.folder + "clefts/openChannels.csv"
        if (not os.path.exists(outputname)):
            self.saveOpenChannels()
        
        df = pd.read_csv(outputname, sep=" ")
        df = df.rename(index=str, columns={"Unnamed: 0": "counter0", "Unnamed: 1": "counter1"})
        return df

    def saveOpenChannels(self):
        '''
        save time steps, number of open RyR and LCC to file in clefts/ dir
        
        ATTENTION: works properly only with not merged cleftdyn branch (addition_CRU_debug_output)
        (2018/10/13)
        where cleft.log output has been changed (num of open channels outputted)
        '''

        outputname = self.folder + "clefts/openChannels.csv"
        if (os.path.exists(outputname)):
            sys.exit("cleft ouput for open channels already exists!")

        totDF = pd.DataFrame()
        for i in range(self.crunumber):
            timesteps = []
            openRyR = []
            openLCC = []
            cleftname = "cleft" + str(i)
            # get all lines from cleftlogs
            lines, totalRyR, totalLCC = self.__getCleftLogLines(i)
            # dirty hack: list starts from 1 since the first time step with time=0 (no dot for float) is not covered by regular expression (CHANGE THIS!)
            for line in lines[1:]:
                nums = re.findall("\d+\.\d+", line)
                timesteps.append(nums[0])
                nums = re.findall("\d+", line)
                openRyR.append(int(nums[-2]))
                openLCC.append(int(nums[-1]))
            
            # create pandas dataframe to save file to csv
            outDF = pd.DataFrame({'clefts': cleftname,'time': timesteps, 'openRyR': openRyR,
                                  'openLCC': openLCC, 'totalRyR': totalRyR, 'totalLCC': totalLCC})
            #dict[cleftname] = outDF
            del(timesteps[:], openRyR[:], openRyR[:])
            totDF = pd.concat([totDF, outDF])
        totDF["ratioRyR"] = totDF["openRyR"]/totDF["totalRyR"]
        totDF["ratioLCC"] = totDF["openLCC"]/totDF["totalLCC"]
        outfile = open(outputname, 'a')
        totDF.to_csv(outfile, header=True, sep=" ")
        outfile.close()


    def getChannelConcentrations(self, cleftnr=0):
        """
        get Ca concentrations at channel mouth from cleft logs given the cleftnumber
        
        Args:
        - cleftnumber
        
        Return:
        - numpy arrays with channel concentrations and timesteps
        """
        lines, totR, totL = self.__getCleftLogLines(cleftnr)
        times = np.zeros([1, len(lines)])
        ryr_conc = np.zeros([len(lines), totR])
        lcc_conc = np.zeros([len(lines), totL])
        first_channel_ind = 5
        for ind, line in enumerate(lines):
            line_list = line.replace("  ", " ").split(" ")
            times[:,ind] = line_list[0]
            line_list = line_list[first_channel_ind:]
            ryr_conc[ind,:] = line_list[1:3*totR:3]
            line_list = line_list[3*totR:]
            lcc_conc[ind,:] = line_list[1:3*totL:3]
        
        return times, ryr_conc, lcc_conc

    def __getCleftLogLines(self, crunum):
        """
        private fct returning the lines of the cleftlogs and the total number of RyR LCCs
        without the first line
        used by getChannelConcentration() and saveOpenChannels() fcts.
        """
        crufile = open(self.folder + "clefts/cleft" + str(crunum) + ".log")
        lines = [line.rstrip('\n') for line in crufile]
        crufile.close()
        totalRyR = int(lines[0].split(" ")[0])
        totalLCC = int(lines[0].split(" ")[2])        
        del(lines[0])
        return lines, totalRyR, totalLCC
    
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
        # sentence to be matched
        match = " number of open crus: "
        with open(out_file, "r") as fp:
            time_lines = [line for line in fp if match in line]
        for line in time_lines:
            nums = re.findall("\d+", line)
            # the first two ints are the time
            time.append(float(nums[0] + "." + nums[1]))
            # the third int found is the num of open crus
            crus.append(int(nums[2]))
        return time, crus
    
    def get_linescan_cru_info(self, start_time = -np.infty, end_time = np.infty, min_z = -np.infty, max_z = np.infty):
        '''
            extracts cru_info_rank*.csv files from linescan folder
        '''
        # find files in subfolder
        
        linescan_folder = self.folder + "linescan/"
        path, dirs, files = next(os.walk(linescan_folder))
        
        frames_cru_info = []
        
        for filename in files:
            
            if filename.find("cru_info") == -1:
                continue
            
            temp_cru_info = pd.read_csv(linescan_folder + filename)
    
            #print temp_cru_info
            #
            '''
            z-disc selection
            '''
            
            temp_cru_info = temp_cru_info[temp_cru_info['cru_z'] >= min_z]
            temp_cru_info = temp_cru_info[temp_cru_info['cru_z'] <= max_z]
            
            '''
            time selection
            '''
            
            temp_cru_info = temp_cru_info[temp_cru_info['time'] >= start_time]
            
        
            '''
            sort values
            '''

            temp_cru_info = temp_cru_info.sort_values(by=["time", "cru_id"])
            
       
            '''
            round values to needed precission
            '''
            
            temp_cru_info = temp_cru_info.round({'time': 1, 'cru_flux': 2, 'cyto_ca2+': 2, 'cyto_bm': 2, 'cyto_bs':2, 'sarco_ca': 0 , 'fluo4': 2})
            
            '''
            delete round postion values and delete duplicates to reduce memory cost
            '''
            temp_cru_info = temp_cru_info.drop_duplicates()

            # uncomment the following to see memory usage

            '''
            append to data frame
            '''
            
            frames_cru_info.append(temp_cru_info)

        cru_info_all = pd.concat(frames_cru_info)

        cru_info_all = cru_info_all.sort_values(by=["time", "cru_id"])
        
        return cru_info_all
