"""
Class Script with functions for the Analysis (later may be also for the evaluation)

import with:
from CleftAnalysis import functions
"""

# Imports
import numpy as np
import glob
import os
import re as re
import pandas as pd
import scipy as sp
import scipy.stats
from scipy.signal import argrelextrema
from sklearn.linear_model import *
from sklearn.model_selection import cross_val_score
import chaospy as cp

import matplotlib.pyplot as plt

class SparkAnalysis:
    """
    store analysis functions for the spark analysis
    """
    
    def __init__(self):
        '''
            initialization for the spark cascade analysis
        '''
        self.datafr = pd.DataFrame()
        self.list_of_sparks = []
        self.cascade = []
        #self.folder = folder
        #self.processes = self.get_used_processes()
        #self.cru_info = self.get_cru_info()

## ----------------------------------------------------------------------------------------------

    #def get_cru_info(self):
        #'''
        #Function collecting information on the CRUs saved in the linescan
        #'''
        #frames_cru_info = []
        #for i in range(0, self.processes):        
            #filename = self.folder + 'linescan/cru_info_rank' + str(i) + ".csv"
            #try:
                #temp_cru_info = pd.read_csv(filename)
            #except FileNotFoundError:
                #print("file not found, probably no linescan activated!")
                #break
            ##    sort values
            #temp_cru_info = temp_cru_info.sort_values(by=["time", "cru_id"])
            ## round values to needed precission
            ##temp_cru_info = temp_cru_info.round({'time': 1, 'cru_flux': 2, 'cyto_ca2+': 4, 'cyto_bm': 2, 'cyto_bs':2, 'sarco_ca': 0 , 'fluo4': 2})        
            ## delete round postion values and delete duplicates to reduce memory cost       
            #temp_cru_info = temp_cru_info.drop_duplicates()
            ##    append to data frame
            #frames_cru_info.append(temp_cru_info)
        #try:
            ##print(frames_cru_info)
            #frames_info = pd.concat(frames_cru_info)
            ##print(frames_cru_info)
            #frames_info = frames_info.sort_values(by=["time", "cru_id"])
        #except ValueError or AttributeError:
            #frames_info = []
            #pass
        #return frames_info

## ----------------------------------------------------------------------------------------------

    #def get_linescan_info(self, spark_no = 1, duration = 1):
        #'''
        #return time, cru_id, scan_direction and scan_direction_pos of a linescan
        
        #Args
        #----------
        #spark_no: spark id
        #duration: time interval in ms from which to take the spark data
        #'''
        ## check spark number is lower than len(spark_candidates_startTimes)
        #frames_linescan = []
        #spark_candidates_startTimes, spark_candidates_endTimes = self.getSparkStartEndTimes(duration)
        #for i in range(0, self.processes):        
            #filename = self.folder + "linescan/" + 'linescan_rank' + str(i) + ".csv"
            ## get linescans and adjust them to spark
            #temp_linescan = pd.read_csv(filename)
            #temp_linescan = temp_linescan[temp_linescan['time'] >= spark_candidates_startTimes[spark_no-1]]
            #temp_linescan = temp_linescan[temp_linescan['time'] < spark_candidates_endTimes[spark_no-1]]
            ## z-disc selection         
            ##temp_cru_info = temp_cru_info[temp_cru_info['cru_z'] >= 3.0]
            ##temp_cru_info = temp_cru_info[temp_cru_info['cru_z'] <= 5.0]    
            
            #temp_linescan = temp_linescan.sort_values(by=["time", "cru_id" , "scan_direction", "scan_direction_pos"])
            ## round values to needed precission
            ## temp_linescan = temp_linescan.round({'time': 1, 'scan_direction_pos':1, 'cyto_ca2+': 2, 'cyto_bm': 2, 'cyto_bs':2, 'sarco_ca': 0 , 'fluo4': 2})
            #temp_linescan = temp_linescan.drop_duplicates()
            ## append to tmp list
            #frames_linescan.append(temp_linescan)            
        #linescan_all = pd.concat(frames_linescan)        
        #linescan_all = linescan_all.sort_values(by=["time", "cru_id" , "scan_direction", "scan_direction_pos"])
        #return linescan_all
## ----------------------------------------------------------------------------------------------

    #def get_used_processes(self):
        #'''
        #Function that reads the output file and returns the number of used processes
        #'''
        #with open(glob.glob(self.folder + "20*")[0]) as f:
            #first = f.readline()
            #f.close()
        #return (int(re.search(r'\d+', first).group()))

# ----------------------------------------------------------------------------------------------

    def plotLineScan(self, species = "cyto_ca2+", spark_no = 1, duration = 1):
        """
        plot the linescan given a species, a spark number and the spark duration
        DISCLAIMER: depending on file size could take long to plot!
        """
        linescan_all = self.get_linescan_info(spark_no, duration)
        plt.figure(figsize=(len(linescan_all['scan_direction'].unique())*6, len(linescan_all['cru_id'].unique())*5))
        
        plot_nr = 0
        
        for cru_id in linescan_all['cru_id'].unique():
            single_cru_df = linescan_all[linescan_all['cru_id'] == cru_id]
            max_val = max(single_cru_df[species])
            min_val = min(single_cru_df[species])
            
            temp = single_cru_df[single_cru_df[species] == max_val]
            time_peak = max(temp['time'])
            
            
            for scan_direction in single_cru_df['scan_direction'].unique():
                plot_nr += 1
                temp_dir = single_cru_df[single_cru_df['time'] == time_peak]
                temp_dir = temp_dir[temp_dir['scan_direction'] == scan_direction]
                positions = temp_dir['scan_direction_pos'].unique()
                species_val = []
                for position in positions:
                    temp_val = temp_dir[temp_dir['scan_direction_pos'] == position ]
                    species_val.append(temp_val[species].max())
        
                plt.subplot( len(linescan_all['cru_id'].unique()), len(single_cru_df['scan_direction'].unique()), plot_nr)
                plt.title("CRU ID: %d"  % cru_id )
                plt.ylabel(species)
                plt.xlabel(str(scan_direction) + "-pos [$\mu$m]")
                plt.ylim(min_val*0.9,max_val*1.1)
                plt.tight_layout()
                plt.plot(positions,species_val)   
        plt.show()

# ----------------------------------------------------------------------------------------------

    def printFW_FWHM(self, species = 'cyto_ca2+', spark_no = 1, duration = 1):
        """
        Function to print the FW and FWHM given: species, spark number and duration
        """
        linescan_all = self.get_linescan_info(spark_no, duration)
        for cru_id in linescan_all['cru_id'].unique():
            single_cru_df = linescan_all[linescan_all['cru_id'] == cru_id]
            max_fluo4 = max(single_cru_df[species])
            base_fluo4 = min(single_cru_df[species])
            
            delta_F_over_F0 = (max_fluo4-base_fluo4)/base_fluo4
            print("(Delta F)/F0: %f" % (delta_F_over_F0))
            
            temp = single_cru_df[single_cru_df[species] == max_fluo4]
            time_peak = max(temp['time'])
            
            
            for scan_direction in single_cru_df['scan_direction'].unique():
                FW = 0.0
                FWHM = 0.0
                
                temp_dir = single_cru_df[single_cru_df['time'] == time_peak]
                temp_dir = temp_dir[temp_dir['scan_direction'] == scan_direction]
                positions = temp_dir['scan_direction_pos'].unique()
                species_val = []
                
                start_pos = positions[0] 
                last_pos = start_pos
                
                for position in positions:
                    temp_val = temp_dir[temp_dir['scan_direction_pos'] == position ]
                    fluo4_val = temp_val[species].mean()
                    
                    #print fluo4_val
                    if ((fluo4_val - base_fluo4) > (max_fluo4 - base_fluo4)/10.0):
                        FW += position - last_pos
                        
                    if ((fluo4_val - base_fluo4) > (max_fluo4 - base_fluo4)/2.0):
                        FWHM += position - last_pos
                        
                    last_pos = position 

                print("CRU ID: %d" % cru_id)
                print("direction: " + str(scan_direction))
                print("FW: %f um" % FW)
                print("FWHM: %f um" % FWHM)
                print("#"*23)  

# ----------------------------------------------------------------------------------------------

    # THIS FUNCTION CAN BE IMPROVED AND IS PROBABLY WRONG!
    def getSparkStartEndTimes(self, duration = 1):
        """
        Function to get the list of the ending time point of sparks
        Input:
        - spark duration in ms
        """
        time_points = self.cru_info['time'].unique()
        spark_candidates_startTimes = [time_points[0]]
        spark_candidates_endTimes = []
        last_time = time_points[0]
        
        # look fo all time lapses longer than 1ms
        for time in time_points:
            if(abs(time - last_time) > duration):
                spark_candidates_endTimes.append(last_time)
                spark_candidates_startTimes.append(time)
            last_time = time
        spark_candidates_endTimes.append(time) 
        return spark_candidates_startTimes, spark_candidates_endTimes

    def mem_usage(self, pandas_obj):
        '''
        return mem_usage given a pandas object
        '''
        if isinstance(pandas_obj,pd.DataFrame):
            usage_b = pandas_obj.memory_usage(deep=True).sum()
        else: # we assume if not a df it's a series
            usage_b = pandas_obj.memory_usage(deep=True)
            usage_mb = usage_b / 1024 ** 2 # convert bytes to megabytes
        return usage_mb #"{:03.2f} MB".format(usage_mb)

    def printFD_FDHW(self, duration = 1):
        spark_no = 0
        spark_candidates_startTimes, spark_candidates_endTimes = self.getSparkStartEndTimes(duration)
        for start_time, end_time in zip(spark_candidates_startTimes, spark_candidates_endTimes):
            spark_no += 1
            print("Spark no: %d" % spark_no)
            spark_df = self.cru_info[self.cru_info['time'] >= start_time]
            spark_df = spark_df[spark_df['time'] <= end_time]
            max_fluo4 = max(spark_df['fluo4']) 
            base_fluo4 = min(spark_df['fluo4'])
            delta_F_over_F0 = (max_fluo4-base_fluo4)/base_fluo4
            print("(Delta F)/F0: %f" % (delta_F_over_F0))
            
            FD = 0.0
            FDHM = 0.0
            
            old_time = start_time
            for time in spark_df['time'].unique():
                
                temp = spark_df[spark_df['time']==time]
                fluo4_val = max(temp['fluo4'])
                
                if ((fluo4_val - base_fluo4) > (max_fluo4 - base_fluo4)/10.0):
                    FD += time - old_time
                    
                if ((fluo4_val - base_fluo4) > (max_fluo4 - base_fluo4)/2.0):
                    FDHM += time - old_time                
                old_time = time
            
            print("FD: %f ms" % FD)
            print("FDHM: %f ms" % FDHM)
            print("#"*23)    
            
    def get_cru_positions(self,cru_info):
        
        crus_x = []
        crus_y = []
        crus_z = []
        crus_start_time = []

        for cru_id in cru_info['cru_id'].unique():
            cru_df = cru_info[cru_info['cru_id'] == cru_id]
            temp = cru_df.iloc[0]
    
            crus_x.append(temp['cru_x'])
            crus_y.append(temp['cru_y'])
            crus_z.append(temp['cru_z'])
            crus_start_time.append(temp['time'])
            
        pos = np.array(zip(crus_x,crus_y,crus_z,crus_start_time))
        
        return pos
    
    def get_cru_cluster(self, positions, scale_time=0.1, threshold=3.0):
        
        
        from sklearn.cluster import MeanShift, estimate_bandwidth
        
        pos_scaled = np.copy(positions)
        #pos_scaled[3:] = pos_scaled[3:]*scale_time
        
        labels = np.array([0])

        if len(pos_scaled) > 2:
            #X = np.array(zip(crus_x,crus_y,crus_z,crus_time))
            # #############################################################################
            # Compute clustering with MeanShift

            bandwidth = estimate_bandwidth(pos_scaled, quantile=0.2, n_samples=len(pos_scaled))
    
            if (bandwidth <0.5):
                bandwidth = threshold
                        
            ms = MeanShift(bandwidth=bandwidth, bin_seeding=True)
            ms.fit(pos_scaled)
            labels = ms.labels_
            cluster_centers = ms.cluster_centers_

        elif len(pos_scaled)==2:
            temp = pos_scaled[0]
            temp -= pos_scaled[1]
            dist1 = np.linalg.norm(temp,2)
            if dist1 > threshold:
                labels = np.array([0,1])
            else:
                labels = np.array([0,0])
                
        crus_per_cluster = np.array([np.count_nonzero(labels == i) for i in range(labels.max()+1)])
        #crus_per_cluster = np.zeros(labels.max())
        
        return labels, crus_per_cluster 
    
    def get_spark_max(self, labels, cru_info, quantity = "cyto_ca2+"):
    
        values = np.zeros(len(set(labels)))
        
        counter = 0

        for cru_id in cru_info['cru_id'].unique():
            
            cru_df = cru_info[cru_info['cru_id'] == cru_id]
            label = labels[counter]
            if values[label] < cru_df[quantity].max():
                values[label] = cru_df[quantity].max()
            counter += 1
            
        return values
    
    def get_FDX(self, labels, cru_info, quantity_max, quantity = "cyto_ca2+", percent = 50.0):
    
        spark_times = np.zeros(len(set(labels)))
        
        counter = 0
        
        for cru_id in cru_info['cru_id'].unique():
            
            label = labels[counter]
            
            cru_df = cru_info[cru_info['cru_id'] == cru_id]
            
            cru_dfX = cru_df[cru_df["cyto_ca2+"] > (float(percent)/100.0)*quantity_max[label]]
            timeX = cru_dfX['time'].max() - cru_dfX['time'].min()
            
            if spark_times[label] < timeX:
                spark_times[label] = timeX
                
        return spark_times

    def __eventSQcounter(self, dataframe, eventduration = 20, getTimeIntervals = False):
        """ prvt fct
        returns the counter for events happened during opening times of eventduration
        from dataframes that are given from the cleftlogs
        if getTimeIntervals is true also write out time intervals and cleftnames to list
        
        Returns:
        counter: int of number of events (Quarks or Sparks)
        """
        counter = 0
        timeintervals = []
        # go through the clefts
        for cleft in dataframe.clefts.drop_duplicates():
            movingtime = min(dataframe[dataframe.clefts == cleft].time)
            begin = movingtime
            counter += 1
            # check for each time step if it is within eventduration
            for times in dataframe[dataframe.clefts == cleft].time:
                if (times < movingtime + eventduration):
                    movingtime = times
                    end = times
                else:
                    timeintervals.append((cleft, begin, end))
                    begin = times
                    movingtime = times
                    end = times
                    counter += 1
            timeintervals.append((cleft, begin, end))
        
        if getTimeIntervals:
            return counter, timeintervals
        else:
            return counter

    def __cleanOfeventSQcounter(self, oneRyR, sparks, timeinter, eventduration = 20):
        """
        if there are more sparks than events
        cleans out the sparks which happen during one event
        returns the number of events with at least one spark
        """
        allev, timeev = self.__eventSQcounter(oneRyR, eventduration=eventduration, getTimeIntervals=True)
        
        count = np.zeros(allev)
        for sprk in timeinter:
            pos = 0
            for ev in timeev:
                if (ev[0] != sprk[0]) and (ev[1] <= sprk[1]):
                    pos += 1
            
            count[pos] = 1
        return int(sum(count))        

    def getQuarkSparkInfo(self, channelDF, eventduration = 20, getCleftnumber = False,
                          getSparkTimeIntervals = False, getMaxChannelsDurations = False):
        """
        gets Quarks and Sparks given cleftlog DataFrame and optional number of clefts firing
        also returns np.array of peak bulk Ca for each cleft
        
        Returns:
        quarks, sparks: int of number of quarks/sparks in simulation and cell
        peak: array of peak Calcium
        FDHM: np.array full duration at half maximum
        FD90: np.array full duration at 90% of max (similar to APD90)
        celftn. : optional return: opening cleft number
        """
        oneRyR = channelDF[channelDF.openRyR > 0]
        moreRyR = oneRyR[oneRyR.openRyR > 1]
        
        allev = self.__eventSQcounter(oneRyR, eventduration=eventduration)
        sparks, timeinter = self.__eventSQcounter(moreRyR, eventduration=eventduration,
                                       getTimeIntervals=True)
        if allev < sparks:
            sparkgroups = self.__cleanOfeventSQcounter(oneRyR, sparks, timeinter, eventduration=eventduration)
            quarks = allev - sparkgroups
        else:
            quarks = allev - sparks
        
        peaks = np.zeros(len(timeinter))
        peak_times = np.zeros(len(timeinter))
        
        FDHM = np.array([])
        FD90 = np.array([])
        
        maxryr =  np.array([])
        durations = np.array([])
        
        ana = Analysis()
        
        # add possibility of getting rid of low FDWM values e.g.
        # FDHM = FDHM[FDHM > 0.9]
        for ind, t in enumerate(timeinter):
            # filter Ca time series 2ms before event to resolve upstroke
            series = moreRyR[(moreRyR.clefts == t[0]) & (moreRyR.time>=t[1] - 2) &
                                     (moreRyR.time<=t[2])].bulkCa
            times = moreRyR[(moreRyR.clefts == t[0]) & (moreRyR.time>=t[1] - 2) &
                            (moreRyR.time<=t[2])].time
            
            peaks[ind] = max(series)
            peak_times[ind] = moreRyR.loc[(moreRyR[(moreRyR.clefts == t[0]) & (moreRyR.time>=t[1])
                                                   & (moreRyR.time<=t[2])].bulkCa.idxmax())].time
            
            timediff = t[2] - t[1]
            maxnumryrs = max(moreRyR[(moreRyR.clefts==t[0]) & (moreRyR.time>=t[1] ) &
                                   (moreRyR.time<=t[2])].openRyR )
            
            durations = np.append(durations, timediff)
            maxryr = np.append(maxryr, maxnumryrs)
            
            # if the series is too short, sometimes single opening
            # events see more than 1 open RyR, -> don't compute FDHM then
            if len(series) <= 2:
                print("Attention: skipping FDHM computation, series too short!")
                continue
            
            # some fast events have the beginning of the spark upstroke after
            # the HM has been computed (IndexError thrown there) -> skip computation then
            # this problem should be fixed by the subtraction of 2ms to time threshold
            try:
                FDHM = np.concatenate((FDHM, ana.APD(times, series, start_time=t[1],
                                                     end_time=t[2])))
                FD90 = np.concatenate((FD90, ana.APD(times, series, start_time=t[1],
                                                     end_time=t[2], percent=90)))
            except IndexError:
                continue
        
        if getCleftnumber:
            return quarks, sparks, peaks, FDHM, FD90, len(oneRyR.clefts.drop_duplicates())
        elif getSparkTimeIntervals:
            # this option is meant for debugging/analysis option
            return quarks, sparks, peaks, FDHM, FD90, timeinter
        elif getMaxChannelsDurations:
            return quarks, sparks, peaks, FDHM, FD90, maxryr, durations
        else:
            return quarks, sparks, peaks, FDHM, FD90 

    def computeSparkBiomarkers(self, sim_dirs, params_vari, start_time, end_time,
                               sampling_dir = "sampling", dropShortFDHM=False):
        """Computes several Biomarkers for the Spark simulations
        (max_Vm, rest_Vm, max_dVdt, dome_Vm,
        mean and std of: APD50, APD90, Ca_peaks, Ca_dia,
        Ca_time_to_peak) 
        
        given a list of parameters params_vari, a list of simulation
        folders and starting and ending times and a sampling dir name
        
        if dropShortAPD True, then all APDs < 50 are dropped, if numerical value,
        then APDs < dropShortAPD are dropped.
        
        Args:
        - sim_dirs: list of simulation directories (usually from sampling)
        - params_vari: list of parameters to be plotted against
        - start_time
        - end_time
        - sampling_dir
        - dropShortFDHM: default False; drops FDHM<0.9 ms if True or <dropShortFDHM if int/float
        
        Returns:
        - spark_data_df: pandas DF
        """
        from . import extract
        data = np.zeros((0,len(params_vari)))
        spark_data_df = pd.DataFrame(data=data, index=[], columns=params_vari)
        
        for sim_dir in sim_dirs:
            
            f = extract.extract(sim_dir,sampling_dir)
            vari = f.readIonic()
            varm = f.readMass()
            para = f.getParameters()
            channelDF = f.getCleftChannelInformation()
            
            if (max(vari['time']) < end_time):
                import warnings
                warnings.warn("SimNr {} ends at {}, while end_time is {}. Skip!".
                              format(sim_dir, max(vari['time']), end_time)) 
                continue
            
            counter = len(spark_data_df) + 1
            spark_data_df.at[counter,'sim'] = sim_dir
            for param in para:
                if param in params_vari:
                    spark_data_df.at[counter,param] = float(para.iloc[0][param])
            
            # may be not needed
            # channels = f.crusInfo()
            # spark_data_df.at[counter,'RyRtot'] = sum(channels[0])
            # spark_data_df.at[counter,'LCCtot'] = sum(channels[1])
            
            # get sparks, quarks counts and mean peak Calcium
            quarks, sparks, peaks, FDHM, FD90=self.getQuarkSparkInfo(channelDF, eventduration=20)
            spark_data_df.at[counter, 'quarks'] = int(quarks)
            spark_data_df.at[counter, 'sparks'] = int(sparks)
            
            spark_data_df.at[counter, 'quarkspark_ratio'] = quarks/sparks
            spark_data_df.at[counter, 'Ca_peak_mean'] = peaks.mean()
            
            # FDHM = FDHM[FDHM > 0.9]
            if dropShortFDHM == True:
                FDHM = FDHM[FDHM > 0.9]
                FD90 = FD90[FD90 > 0.9]
            elif type(dropShortFDHM) == int or type(dropShortFDHM) == float:
                FDHM = FDHM[FDHM > dropShortFDHM]
                FD90 = FD90[FD90 > dropShortFDHM]
            
            spark_data_df.at[counter, 'FDHM_mean'] = FDHM.mean()
            spark_data_df.at[counter, 'FDHM_mean_std'] = FDHM.std()
            
            spark_data_df.at[counter, 'FD90_mean'] = FD90.mean()
            spark_data_df.at[counter, 'FD90_mean_std'] = FD90.std()            
            
            # for later, when processChannelInfo is done add an if condition to check
            # if channelInfo.csv exists in clefts directory to speed up the loading!
            # might not be needed
            chInfo_df = f.processCRUInfo()["openRyR_per_ms"]
            spark_data_df.at[counter, "openRyR_per_ms"] = float(chInfo_df)
            
            spark_data_df.at[counter, 'folder'] = "../{}/{}/".format(sampling_dir, sim_dir)
        
        
        return spark_data_df

    def plotOpenRyRtimeSeries(self, sim_dir, sampling_dir, saveFig=False):
        """
        Plot all time series of open RyRs of one simulation given:
        - sim_dir: number of simulation directory
        - sampling_dir: sampling subfolder
        """
        from . import extract
        
        f = extract.extract(sim_dir,sampling_dir)
        channelDF = f.getCleftChannelInformation()
        oneRyR = channelDF[channelDF.openRyR > 0]
        openclefts = oneRyR.clefts.drop_duplicates()
        
        numopenclefts = len(openclefts)
        
        fig, ax = plt.subplots(numopenclefts, 1, tight_layout=True, figsize=(14, 3*numopenclefts))
        for ind, cl in enumerate(openclefts):
            ax[ind].plot(oneRyR[oneRyR.clefts == cl].time, oneRyR[oneRyR.clefts == cl].openRyR)
            ax[ind].set_title(cl)
        
        plt.tight_layout()
        if saveFig == True:
            print("Attention: no saving path specified, saving into cwd")
            fig.savefig("./time-series-openryr_sim{}_from{}.pdf".format(sim_dir, sampling_dir))
        elif type(saveFig) == str:
            fig.savefig(saveFig)
        else:
            fig.show()        
        
        
    def sparkAroundCleft(self, spark_info, waitingTime = 30, radius = 1, channelDF = pd.DataFrame(), getArea = False):
        '''
        - gets sparks info with [cleft, start time, end time]
        - gets a waiting time 
        - can get a dataframe if this function is only used, for that getArea = True
        
        
        Calculates if around a certain cleft with a spark, appear other sparks in the Clefts in a sphere around them in the time 
        after the first spark
        
    
        '''
        cleftNo = spark_info[0]
        timeOfSpark = spark_info[1]
        
        
        if getArea:
            self.list_of_sparks = self.getQuarkSparkInfo(channelDF, getSparkTimeIntervals = True)
            
                
        
        cru_x = self.datafr[self.datafr.clefts == cleftNo].xcru[0]
        cru_y = self.datafr[self.datafr.clefts == cleftNo].ycru[0]
        cru_z = self.datafr[self.datafr.clefts == cleftNo].zcru[0]
        df = self.datafr[self.datafr.clefts != cleftNo]
        df = df[((df.xcru - cru_x)**2 + (df.ycru - cru_y)**2 + (df.zcru - cru_z)**2 <= radius**2) & 
                (df.time >= timeOfSpark) & (df.time <= timeOfSpark + waitingTime)]

        clef = df.clefts.drop_duplicates().values
        time = df.time.values
        sparks = [sp for sp in self.list_of_sparks if ((sp[0] in clef) & (sp[1] in time))]
                    
        return sparks
    
    def nextNeighbours(self, spark_info, waitingTime = 30, channelDF = pd.DataFrame(), getArea = False):
        '''
        - gets sparks info with [cleft, start time, end time]
        - gets a waiting time 
        - can get a dataframe if this function is only used, for that getArea = True
        
        Calculates if around a certain cleft with a spark, appear other sparks in the nearest Neighbourhood around them in the time 
        after the first spark
    
        '''
        cleftNo = spark_info[0]
        timeOfSpark = spark_info[1] 
        cleft_in_z = np.array([])
        dist_z = 0
        
        variance = 1.01
        radius = abs(self.datafr.xcru.drop_duplicates()[0] - self.datafr.xcru.drop_duplicates()[1]) * variance
        if len(self.datafr.zcru.drop_duplicates().values) > 1 : 
            dist_z = abs(self.datafr.zcru.drop_duplicates()[0] - self.datafr.zcru.drop_duplicates()[1]) * variance
        
        if getArea:
            self.list_of_sparks = self.getQuarkSparkInfo(channelDF, getAllClefts = True)
        
        cru_x = self.datafr[self.datafr.clefts == cleftNo].xcru[0]
        cru_y = self.datafr[self.datafr.clefts == cleftNo].ycru[0]
        cru_z = self.datafr[self.datafr.clefts == cleftNo].zcru[0]
        df = self.datafr[self.datafr.clefts != cleftNo]

        if dist_z > radius :
            cleft_in_z = df[(df.xcru == cru_x) & (df.ycru == cru_y) & ((df.zcru - cru_z)**2 <= dist_z**2) & 
                            (df.time >= timeOfSpark) & (df.time <= timeOfSpark + waitingTime)].clefts.drop_duplicates().values
            
        df = df[((df.xcru - cru_x)**2 + (df.ycru - cru_y)**2 + (df.zcru - cru_z)**2 <= radius**2) & 
                (df.time >= timeOfSpark) & (df.time <= timeOfSpark + waitingTime)]

        clef = np.concatenate((df.clefts.drop_duplicates().values, cleft_in_z), axis = 0)
        time = df.time.values
        sparks = [sp for sp in self.list_of_sparks if ((sp[0] in clef) & (sp[1] in time))]
                    
        return sparks

    def __recursive_cascade(self, spark_info, layer, radius, neighbour = False):
        '''
        - gets spark info with [cleft, start time, end time]
        - gets the layer of the cascade
        - gets a radius for the neighbouring clefts or the next neighbourhood
        
        
        calculates the cascade of sparks in a recursive way
        '''

        if neighbour == False:
            sparks = self.sparkAroundCleft(spark_info,radius = radius)
        else:
            sparks = self.nextNeighbours(spark_info)
        
        self.cascade.append([layer, spark_info])

        self.list_of_sparks = [sp for sp in self.list_of_sparks if not sp in sparks]
        for sparkinf in sparks:
            self.__recursive_cascade(spark_info = sparkinf, layer = layer + 1, 
                                                radius = radius, 
                                                neighbour = neighbour)
        return 1
    
    def sparkCascade(self, channelDF, radius = 1, neighbour = False):
        '''
         - gets a dataframe from a simulation folder
         - calculates all sparks 
         - calculates if there is a cascade of sparks and returns the sparks
         
         returns:
             array of:
         [layer in the cascade,(cleft, start time, end time of spark)]
         
         one cascade is from one [0,(...)] to the next [0,(...)]
         the ones marked with [1,(...)] are triggered by the last [0,(...)] and so on ...
        '''
        self.datafr = channelDF
        quarks, sparks, peaks, FDHM, FD90, timeinter = self.getQuarkSparkInfo(self.datafr, getSparkTimeIntervals = True)
        self.list_of_sparks = timeinter
        self.cascade = []
        length = len(self.list_of_sparks)
        layer = 0
        
        while self.list_of_sparks:
            start_time = [sp[1] for sp in self.list_of_sparks]
            spark = self.list_of_sparks[start_time.index(min(start_time))]
            self.list_of_sparks = [sp for sp in self.list_of_sparks if ((sp[0] != spark[0]) or 
                                                                        (sp[1] != spark[1]))]
            self.__recursive_cascade(spark, layer, radius, neighbour = neighbour)
            
        if length == len(self.cascade):    
            return self.cascade
        else:
            return -1
        
    def cascadeInCleft(self, channelDF, cleft, start_time, end_time, radius = 1, neighbour = False):
        '''
        - gets a dataframe 
        - gets a specfic cleft with a spark starting and ending time
        - gets a radius or the next neighbourhood
        
        Calculates the cascade starting from a certain cleft at a certain time
        
        returns:
             array of:
         [layer in the cascade,(cleft, start time, end time of spark)]
         
         one cascade is from one [0,(...)] to the next [0,(...)]
         the ones marked with [1,(...)] are triggered by the last [0,(...)] and so on ...
        '''
        
        self.list_of_sparks = self.getQuarkSparkInfo(self.datafr, getSparkTimeIntervals = True)
        self.cascade = list(self.list_of_sparks)
        self.datafr = channelDF
        
        self.__recursive([cleft, start_time, end_time], layer = 0, 
                               radius = radius, neighbour = neighbour)
        
        return self.cascade
    
    def makroSpark(self, cascade):
        '''
        - gets a cascade
        
        calculates if there is a macro spark 
        a macro spark is a cascade with more than two sparks 
        returns a list of macro cascades with their first spark as macro spark
        '''
        macro = []
        micro = []
        for spark in cascade:
            if spark[0] == 0:
                if len(micro) > 2:
                    macro.append(micro)
                micro = [spark]
            if spark[0] != 0:
                micro.append(spark)
        if len(micro) > 2:
            macro.append(micro)
        return macro 

# ----------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------

class Analysis:
    """
    Analyse AP simulations
    """

    def __init__ (self):
        """
        Constructor of the main analysis class (still empty)
        """

    def threshold_crossings(self, series, threshold, return_derivatives=False):
        """Find all up and down threshold crossings within the series using
        linear interpolation. (adapted from caostk)
        
        Args
        ----------
        series : {N,} ndarray
        Time series (or similar) in which to find threshold crossings.
        threshold : float
        The threshold value to search.
        return_derivatives : True
        If True, also returns two additional arrays with the first order
        local derivative approximation at the points.
        
        Returns
        -------
        up : ndarray
        Array of crossings from low to high
        down : ndarray
        Array of crossings going from high to low
        """
        series = np.asarray(series)
        if series.ndim != 1:
            raise ValueError("series must be one dimensional.")
        high = series >= threshold
        low = ~high
        up = np.flatnonzero(low[:-1] & high[1:])
        down = np.flatnonzero(high[:-1] & low[1:])
        del high, low
        
        diff_up = series[up+1] - series[up]
        diff_down = series[down+1] - series[down]
        up = up + (threshold - series[up]) / diff_up
        down = down + (threshold - series[down]) /  diff_down
        
        if not return_derivatives:
            return up.astype(int), down.astype(int)
        else:
            return up, down, diff_up, diff_down

# ----------------------------------------------------------------------------------------------

    def APDHM(self, times, series):
        """
        Compute the APD at half maximum of a given series using threshold_crossings and APD fct
        given the times and the series
        
        Parameters
        ----------
        - times
        - series
        
        Returns
        ----------
        - apd: ndarray 
        Array of the APDs
        """
        return self.APD(times, series, percent = 50)

# ----------------------------------------------------------------------------------------------
    
    def APD(self, times, series, percent = 50.0, start_time = -np.infty, end_time = np.infty,
            forceNoDAD = False, thresholdDAD = 50):
        """
        Compute the APD at (percent) % of a given series using the threshold_crossings fct
        given the times and the series. For instance, for percent = 90.0 the APD_90 values are computed.
        
        Parameters
        ----------
        - times
        - series
        - percent (default is 50%)
        - start/end _time: starting and ending point for APD computation
        - forceNoDAD: if DAD present, drop APDs < thresholdDAD ms
        - thresholdDAD: threshold where to drop DADs
        
        Returns
        ----------
        - apd: ndarray 
        Array of the APDs
        """        
        series = np.array(series)
        times = np.array(times)
        
        #extract values of series between start_time and end_time 
        series = series[times > start_time]
        times = times[times > start_time]
        series = series[times < end_time] 
        times = times[times < end_time]
        
        
        xm = np.max(series) - (np.max(series) - np.min(series))*percent/100.0
        up, down = self.threshold_crossings(series, xm)       
        # remove all array entries exceeding other array
        # first a down is not allowed
        if (down[0] < up[0]):
            down = down[1:]        
        while(up.size != down.size):
            if (up.size > down.size):
                up = up[:-1]
            else:
                import warnings
                warnings.warn("The time series might be discontinous.\nFound downstroke with no preceeding upstroke. \t This function might be unsuited for the APD computation")                
                down = down[:-1]
        apd = times[down] - times[up]
        if forceNoDAD:
            print("Dropping APDs lower than {} ms".format(thresholdDAD))
            apd = apd[apd > thresholdDAD]
        return apd

# ----------------------------------------------------------------------------------------------

    def get_DI(self, times, series, start_time = -np.infty, end_time = np.infty):
        """
        Compute the DI of a given series using the threshold_crossings fct
        given the times and the series.
        
        Parameters
        ----------
        - times
        - series
          
        Returns
        ----------
        - apd: ndarray 
        Array of the DIs
        """
        percent = 95.0
        series = -np.array(series)
        times = np.array(times)
        
        #extract values of series between start_time and end_time 
        series = series[times > start_time]
        times = times[times > start_time]
        series = series[times < end_time] 
        times = times[times < end_time]
        
        
        xm = (np.max(series) + np.min(series))*percent/100.0
        up, down = self.threshold_crossings(series, xm)       
        # remove all array entries exceeding other array
        # first a down is not allowed
        if (down[0] < up[0]):
            down = down[1:]        
        while(up.size != down.size):
            if (up.size > down.size):
                up = up[:-1]
            else:
                import warnings
                warnings.warn("The time series might be discontinous.\nFound downstroke with no preceeding upstroke. \t This function might be unsuited for the APD computation")                
                down = down[:-1]
        return (times[down] - times[up])

# ----------------------------------------------------------------------------------------------

    def get_max_Vm(self, times, series, start_time = -np.infty, end_time = np.infty):
        """
        Computes the maximum voltage of the membrane potential from a series of time points and corresponding membrane potential values
        
        Parameters
        ----------
        - times
        - series
        
        Returns
        ----------
        - maximum of series (scalar float value)
        """
        
        series = np.array(series)
        times = np.array(times)
        
        #extract values of series between start_time and end_time 
        series = series[times > start_time]
        times = times[times > start_time]
        series = series[times < end_time] 
        times = times[times < end_time]
        
        
        return np.max(series)

# ----------------------------------------------------------------------------------------------

    def get_rest_Vm(self, times, series, start_time = -np.infty, end_time = np.infty):
        """
        Computes the resting voltage of the membrane potential from a series of time points and corresponding membrane potential values
        
        Parameters
        ----------
        - times
        - series
        
        Returns
        ----------
        - minimum of series (scalar float value)
        """
        
        series = np.array(series)
        times = np.array(times)
        
        #extract values of series between start_time and end_time 
        series = series[times > start_time]
        times = times[times > start_time]
        series = series[times < end_time] 
        times = times[times < end_time]
        
        return np.min(series)    

# ----------------------------------------------------------------------------------------------
    
    def get_max_dVdt(self, times, series, scheme = 3, start_time = -np.infty, end_time = np.infty):
        """
        Compute maximum dVdt of a given series by using a first or second order approximation of the derivative
        
        Parameters
        ----------
        - times
        - series
        - order (optional)
          
        Returns
        ----------
        - the max direvative of the membrane potential dVdt
        """
        
        series = np.array(series)
        times = np.array(times)
        
        #extract values of series between start_time and end_time 
        series = series[times > start_time]
        times = times[times > start_time]
        series = series[times < end_time] 
        times = times[times < end_time]
        
        max_dVdt = -1.0
                
        for i in range(len(series) - 2):
            
            dVdt = -1.0
            
            if scheme == 1: # first order scheme
                dVdt = (series[i+1] -series[i])/(times[i+1] - times[i])
            if scheme == 2: # smoothed, first order scheme
                dVdt = (series[i+2] -series[i])/(times[i+2] - times[i])
            if scheme == 3: # second order scheme
                dVdt = ( - 3.0*series[i] + 4.0*series[i+1] - series[i+2])/(times[i+2] - times[i])
            
            if (dVdt > max_dVdt):
                max_dVdt = dVdt
  
        return max_dVdt
    
# ----------------------------------------------------------------------------------------------

    def get_dome_Vm(self, times, series, time_threshold = 100.0, smooth = 50, start_time = -np.infty, end_time = np.infty):
        """
        Computes domes in a membrane potential series and returns the corresponding Vm values and timepoints
        
        Parameters
        ----------
        - times
        - series
        - smooth (optional, determines the size of neighborhood for local peak search)
        
        Returns
        ----------
        - value of dome Vm
        """
        
        series = np.array(series)
        times = np.array(times)
        
        #extract values of series between start_time and end_time 
        series = series[times > start_time]
        times = times[times > start_time]
        series = series[times < end_time] 
        times = times[times < end_time]
        
        times_maxima = times[argrelextrema(series, np.greater_equal, order=smooth)]
        Vm_maxima = series[argrelextrema(series, np.greater_equal, order=smooth)]
        
        #Vm_maxima = Vm_maxima[Vm_maxima>0.0]
        
        if Vm_maxima[0] < Vm_maxima[1]:
            times_maxima = np.delete(times_maxima, 0)
            Vm_maxima = np.delete(Vm_maxima, 0)
        
        Vm_domes = []
        t_domes = []
        
        factor = 0.7
        
        #for i in range(len(Vm_maxima)/2 - 2):
            #if (factor*Vm_maxima[2*i] > Vm_maxima[2*i+1]) and (Vm_maxima[2*i+1] < factor*Vm_maxima[2*i+2]):
                #if (Vm_maxima[2*i+1] > 0.0) and (times_maxima[2*i] - times_maxima[2*i+1] < time_threshold):
                    #Vm_domes.append(Vm_maxima[2*i+1])
                    #t_domes.append(times_maxima[2*i+1])
                
                    
        for i in range(len(Vm_maxima)):
            if (factor*Vm_maxima.max() > Vm_maxima[i]) and (Vm_maxima[i] > 0.0):
                if (len(Vm_domes)==0):
                    Vm_domes.append(Vm_maxima[i])
                    t_domes.append(times_maxima[i])
                if (len(Vm_domes)>0):
                    if (times_maxima[i] - t_domes[-1] > time_threshold):
                        Vm_domes.append(Vm_maxima[i])
                        t_domes.append(times_maxima[i])
                else:
                    print("Time interval between domes is too small, increase smoothing factor!")
            else:
                import warnings
                warnings.warn("Implement exception here!")   
        

        #TODO: revise that code and include exceptions
        Vm_domes = np.array(Vm_domes)
        t_domes = np.array(t_domes)

        return t_domes, Vm_domes

# ----------------------------------------------------------------------------------------------

    def get_Ca_peaks(self, times, series, smooth = 1000, time_threshold = 100, start_time = -np.infty, end_time = np.infty):
        """
        Computes peaks in series and returns the corresponding peak values and timepoints
        
        Parameters
        ----------
        - times
        - series
        - smooth (optional, determines the size of neighborhood for local peak search)
        - threshold (optional, minimum time between two peaks)
        
        Returns
        ----------
        Tuple of arrays of peak times and peak values
        """
        
        series = np.array(series)
        times = np.array(times)
        
        #extract values of series between start_time and end_time 
        series = series[times > start_time]
        times = times[times > start_time]
        series = series[times < end_time] 
        times = times[times < end_time]
        
        times_maxima = times[argrelextrema(series, np.greater_equal, order=smooth)]
        peaks_maxima = series[argrelextrema(series, np.greater_equal, order=smooth)]
        #times_minima = times[argrelextrema(series, np.less_equal, order=smooth)]
        #peaks_minima = series[argrelextrema(series, np.less_equal, order=smooth)]
     
     
        min_diff = min([abs(j-i) for i,j in zip(times_maxima, times_maxima[1:])]) 

        if (min_diff > time_threshold):
            import warnings
            warnings.warn("The time series might be discontinous.\nFound more than peak in a small time intervall. \t Please increase the smoothing factor!")                
        return times_maxima, peaks_maxima


    def computeBiomarkers(self, sim_dirs, params_vari, start_time, end_time,
                          sampling_dir = "sampling", dropShortAPD=False):
        """Computes several Biomarkers 
        (max_Vm, rest_Vm, max_dVdt, dome_Vm,
        mean and std of: APD50, APD90, Ca_peaks, Ca_dia,
        Ca_time_to_peak) 
        
        given a list of parameters params_vari, a list of simulation
        folders and starting and ending times and a sampling dir name
        
        if dropShortAPD True, then all APDs < 50 are dropped, if numerical value,
        then APDs < dropShortAPD are dropped.
        
        Args:
        - sim_dirs: list of simulation directories (usually from sampling)
        - params_vari: list of parameters to be plotted against
        - start_time
        - end_time
        - sampling_dir
        - dropShortAPD: default is False, can be true or a ms numerical value
        
        Returns:
        - sim_data_df: pandas DF
        """

        from . import extract

        data = np.zeros((0,len(params_vari)))
        sim_data_df = pd.DataFrame(data=data, index=[], columns=params_vari)
        counter = 0

        for sim_dir in sim_dirs:
            
            f = extract.extract(sim_dir,sampling_dir)
            vari = f.readIonic()
            varm = f.readMass()
            para = f.getParameters()
            
            if (max(vari['time']) < end_time):
                import warnings
                warnings.warn("SimNr %s ends at %s, while end_time is %s. Skip!" % (sim_dir, max(vari['time']), end_time)) 
                continue
                
            counter = len(sim_data_df) + 1
            sim_data_df.at[counter,'sim'] = sim_dir
            for param in para:
                if param in params_vari:
                    sim_data_df.at[counter,param] = float(para.iloc[0][param])

            sim_data_df.at[counter,'max_Vm'] = self.get_max_Vm(list(vari['time']),list(vari['Vm']),start_time=start_time,end_time=end_time)
            sim_data_df.at[counter,'rest_Vm'] = self.get_rest_Vm(list(vari['time']),list(vari['Vm']),start_time=start_time,end_time=end_time)
            sim_data_df.at[counter,'max_dVdt'] = self.get_max_dVdt(list(vari['time']),list(vari['Vm']),3,start_time=start_time,end_time=end_time)
            times_dome, Vm_domes = self.get_dome_Vm(list(vari['time']),list(vari['Vm']),start_time=start_time,end_time=end_time)
    
            if (len(Vm_domes)>0):
                sim_data_df.at[counter,'dome_Vm'] = Vm_domes.mean()
            else:
                print("No dome found!")

            # compute APDs according to dropping or not of the short APDs
            if not dropShortAPD:
                APD50 = self.APD(list(vari['time']),list(vari['Vm']),percent=50,start_time=start_time,end_time=end_time, forceNoDAD=False)
                APD90 = self.APD(list(vari['time']),list(vari['Vm']),percent=90,start_time=start_time,end_time=end_time, forceNoDAD=False)
                
            elif dropShortAPD == True:
                APD50 = self.APD(list(vari['time']),list(vari['Vm']),percent=50,start_time=start_time,end_time=end_time, forceNoDAD=True)
                APD90 = self.APD(list(vari['time']),list(vari['Vm']),percent=90,start_time=start_time,end_time=end_time, forceNoDAD=True)
            elif type(dropShortAPD) == int or type(dropShortAPD) == float:
                APD50 = self.APD(list(vari['time']),list(vari['Vm']),percent=50,start_time=start_time,end_time=end_time, forceNoDAD=True, thresholdDAD=dropShortAPD)
                APD90 = self.APD(list(vari['time']),list(vari['Vm']),percent=90,start_time=start_time,end_time=end_time, forceNoDAD=True, thresholdDAD=dropShortAPD)
            
            sim_data_df.at[counter,'APD50_mean'] = APD50.mean()
            sim_data_df.at[counter,'APD50_std'] = APD50.std()
            
            sim_data_df.at[counter,'APD90_mean'] = APD90.mean()
            sim_data_df.at[counter,'APD90_std'] = APD90.std()
    
            times_peaks, Ca_peaks = self.get_Ca_peaks(list(vari['time']),list(vari['Ca_i']),start_time=start_time,end_time=end_time)
            sim_data_df.at[counter,'Ca_peak_mean'] = Ca_peaks.mean()
            sim_data_df.at[counter,'Ca_peak_std'] = Ca_peaks.std()
    
            times_peaks, Ca_minima = self.get_Ca_sys_minima(list(vari['time']),list(vari['Ca_i']),start_time=start_time,end_time=end_time)
            sim_data_df.at[counter,'Ca_dia_mean'] = Ca_minima.mean()
            sim_data_df.at[counter,'Ca_dia_std'] = Ca_minima.std()
            
            times_to_peak = self.get_Ca_times_to_peak(list(vari['time']),list(vari['Ca_i']),start_time=start_time, end_time=end_time)
    
            sim_data_df.at[counter,'Ca_time_to_peak_mean'] = times_to_peak.mean()
            sim_data_df.at[counter,'Ca_time_to_peak_std'] = times_to_peak.std()
            
            sim_data_df.at[counter,'Na_i'] = np.array(vari['Na_i'])[-100:].mean()

            cyto_vol = float(para['xmax']*para['ymax']*para['zmax'])
            bspecial_conc = np.array(varm['SpecialBm'])/(float(cyto_vol))
            ca_exp =  self.get_experimental_Ca(bspecial_conc,para)
            
            times_peaks, Ca_exp_peaks = self.get_Ca_peaks(list(varm['Time']),ca_exp,start_time=start_time,end_time=end_time)
            sim_data_df.at[counter,'Ca_exp_peak_mean'] = Ca_exp_peaks.mean()
            sim_data_df.at[counter,'Ca_exp_peak_std'] = Ca_exp_peaks.std()
    
            times_peaks, Ca_exp_minima = self.get_Ca_sys_minima(list(varm['Time']),ca_exp,start_time=start_time,end_time=end_time)
            sim_data_df.at[counter,'Ca_exp_dia_mean'] = Ca_exp_minima.mean()
            sim_data_df.at[counter,'Ca_exp_dia_std'] = Ca_exp_minima.std()
            
            sim_data_df.at[counter, 'folder'] = "../{}/{}/".format(sampling_dir, sim_dir)

        
        return sim_data_df
    
    def computeCleftGeometry(self, input_df):
        
        """Computes geometric cleft properties and writes them into a given biomarker data frame
    
        
        Args:
        - input_df that contains simulation folders
        
        Returns:
        - sim_data_df: pandas DF
        """

        from . import extract
        
        sim_data_df = input_df.copy()

        for index, row in sim_data_df.iterrows():
    
            folder = row['folder'][3:]
    
            ex = extract.extract(folder)
            ryr_numbers, lcc_numbers, ryr_locations, lcc_locations = ex.crusInfo(getLocations = True)
            
            dists_nn = []
            dists_4nn = []
            for ryr_location in ryr_locations:
                
                for i in range(0,len(ryr_location[0])):
                    dists = []
                    x1 = ryr_location[0][i]
                    y1 = ryr_location[1][i]
                    distance = np.infty
                    for j in range(0,len(ryr_location[0])):
                        if (i==j):
                            continue
                        x2 = ryr_location[0][j]
                        y2 = ryr_location[1][j]
                        distance = np.sqrt((x1-x2)**2 + (y1-y2)**2)
                        dists.append(distance)
                    
                    dists.sort()
                    

                    dists_nn.append(dists[0])
                    mean_4nn = sum(dists[0:4])/4.0
                    dists_4nn.append(mean_4nn)
    
    
            dists_nn = np.array(dists_nn)
            dists_4nn = np.array(dists_4nn)

            sim_data_df.at[index,"nn_mean"] = dists_nn.mean()*1000.0
            sim_data_df.at[index,"nn_std"] = dists_nn.std()*1000.0
            sim_data_df.at[index,"4nn_mean"] = dists_4nn.mean()*1000.0
        
        return sim_data_df
    
    def get_experimental_Ca(self, fluo4_conc, parameters):
        """Computes the experimental calcium concentration, which is inferred from a fluorescent dye as fluo4. 
        In cleftdyn fluo4 is represented by "Bspecial". 
        
        Args:
        - series of the buffer concentrations in this case 'Bspecial'.
        - parameters that have been used for the simulation 
        
        Returns:
        - series of experimental calcium
        """
        
        fluo4_conc = np.array(fluo4_conc)
        
        Fmax = float(parameters['Bspecial_tot'])
        Fmin = 0.0
        Kd = float(parameters['kspecial_minus'])/float(parameters['kspecial_plus'])
                
        ca_exp = Kd*(fluo4_conc- Fmin)/(Fmax - fluo4_conc)
        
        return ca_exp

    def get_Ca_sys_minima(self, times, series, smooth = 1000, time_threshold = 100, start_time = -np.infty, end_time = np.infty):
        """
        Computes minima during systole in series and returns the corresponding values and timepoints
        
        Parameters
        ----------
        - times
        - series
        - smooth (optional, determines the size of neighborhood for local minima search)
        - threshold (optional, minimum time between two minima)
        
        Returns
        ----------
        Tuple of arrays of peak times and minima values
        """
        
        series = np.array(series)
        times = np.array(times)
        
            
        times_, minima_ = self.get_Ca_peaks(times, -series, smooth = smooth, time_threshold = time_threshold, start_time = start_time, end_time = end_time)
        return times_, -minima_

# ----------------------------------------------------------------------------------------------

    def get_Ca_times_to_peak(self, times, series, smooth = 1000, time_threshold = 100, start_time = -np.infty, end_time = np.infty):
        """
        Computes the time needed from a calcium in diastole minima to the following maxima
        
        Parameters
        ----------
        - times
        - series
        - smooth (optional, determines the size of neighborhood for local minima search)
        - threshold (optional, minimum time between two minima/peaks)
        
        Returns
        ----------
        Arrays of with times to peak
        """
        
        series = np.array(series)
        times = np.array(times)
        
        t_peaks, peaks = self.get_Ca_peaks(times, series, smooth = smooth, time_threshold = time_threshold, start_time = start_time, end_time = end_time)
        t_min, minima = self.get_Ca_sys_minima(times, series, smooth = smooth, time_threshold = time_threshold, start_time = start_time, end_time = end_time)
        
        if t_peaks[0] < t_min[0]:
            t_peaks = np.delete(t_peaks, 0)
            
        times_to_peak = []
        
        for i in range(min(len(t_peaks),len(t_min)-1)):
            if ((t_peaks[i] > t_min[i]) and (t_peaks[i] < t_min[i+1])):
                times_to_peak.append(t_peaks[i] -t_min[i])
            else:
                import warnings
                warnings.warn("Implement exception here!")   
                
        return np.array(times_to_peak)

    def simFoldersHighSTD(self, sim_df, obj = "APD90_std", mid_val=None):
        """
        Returns np.array of simulation folders with obj values greater than mid value
        
        Parameters
        ----------
        sim_df: simulation DF
        obj: objective (default the APD90_std, since fct meant to study alternans)
        mid_val: optional externally fixed threshold value
        
        Returns
        ----------
        np.array with folder names
        """
        if not mid_val:
            mid_val = (max(sim_df[obj]) + min(sim_df[obj]) )/2
        else:
            mid_val = mid_val
        
        greater = sim_df[sim_df[obj] > mid_val]                
        sim_folders = greater["folder"]
        
        return np.array(sim_folders)

    def paramMeanUniformDistribution(self, params_vari, params_ranges):
        '''
        Returns dictionary of mean values and a chaospy uniform distribution
        
        Parameters
        ----------
        - params_vari: list of varied parameters
        - params_ranges: dictionary of parameter ranges        
        
        Returns
        ----------
        - dictionary of parameter mean values
        - chaospy uniform distribution             
        '''
        params_mean = dict()
        dist_args = []
        for pars in params_vari:
            params_mean[pars] = 0.5*(params_ranges[pars][0]+params_ranges[pars][1])
            dist_args.append(cp.Uniform(params_ranges[pars][0],params_ranges[pars][1]))
        
        dist = cp.J(*dist_args)
        return params_mean, dist

    def regressionFit(self, sim_data_df, params_vari, distr, polyOrder = 3,
                      objective = "APD50_mean", printInfo = False):
        '''
        Returns approx fct from regression fit given the biomarker DF,
        the varied parameters, the parameter cp distribution the polynomial order 
        and the fit objective.
        
        Parameters
        ----------
        - sim_data_df: Biomarker DF
        - params_vari: list of varied parameters
        - distr: uniform distribution of the varied parameter range
        - polyOrder: order of the polynomial to be fitted
        - objective: which value to analyse regression with
        - printInfo: default False, if True prints additional informations
        
        Returns
        ----------
        - approximated fit function
        '''
        # set objective on y and parameters on x
        ydata = np.array(list(sim_data_df[objective]))
        xdata = np.array(sim_data_df[params_vari])
        if printInfo:
            print("""Information about fitting:
            Number of data points: {}            
            Number of Parameters: {}
            Polynomial degree: {}
            Number of polynomial coefficients that have to be calculated: {}"""
                  .format(len(ydata), len(xdata.T), polyOrder,
                          int(sp.special.binom(len(xdata.T) + polyOrder,polyOrder))))
        # compute regression and fit
        orth_poly = cp.orth_ttr(polyOrder,distr)
	
        # compute experimental matrix
        A_exp = np.zeros((len(xdata),len(orth_poly)))
        
        for i in range(len(xdata)):
            for j in range(len(orth_poly)):
                A_exp[i,j] = orth_poly[j](*xdata[i])
        
        # calculation of the pseudo inverse
        A_inv = np.linalg.pinv(A_exp)
        coeffs = np.dot(A_inv,ydata)
        
        
        # calculation of the approximate function
        func_approx = coeffs[0]*orth_poly[0]
        for i in range(1,len(orth_poly)):
            func_approx += coeffs[i]*orth_poly[i]
        
        return func_approx


    def regressionFit_cp(self, sim_data_df, params_vari, distr, polyOrder = 3,
                      objective = "APD50_mean", printInfo = False):
        '''
        Returns approx fct from regression fit of chaospy given the biomarker DF,
        the varied parameters, the parameter cp distribution the polynomial order 
        and the fit objective.
        
        Parameters
        ----------
        - sim_data_df: Biomarker DF
        - params_vari: list of varied parameters
        - distr: uniform distribution of the varied parameter range
        - polyOrder: order of the polynomial to be fitted
        - objective: which value to analyse regression with
        - printInfo: default False, if True prints additional informations
        
        Returns
        ----------
        - chaospy approximated fit function
        '''
        # set objective on y and parameters on x
        ydata = np.array(list(sim_data_df[objective]))
        xdata = np.array(sim_data_df[params_vari].T)
        if printInfo:
            print("""Information about fitting:
            Number of data points: {}            
            Number of Parameters: {}
            Polynomial degree: {}
            Number of polynomial coefficients that have to be calculated: {}"""
                  .format(len(ydata), len(xdata), polyOrder,
                          int(sp.special.binom(len(xdata) + polyOrder,polyOrder))))
        # compute regression and fit
        # rule='T': Ridge Regression/Tikhonov Regularization Order
        orth_poly = cp.orth_ttr(polyOrder,distr)
        func_approx = cp.fit_regression(orth_poly, xdata, ydata, rule = 'T')
        return func_approx
    
    def lasso_regression(self, sim_data_df, params_vari, distr, alphas, polyOrder = 3, 
                      objective = "APD50_mean", kfold = 10, printInfo = False):
        '''
        Returns approx fct from regression fit with the lasso fit method given the biomarker DF,
        the varied parameters, the parameter cp distribution the polynomial order 
        and the fit objective.
        
        Parameters
        ----------
        - sim_data_df: Biomarker DF
        - params_vari: list of varied parameters
        - distr: uniform distribution of the varied parameter range
        - polyOrder: order of the polynomial to be fitted
        - alphas: list of alphas, alpha penalizes the L1-norm of polynomial coefficients
        - objective: which value to analyse regression with
        - printInfo: default False, if True prints additional informations
        
        Returns
        ----------
        - a disctionary with
            - chaospy approximated fit function
            - alpha that maximazes the cross validation score
            - cross validation score
            - coefficient of determination (R2)
        '''
        
        
        
        orth_poly = cp.orth_ttr(polyOrder,distr)
               
        ydata = np.array(list(sim_data_df[objective]))
        xdata = np.array(sim_data_df[params_vari])
        
        res = dict()
        
        if printInfo:
            print("""Information about fitting:
            Number of data points: {}            
            Number of Parameters: {}
            Polynomial degree: {}
            Number of polynomial coefficients that have to be calculated: {}"""
                  .format(len(ydata), len(xdata.T), polyOrder,
                          int(sp.special.binom(len(xdata.T) + polyOrder,polyOrder))))
        
        data = sim_data_df[[*params_vari,objective]]
        
        for i in range(0,len(orth_poly)):
            colname = 'phi_%d'%i

            #power of 1 is already there
            #new var will be x_power
            temp = np.zeros(len(xdata))
            for j in range(len(xdata)):
                temp[j] = orth_poly[i](*xdata[j])
            data[colname] = temp
       
        predictors = ['phi_%d'%i for i in range(len(orth_poly))]
        
        best_alpha = 0.0
        best_cross_val_score = -np.infty
        r2_score = -np.infty
  
        for alpha in alphas:
            val = 0.0
            r2_val = 0.0
            if alpha == 0.0:
                linreg = LinearRegression(normalize=True)
                reg = linreg.fit(data[predictors],data[objective])
                val = cross_val_score(linreg,data[predictors],data[objective],cv=kfold).mean()
                r2_val = reg.score(data[predictors],data[objective])
            else:
                lassoreg = Lasso(alpha=alpha,normalize=True, max_iter=1e5)
                lassoreg.fit(data[predictors],data[objective])
                val = cross_val_score(lassoreg,data[predictors],data[objective],cv=kfold).mean()
                r2_val = lassoreg.score(data[predictors],data[objective])
            
            if (val > best_cross_val_score):
                best_cross_val_score = val
                best_alpha = alpha
                r2_score = r2_val
            
        res['alpha'] = best_alpha
        res['cross_val_score'] = best_cross_val_score
        res['r2_score'] = r2_score
        
        lassoreg = Lasso(alpha = best_alpha, normalize = True, max_iter = 1e5)   
        lassoreg.fit(data[predictors],data[objective])
        coeffs = lassoreg.coef_
        
        func_approx = coeffs[0]*orth_poly[0] + lassoreg.intercept_*orth_poly[0]
        for i in range(1,len(orth_poly)):
            func_approx += coeffs[i]*orth_poly[i]
        
        res['func_approx'] = func_approx
        
        return res
    

    def calculate_regression_error(self, sim_data_df, params_vari, func_approx,
                      objective = "APD50_mean", printInfo = False):
        
        '''
        Returns the L^2, L^{\infty} and the relative error of the approx fct from regression fit of chaospy given
        with resprect to the datapoints of the objective.
        
        Parameters
        ----------
        - sim_data_df: Biomarker DF
        - params_vari: list of varied parameters
        - func_approx: function that was used for model approximation
        - objective: which value to analyse regression with
        
        Returns
        ----------
        - L2 and Linf error
        '''
        # set objective on y and parameters on x
        ydata = np.array(list(sim_data_df[objective]))
        xdata = np.array(sim_data_df[params_vari])
        
        L2_err = 0.0
        Linf_err = 0.0
        rel_err = 0.0
        
        for i in range(0,len(xdata)):
            L2_err += (func_approx(*xdata[i]) - ydata[i])**2
            temp = abs(func_approx(*xdata[i]) - ydata[i])
            rel_err += abs((func_approx(*xdata[i]) - ydata[i])/ydata[i])
            if (Linf_err < temp):
                Linf_err = temp
        L2_err = np.sqrt(L2_err/float(len(ydata)))
        rel_err = rel_err/float(len(ydata))
        
        return L2_err, Linf_err, rel_err 
    
    def calculate_LOO_error(self, sim_data_df, params_vari, distr, polyOrder = 3,
                            objective = "APD50_mean", printInfo = False, analytical = True):
        '''
        Returns the leave-one out (LOO) error for cross validation as well as the Q2 regression error estimator
        
        with resprect to the datapoints of the objective.
        
        Parameters
        ----------
        - sim_data_df: Biomarker DF
        - params_vari: list of varied parameters
        - distr: uniform distribution of the varied parameter range
        - polyOrder: order of the polynomial to be fitted
        - objective: which value to analyse regression with
        
        Returns
        ----------
        - L2 leave-one out error
        - Linf leave-one out error 
        - relative leave-one out error 
        - Q2 regression error estimator
        '''
        
        # set objective on y and parameters on x
        ydata = np.array(list(sim_data_df[objective]))
        xdata = np.array(sim_data_df[params_vari])
        
        # construct orthonormal basis
        orth_poly = cp.orth_ttr(polyOrder,distr)
        
        # compute experimental matrix
        A_exp = np.zeros((len(xdata),len(orth_poly)))
        
        for i in range(len(xdata)):
            for j in range(len(orth_poly)):
                A_exp[i,j] = orth_poly[j](*xdata[i])
                
        # calculation of the pseudo inverse
        A_inv = np.linalg.pinv(A_exp)
        coeffs = np.dot(A_inv,ydata)
        
        # calculation of the approximate function 
        
        func_approx = coeffs[0]*orth_poly[0]
        for i in range(1,len(orth_poly)):
            func_approx += coeffs[i]*orth_poly[i]
            
        rel_err = 0.0
        L2_err = 0.0
        Linf_err = 0.0
       
        if (analytical):
            # fast analytical computation
            
            H = np.dot(A_exp,A_inv)
            for i in range(0,len(xdata)):
                             
                h = H[i,i]                         
                L2_err += ((func_approx(*xdata[i]) - ydata[i])/(1.0 - h))**2
                rel_err += abs((func_approx(*xdata[i]) - ydata[i])/(ydata[i] - ydata[i]*h))
                abs_err = abs((func_approx(*xdata[i]) - ydata[i])/(1.0 - h))
                if (Linf_err < abs_err):
                    Linf_err = abs_err
                            
            L2_err = np.sqrt(L2_err/float(len(ydata)))
            rel_err = rel_err/float(len(ydata))
            
        else:
            # slow direct computation for comparison
             
            for k in range(0,len(xdata)):
                
                ydata_oo = np.copy(ydata)
                ydata_oo = np.delete(ydata_oo,k,0)
                xdata_oo = np.copy(xdata)
                xdata_oo = np.delete(xdata_oo,k,0)
                
                A_exp_oo = np.zeros((len(xdata_oo),len(orth_poly)))
                
                for i in range(len(xdata_oo)):
                    for j in range(len(orth_poly)):
                        A_exp_oo[i,j] = orth_poly[j](*xdata_oo[i])
                        
                # calculation of the pseudo inverse
                A_inv_oo = np.linalg.pinv(A_exp_oo)
                coeffs_oo = np.dot(A_inv_oo,ydata_oo)
                
                # calculation of the approximate function
                func_approx_oo = coeffs_oo[0]*orth_poly[0]
                for i in range(1,len(orth_poly)):
                    func_approx_oo += coeffs_oo[i]*orth_poly[i]
                
                
                L2_err += (func_approx_oo(*xdata[k]) - ydata[k])**2
                abs_err = abs(func_approx_oo(*xdata[k]) - ydata[k])
                rel_err += abs((func_approx_oo(*xdata[k]) - ydata[k])/ydata[k])
                if (Linf_err < abs_err):
                    Linf_err = abs_err
                    
            L2_err = np.sqrt(L2_err/float(len(ydata)))
            rel_err = rel_err/float(len(ydata))
        
        # computation of the Q2 error estimator
        Q2 = 0.0
        var_obj = ydata.var()
        if (var_obj > 0.0):
            Q2 = 1.0 -  L2_err*L2_err/var_obj
        else:
            import warnings
            warnings.warn("""\n Variance of the data is zero!""")
            
        return L2_err, Linf_err, rel_err, Q2 

    def slope_error1DFit(self, sim_df, params_vari, params_ranges, obj, getaxes = False):
        '''
        Returns slope and error of slope fitting from 1D first order polynomial fit to objective 
        given the simulation data, the objective of the fit and 
        the parameter vs which to fit.
        Default arg getaxes allows to get the lists of x, y and fit values for each parameter
        
        Args
        ----------
        - sim_df: simulation DF
        - params_vari: list of parameters
        - params_ranges: ranges of parameter (for normalisation)
        - obj: objective of fit
        - getaxes=False: default arg to output also axes lists (for plotting)
        
        Returns
        ----------
        slope: dict of fitted slopes
        err: dict of errors
        x, y, fit: np arrays of plotting values
        '''
        slope = {}
        err = {}
        x={}
        y={}
        fit={}
        
        for pars in params_vari:
            xtmp=sim_df[pars]/max(params_ranges[pars])
            ytmp=sim_df[obj]
            f, cov = np.polyfit(xtmp, ytmp, deg=1, cov=True)
            
            slope[pars] = f[0]
            err[pars] = np.sqrt(np.diag(cov))[0]
            x[pars] = xtmp
            y[pars] = ytmp
            fit[pars] = np.poly1d(f)
        
        slope = pd.DataFrame(slope, index=[0])
        err = pd.DataFrame(err, index=[0])
        x = pd.DataFrame(x)
        y = pd.DataFrame(y)
        if slope.isnull().values.any():
            print("Attention: the fit has some NaN value")
        if getaxes:
            return slope, err, x, y, fit
        else:
            return slope, err

    def kNNRegression(self, sim_data_df, obj, param, k=15, get_error=False, steps=500):
        '''
        Returns knn prediction on np.linspace() and eventually error 
        given x and y training data from objective and parameter
        
        Parameters
        ----------
        - sim_data_df: dataframe to analyse
        - obj: objective to be analysed
        - param: parameter to consider
        - k: number of NN to consider (can be set to "optimal" or "optimal-force")
        - get_error: output also the root mean square error
        - steps: adjust number of steps of the linspace
        
        Returns
        ----------
        x_predict, y_predict, error (optional)
        '''
        from sklearn import neighbors
        from sklearn.metrics import mean_squared_error
        
        y = sim_data_df[obj]
        x = np.array(sim_data_df[param])
        x = x.reshape(len(sim_data_df[param]), 1)
        T = np.linspace(min(x), max(x), steps)
        
        if "opt" in str(k):
            err = []
            for i in range(1,150):
                knn = neighbors.KNeighborsRegressor(i, weights="uniform", algorithm="kd_tree")
                y_ = knn.fit(x, y).predict(T)
                err.append(np.sqrt(mean_squared_error(T, y_)))
            optk = err.index(min(err))+1
            
            # check if optimisation is failing and if it was forced 
            if optk > 40 and "force" in k:
                k = optk 
                print("forcing k to {}".format(k))               
            elif optk > 40 and not "force" in k:
                k = 15
                print("keep k at {} (avoided forcing to {})".format(k, optk))
            else:
                k = optk
                print("k optimised to {}".format(k))
        
        knn = neighbors.KNeighborsRegressor(k, weights="uniform",
                                            algorithm="kd_tree", leaf_size=30)
        y_ = knn.fit(x, y).predict(T)
        error = np.sqrt(mean_squared_error(T, y_))
        
        if get_error:
            return T, y_, error
        else:
            return T, y_
        

    def xy4simplePCAplot(self, func_approx, params_vari, params_ranges,
                         param_x = False, steps = 100):
        '''
        Returns X, Y for 1D plot given fitted chaospy fct, varied parameters, 
        parameter ranges, the parameter of interest and optional the number of steps
        
        Parameters
        ----------
        - func_approx: chaospy approximated fit function
        - params_vari: list of varied parameters
        - params_ranges: dictionary of parameter ranges          
        - param_x: first parameter of interest
        - steps: default to 100 (increase for more accuracy in plot)
        
        Returns
        ----------
        X, Y for 1D plot
        '''
        x_range = np.linspace(params_ranges[param_x][0],params_ranges[param_x][1],steps)
        params_mean = self.paramMeanUniformDistribution(params_vari, params_ranges)[0]
        
        argslist = []
        for pars in params_vari:
            if pars==param_x:
                argslist.append(x_range)
            else:
                argslist.append(params_mean[pars])
        
        y_range = func_approx(*argslist)
        return x_range, y_range

    def meshgrid4PCEplot(self, func_approx, params_vari, params_ranges,
                         param_x = False, param_y = False, steps = 100,
                         cutplane = False):
        '''
        Returns X, Y, Z of a meshgrid for 2D colorcoded plot given fitted chaospy fct,
        varied parameters, parameter ranges, the two parameters of interest and
        optional the number of steps
        
        Parameters
        ----------
        - params_vari: list of varied parameters
        - params_ranges: dictionary of parameter ranges  
        - func_approx: chaospy approximated fit function
        - param_x: first parameter of interest
        - param_y: second parameter of interest
        - steps: default to 100 (increase for more accuracy in plot)
        
        Returns
        ----------
        X, Y, Z as meshgrid for colorcoded 2D plot
        '''
        if (not param_x) or (not param_y):
            param_x = params_vari[0]
            param_y = params_vari[1]
            print("Warning: given paramters not specified, taking {} and {}"
                  .format(param_x, param_y))
        elif (not param_x in params_vari) or (not param_y in params_vari):
            print("{} or {} not in parameter list, aborting".format(param_x, param_y))
            import warnings
            warnings.warn("""\nSome parameter value is not matching the parameter list""")
        
        params_mean = self.paramMeanUniformDistribution(params_vari, params_ranges)[0]

        # checks if a cutting plane of the other parameters is given
        if type(cutplane) == dict:
            for param in cutplane.keys():
                params_mean[param] = cutplane[param]
        elif (cutplane and type(cutplane) != dict):
            print("""A cutting plane argument different from the mean value is given, 
            but can't be evaluated\n\t\tplease use a dictionary!""")
        
        x1 = np.linspace(params_ranges[param_x][0],params_ranges[param_x][1],steps)
        x2 = np.linspace(params_ranges[param_y][0],params_ranges[param_y][1],steps)
        X, Y = np.meshgrid(x1,x2)
        
        argslist = []
        for pars in params_vari:
            if pars==param_x:
                argslist.append(X)
            elif pars==param_y:
                argslist.append(Y)
            else:
                argslist.append(params_mean[pars])
        
        Z = func_approx(*argslist)
        return X, Y, Z


    def splitBioDataFramesforXvalidation(self, sim_df, train_ratio=False,
                                           untr_ratio=False, offset=1):
        """
        internal fct splitting the dataframe from biomarker analysis into two chunks
        vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
        NOTE: USE sklearn.model_selection.train_test_split instead OF THIS FCT!
        ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        Args:
        - ratios: 0.0-1.0 or 0%-100% ratios of trained and untrained data
        - offset: the offset where to start picking untrained data
        
        Returns 
        - two df: the 'untrained' and 'trained' DF
        """
        # check for the format of the ratio given
        if train_ratio and (1 <= float(train_ratio) <= 100):
            train_ratio = float(train_ratio)/100
            untr_ratio = 1 - train_ratio
        elif untr_ratio and (1 <= float(untr_ratio) <= 100):
            untr_ratio = float(untr_ratio)/100
            train_ratio = 1 - train_ratio
        elif train_ratio and (0 <= float(train_ratio) <= 1):
            train_ratio = float(train_ratio)
            untr_ratio = 1 - train_ratio
        elif untr_ratio and (0 <= float(untr_ratio) <= 1):
            untr_ratio = float(untr_ratio)
            train_ratio = 1 - train_ratio
        else:
            print("Given ratios are not valid!")
            #return sim_df
        
        totlines = len(sim_df)
        untrlines = int(totlines*untr_ratio)
        #trainedlines = totlines - untrlines
        if offset > (totlines - untrlines):
            offset = totlines - untrlines
        
        untrDF = sim_df[offset: offset+untrlines]
        trDF = sim_df.drop(sim_df.index[offset: offset + untrlines])
        return untrDF, trDF

    def reduceClefts_byRatioThreshold(self, df, threshold = 0.5, reduct = "ratioRyR"):
        '''
        Returns DF with clefts having open ratio > threshold
        
        Parameters
        ----------
        - pandas DataFrame with ratio of RyR, LCC and cleft names
        - threshold (optional, minimum ratio to consider cleft)
        - reduct (ratio has to be considered) std: ratioRyR
        
        Returns
        ----------
        pandas DataFrame with relevant clefts
        '''
        rel_clefts = df[df[reduct] > threshold]["clefts"].drop_duplicates()
        res = pd.DataFrame()
        for cleft in rel_clefts:
            res = pd.concat([res, df[df["clefts"] == cleft]])
        return res

    def plotBiomarkerGrid(self, dataset, params, objectives, saveFig = False):
        """
        plot matrix with biomarkers given dataset with biomarkers, parameters and objectives
        
        Parameters
        ----------
        dataset: DF with biomarkers
        params: list of parameters vs which to plot
        objectives: y-values of plotting matrix
        saveFig: (optional) path where to save plot        
        """
        
        depth = len(objectives)
        width = len(params)
        
        colors = ['red','darkorange', 'gold','royalblue','seagreen', 'black', 'purple']
        plt.rcParams.update({'font.size': 12})
        plt.rcParams.update({'xtick.labelsize': 12})
        plt.rcParams.update({'ytick.labelsize': 12})
        plt.rcParams.update({'axes.labelsize': 12})
        
        fig, axarr = plt.subplots(depth, width, sharey='row', figsize=(3.5*width, 3*depth))
        
        for ind, par in enumerate(params):
            for jnd, obj in enumerate(objectives):
                a = min(dataset[par])
                b = max(dataset[par])
                axarr[jnd,ind].set_xlim((a,b))
                axarr[jnd,ind].locator_params(nbins=5, axis='x')
                axarr[jnd,ind].scatter(dataset[par], dataset[obj], 
                                       marker='o',color=colors[ind],s=100)
                axarr[jnd,ind].set_xlabel(par)
                axarr[jnd,ind].set_ylabel(obj)        
                
        fig.tight_layout()
        if saveFig == True:
            print("Attention: no saving path specified, saving into cwd")
            fig.savefig("./biomarkers_{}_vs_{}.pdf".format(objectives, params))
        elif type(saveFig) == str:
            fig.savefig(saveFig)
        else:
            fig.show()

    def fastOneSimPlot(self, simnum, simfolder, first = "Ca_i", second = "Vm", third = "I_LCC"):
        '''
        Plots three whole time series given the simulation folder and parent folder
        
        Parameters
        ----------
        - sim parent folder
        - first variable to be plotted (default: Ca_i)
        - second variable to be plotted (default: Vm)
        - third variable to be plotted (default: I_LCC)  
        '''
        
        fig1, axs = plt.subplots(3, sharex = True, figsize = (20, 11))
        from . import extract
        f = extract.extract(simnum, simfolder)
        ion = f.readIonic()
        
        axs[0].plot(ion["time"], ion[first], linewidth = 2)
        axs[1].plot(ion["time"], ion[second], linewidth = 2)    
        axs[2].plot(ion["time"], ion[third], linewidth = 2)
        # axis labels
        axs[0].set_ylabel(r"${}$".format(first))
        axs[1].set_ylabel(r"${}$".format(second))
        axs[2].set_ylabel(r"${}$".format(third))
        axs[2].set_xlabel(r"t $[ms]$")
        # -
        axs[2].tick_params(axis = 'both', which = 'major')
        
        fig1.set_tight_layout(True)
        
        fig1.show()
        #fig1.saveplot("figures-remote/three_plots{}".format(simnum))
