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

import matplotlib.pyplot as plt

class SparkAnalysis:
    """
    store analysis functions for the spark analysis
    """
    
    def __init__(self, folder):
        pass
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

# ----------------------------------------------------------------------------------------------

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

# ----------------------------------------------------------------------------------------------

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
    
    def APD(self, times, series, percent = 50.0, start_time = -np.infty, end_time = np.infty):
        """
        Compute the APD at (percent) % of a given series using the threshold_crossings fct
        given the times and the series. For instance, for percent = 90.0 the APD_90 values are computed.
        
        Parameters
        ----------
        - times
        - series
        - percent (default is 50%)
        
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

# ----------------------------------------------------------------------------------------------

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
        
        t_peaks, peaks = self.get_Ca_peaks(times, series, smooth = smooth, time_threshold = time_threshold, start_time=1000, end_time=6000)
        t_min, minima = self.get_Ca_sys_minima(times, series, smooth = smooth, time_threshold = time_threshold, start_time=1000, end_time=6000)
        
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

# ----------------------------------------------------------------------------------------------

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

    def fastOneSimPlot(self, simnum, simfolder, first = "Ca_i", second = "Vm", third = "I_LCC"):
        '''
        Plots three whole time series given the simulation folder and parent folder
        
        Parameters
        ----------
        - simfolder
        - sim parent folder
        - first variable to be plotted (default: Ca_i)
        - second variable to be plotted (default: Vm)
        - third variable to be plotted (default: I_LCC)
        
        Returns
        ----------
        Plot        
        '''
        
        fig1, axs = plt.subplots(3, sharex = True, figsize = (20, 11))
        from classes import extract
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
        
