######################################################################################################################################
#
#                FireSimulator FBP serial and parallel version 1.0 - April 2018
#                Authors: Cristobal Pais, David L. Woodruff
#                example: mpiexec -n X python Path\Simulator1Beta.py  where X is the number of parallel processes
#
######################################################################################################################################

"""
March 2018, DLW. We now reqire a weather file with a header and at least
one row. The header (for now) is somewhat "FBP-friendly".
"""

# Importations
import inspect
import sys
import os
from itertools import repeat
import itertools
import numpy as np
import pandas as pd

# Weather class: Weather objects where all parameters like Wind speed, Wind direction, 
# Dew Point, Temperature, etc., are computed and updated from an initial state
"""
March 2018: A weather data frame is required. It might have only one row.
The new class has a public interface that callers can send in a date-time
or a row number and get all the corresponding weather-related data one
way or the other. (E.g., ffmc might be from the weather file or might
be computed). The weather data can be returned in the form of a named 
attributes with FBP input field names (extras are OK, "Scenario"). Here
is what we are hoping for, but we can live without some...
because they may come in from Data.dat...:
datetime, AACP, TMP, RH, WS, WD, FFMC, DMC, DC, ISI, BUI, FWI
TBD: if the weather data does not have FFMC, compute it.
"""
class Weather:
    
    def __init__ (self,wdf_path):
        """Constructor:
            args: wdf is a path to a csv file with rows over time
            outputs: updates the object so that access to the weather can be 
                     entirely via the object.
        """
        self._wdf = pd.read_csv(wdf_path, sep=",", index_col="datetime") # header implied
        self._wdf.columns = self._wdf.columns.str.lower()
        

    def set_columns(self, df, datetime=None, Row=None):
        """
        Update df (presumably for FBP) using whatever data we have from the weather 
        inputs:
            df: data frame to update
        """
        assert (datetime is None or Row is None)
        if Row != None:
            for dfcolname in df.columns:
                if dfcolname in self._wdf.columns:
                    weatherval = self._wdf.loc[self._wdf.index[Row], dfcolname]
                    df.loc[:,dfcolname] = weatherval
        else:
            raise RuntimeError("Datetime not supported yet")
                        
    
    ##### CP: TO BE MODIFIED!
    def update_Weather_FBP(self,df, WeatherOpt, weatherperiod=None, datetime=None):
        #Updates the current weather in df 
        # Has some code for random weather Weather
        if WeatherOpt == "constant":
            print("weather constant, so why is",inspect.stack()[0][3],"called from",inspect.stack()[1][3])
            return df

        if WeatherOpt != "random":
            #print "Weather is not random\n"
            self.set_columns(df, Row=weatherperiod)
            return df

        else:
            print("WARNING: DLW: random weather needs maintenance")
            print("         CP:  random sampling from file? distribution form hist. data?")
            row = 0            
            for i in df['ws']:
                #print "Row",row," I",i
                WS = np.round(np.random.uniform(-i, 30),2)

                if (i + WS) <= 50:
                    df.set_value(row, 'ws', i+WS) 
                else:
                    df.set_value(row, 'ws', 50.0)         
                row+=1
            #print "WS:",df.loc[:,"ws"]

            row = 0
            for i in df['waz']:
                #print "Row",row," I",i
                WAZ = np.round(np.random.uniform(-15, 15),2)

                if (i + WAZ) <= 359 and (i + WAZ) >= 0:
                    df.set_value(row, 'waz', i+WAZ) 
                elif (i + WAZ) >= 360:
                    df.set_value(row, 'waz', i+WAZ-360)
                elif (i + WAZ) < 0:
                    df.set_value(row, 'waz', i+WAZ+360)
                row+=1
            #print "WAZ:",df.loc[:,"waz"]
            # dlw says: do not change the slope
        
            return df
        """
        def update_Weather(self,period,optweather):
        #Updates the weather in the DF
        # Weather's random coefficients
        # DLW feb 2018: note that file based updates seem to take place in main and/or a different routine here
            
        if optweather == "random":
            WS = round(uniform(-0.15, 0.15),2)
            WD = round(uniform(-0.30, 0.30),2)
            TP = round(uniform(-0.10, 0.10),2)
            DP = round(uniform(-0.05, 0.05),2)
            RA = round(uniform(-0.03, 0.03),2)
            RD = round(uniform(-0.01, 0.01),2)
            
            #Rain probability: pr normal, qr if rain before
            pr = 0.10   
            qr = 0.05
            RAmount = 20
            
            #Update values
            self.WindSpeed = round(self.WindSpeed*(1+WS),2)
            if self.WindDirection*(1+WD) > 360:
                self.WindDirection = round(self.WindDirection*(1+WD)-360,2)
            else: 
                self.WindDirection = round(self.WindDirection*(1+WD),2)
            
            self.Temperature = round(self.Temperature*(1+TP),2)
            self.DPoint = round(self.DPoint*(1+DP),2)
            self.Radiation = round(self.Radiation*(1+RD),2)
            
            prain = round(uniform(0, 1),2)
            if self.Rain == 0:
                if prain < pr:
                    self.Rain = RAmount
                else:
                    self.Rain = 0
            else:
                if prain < qr:
                    self.Rain = round(self.Rain*(1+RA),2)
                else:
                    self.Rain = 0
        """
        
    
    # Prints-out weather report for current period (TBD)
    def print_info(self, period):
        print("Weather Info for weather period ", str(period))
        print("Report TBD xxxxx (need to get the row from the dataframe)")
    

