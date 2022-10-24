#!/usr/bin/env python
# coding: utf-8

# In[1]:


#https://github.com/giswqs/geemap/blob/master/examples/notebooks/11_export_image.ipynb
import ee
import geemap
import numpy as np
import matplotlib.pyplot as plt
import os
#from paths import *
import requests
import pandas as pd
import xarray as xr
from os import listdir
from datetime import datetime, timedelta, date
import contextlib

# Initialize the Earth Engine module
ee.Initialize()


# In[9]:


#########################################################################
############################ USER INPUTS ################################
#########################################################################
# PATHS
# path to temporary folder to store tif files from gee
TIFpath = '/nfs/depot/cce_u1/hill/dfh/op_snowmodel/get_met_data/GEE_Downloads_wa_sq_backfill/'
# path to where you want your output met .dat fime
OUTpath = '/nfs/depot/cce_u1/hill/dfh/op_snowmodel/wa_sq_snowmodel/met/mm_wa_sq.dat'

# DOMAIN
# choose the modeling domain
domain = 'WA_SQ'

# TIME
# choose if want to set 'manual' or 'auto' date 
date_flag = 'manual'
# If you choose 'manual' set your dates below  
# This will start on the 'begin' date at 0:00 and the last iteration will 
# be on the day before the 'end' date below.
st_dt = '2022-10-15'
ed_dt = '2022-10-21'
#########################################################################


# In[12]:


# Date setup function
def set_dates(st_dt,ed_dt,date_flag):
    if date_flag == 'auto':
        # ###automatically select date based on today's date 
        hoy = date.today()
        antes = timedelta(days = 2)
        #end date 3 days before today's date
        fecha = hoy - antes
        eddt = fecha.strftime("%Y-%m-%d") 
        #whole water year
        if (hoy.month == 10) & (hoy.day == 3):
            eddt = fecha.strftime("%Y-%m-%d") 
            stdt = str(hoy.year - 1)+'-10-01'
        #start dates
        elif fecha.month <10:
            stdt = str(fecha.year - 1)+'-10-01'
        else:
            stdt = str(fecha.year)+'-10-01'
    elif date_flag == 'manual':
        stdt = st_dt
        # add one day to end date because GEE ends on date before last date
        eddt = (datetime.strptime(ed_dt, "%Y-%m-%d")+timedelta(days = 1)).strftime("%Y-%m-%d") 
    return stdt, eddt


# Download CFSv2 met data function
def get_cfsv2(domain, TIFpath, stdt, eddt):
    # in GEE the last iteration is on the day before the 'end' date below
    # we adjust this here since it is not intuative
    #eddt = (datetime.strptime(eddt, '%Y-%m-%d')+timedelta(days = 1)).strftime('%Y-%m-%d')
    
    #create directory with initiation date for ensemble if it doesn't exist
    get_ipython().system('mkdir -p $TIFpath')

    #path to CSO domains
    domains_resp = requests.get("https://raw.githubusercontent.com/snowmodel-tools/preprocess_python/master/CSO_domains.json")
    domains = domains_resp.json()

    '''
    // These are the min and max corners of your domain in Lat, Long
    // Western Wyoming
    // Input the minimum lat, lower left corner
    '''
    minLat = domains[domain]['Bbox']['latmin']
    #// Input the minimum long, lower left corner
    minLong = domains[domain]['Bbox']['lonmin']
    #// Input the max lat, upper right corner
    maxLat = domains[domain]['Bbox']['latmax']
    #// Input the max Long, upper right corner
    maxLong = domains[domain]['Bbox']['lonmax']

    #/ These are the min and max corners of your reanalysis in Lat, Long (create a slightly larger box)
    #// Input the minimum lat, lower left corner
    minLatMET = (minLat - 0.25);
    #// print(minLat2);
    #// Input the minimum long, lower left corner
    minLongMET = (minLong - 0.5);
    #// Input the max lat, upper right corner
    maxLatMET = (maxLat + 0.25);
    #// Input the max Long, upper right corner
    maxLongMET = (maxLong + 0.5);

    # This resolution for the NLCD and DEM outputs for the SnowModel domain in meters
    sm_resolution = 100

    '''// Resolution for the PRISM output. This shoud change by Latitude of the domain
    // because the PRISM product spatial resolution is 2.5 minutes, which equals 150 arc seconds.
    // You can use this arc-second calculator to estimate the correct value for the PRISM resolution by latitude
    // https://opendem.info/arc2meters.html
    // This is one arc-second in meters for 43 degrees N Latitude'''
    one_arcsecond = 22.57
    PRISM_resolution = one_arcsecond * 150

    '''// Define the final output projection using EPSG codes'''
    epsg_code = domains[domain]['mod_proj']

    #// Name the DEM output
    dem_name = 'DEM'
    #// Name the Land Cover output
    lc_name = 'NLCD2016'

    my_domain = ee.Geometry.Rectangle(**{'coords':[minLong,minLat,maxLong,maxLat],'proj': 'EPSG:4326','geodesic':True,});
    my_domain_met = ee.Geometry.Rectangle([minLongMET,minLatMET,maxLongMET,maxLatMET])
    
    # download reanalysis data
    cfsv2 = ee.ImageCollection('NOAA/CFSV2/FOR6H')         .filterBounds(my_domain_met)         .filter(ee.Filter.date(stdt,eddt))

    data = cfsv2.select('Temperature_height_above_ground',         'Geopotential_height_surface',         'u-component_of_wind_height_above_ground',         'v-component_of_wind_height_above_ground',         'Pressure_surface',         'Specific_humidity_height_above_ground',         'Precipitation_rate_surface_6_Hour_Average')
    with contextlib.redirect_stdout(None):
        geemap.ee_export_image_collection(data, out_dir=TIFpath,region=my_domain_met,scale=22200,crs=epsg_code)

# In[21]:


# function to check for missing dates
def missing_slice_check(stdt, eddt, TIFpath):
    # create a 6-hourly timeseries with no missing values from the start to end date
    timesin = pd.date_range(start=stdt, end=eddt, freq='6H')[:-1]
    for time in timesin:
        nam = time.strftime('%Y%m%d%H')
        
    # compile list of all tif files downloaded from gee
    gee_times =[]

    for file in listdir(TIFpath):
        if file.endswith("tif"):
            datetmp = datetime.strptime(file[:-4], '%Y%m%d%H')
            gee_times.append(datetmp)
    gee_times = sorted(gee_times)

    # check for to see if all time slices downloaded from GEE
    if len(timesin) != len(gee_times):
    #### on 4/16 Nina edited code to print all missing timeslices
        print('gee is missing timeslice(s):\n',timesin[~timesin.isin(gee_times)].values)
        
        # if 4 or more consecutive timeslices are missing - quit the function
        duration = []
        for i in range(len(gee_times)-1):
            time_delta = gee_times[i+1] - gee_times[i]
            duration.append(time_delta.total_seconds()/60/60)
        if max(duration) >= 48:    
            print('at least two full days of met data are missing - quitting function')

        # if there are less than 4 missing consecutive time slices 
        # repeat the met data from the last valid time slice 
        else:	
            missing_idx = np.where(~timesin.isin(gee_times))	
            missing_dt = timesin[missing_idx]	
            if len(missing_dt)==1:	
                if missing_idx == 0:	
                    pre_dt=TIFpath+timesin[np.squeeze(missing_idx)+1].strftime('%Y%m%d%H')+'.tif'	
                    mis_dt = TIFpath+timesin[np.squeeze(missing_idx)].strftime('%Y%m%d%H')+'.tif' 	
                    get_ipython().system('cp $pre_dt $mis_dt')	
                    print('replaced', timesin[np.squeeze(missing_idx)].strftime('%Y%m%d%H'),' with ', timesin[np.squeeze(missing_idx)-1].strftime('%Y%m%d%H'))	
                else:	
                    pre_dt=TIFpath+timesin[np.squeeze(missing_idx)-1].strftime('%Y%m%d%H')+'.tif'	
                    mis_dt = TIFpath+timesin[np.squeeze(missing_idx)].strftime('%Y%m%d%H')+'.tif' 	
                    get_ipython().system('cp $pre_dt $mis_dt')	
                    print('replaced', timesin[np.squeeze(missing_idx)].strftime('%Y%m%d%H'),' with ', timesin[np.squeeze(missing_idx)-1].strftime('%Y%m%d%H'))	
            else:	
                for j in range(len(missing_dt)):	
                    if np.squeeze(missing_idx)[j] == 0:	
                        print('choose earlier start date so missing time slices can be filled in')	
                    else:	
                        pre_dt=TIFpath+timesin[np.squeeze(missing_idx)[j]-1].strftime('%Y%m%d%H')+'.tif'	
                        mis_dt = TIFpath+timesin[np.squeeze(missing_idx)[j]].strftime('%Y%m%d%H')+'.tif' 	
                        get_ipython().system('cp $pre_dt $mis_dt')	
                        print('replaced', timesin[np.squeeze(missing_idx)[j]].strftime('%Y%m%d%H'),' with ', timesin[np.squeeze(missing_idx)[j]-1].strftime('%Y%m%d%H'))	



# set time parameters
stdt, eddt = set_dates(st_dt,ed_dt,date_flag)


# In[17]:


# download GEE data
get_cfsv2(domain, TIFpath, stdt, eddt)


# In[22]:


# fill in missing time slices or throw error if missing >4 slices
missing_slice_check(stdt, eddt, TIFpath)

# delete directory with tif files 
get_ipython().system('rm -rf $TIFpath')

