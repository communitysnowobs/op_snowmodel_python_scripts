#!/usr/bin/env python
# coding: utf-8
#this python script will pull CSO observations for desired days. It will create a shapefile
#and save the data to a file.

import geopandas as gpd
import pandas as pd
import numpy as np
from datetime import datetime, timedelta, date
import requests
import json
from geojson import Point, Feature, FeatureCollection, dump
from rasterstats import point_query
from shapely import geometry as sgeom
import ulmo
from collections import OrderedDict
from pathlib import Path

#########################################################################
############################ USER INPUTS ################################
#########################################################################

#path for created shape files
OUTpath = '/scratch/ms_shapefiles/'

# TIME - note that fetched data will range from 12 am on st_dt to 12 am
# on ed_dt
# for now, I am only using manual option for a test.  
#st_dt = '2022-10-01'
#ed_dt = '2022-10-05'

#stdt=st_dt
#eddt=ed_dt
#print(stdt, eddt)

# TIME
# choose if want to set 'manual' or 'auto' date 
date_flag = 'auto'
# If you choose 'manual' set your dates below  
# This will start on the 'begin' date at 0:00 and the last iteration will 
# be on the day before the 'end' date below.
st_dt = '2023-02-08'
ed_dt = '2023-02-15'
#########################################################################

# Date setup function
def set_dates(st_dt,ed_dt,date_flag):
    if date_flag == 'auto':
        # ###automatically select date based on today's date 
        hoy = date.today()
        antes = timedelta(days = 2)
        #start date 2 days before today's date
        fecha = hoy - antes
        stdt = fecha.strftime("%Y-%m-%d")
        antes = timedelta(days = 1)
        #end date 1 days before today's date
        fecha = hoy - antes
        eddt = fecha.strftime("%Y-%m-%d")
    elif date_flag == 'manual':
        stdt = st_dt
        eddt = ed_dt
    return stdt, eddt

# Function to build geodataframe of CSO point observations 
def get_cso(st, ed):
    ''' 
    st = start date 'yyyy-mm-dd'
    ed = end date 'yyyy-mm-dd'
    '''
    #path to CSO domains
    domains_resp = requests.get("https://raw.githubusercontent.com/snowmodel-tools/preprocess_python/master/CSO_domains.json")
    domains = domains_resp.json()

    Bbox = {"latmax": 90, "latmin": -90, "lonmin": -180, "lonmax": 180}
    stn_proj = "epsg:4326"
    
    #Issue CSO API observations request and load the results into a GeoDataFrame
    params = {
      "bbox": f"{Bbox['lonmin']},{Bbox['latmax']},{Bbox['lonmax']},{Bbox['latmin']}",
      "startDate": st,
      "endDate": ed,
      "format": "geojson",
      "limit": 10000,
    }

    csodata_resp = requests.get("https://api.communitysnowobs.org/observations", params=params)
    csodatajson = csodata_resp.json()
    
    if not 'features' in csodatajson or len(csodatajson['features']) == 0:
    	#no observations on this day, create empty dataframe
    	print('No CSO measurements')
    	gdf=pd.DataFrame()
    
    else:
    	#turn into geodataframe
    	gdf = gpd.GeoDataFrame.from_features(csodatajson, crs=stn_proj)
    	gdf['timestamp'] = pd.to_datetime(gdf.timestamp)
    	gdf.sort_values('timestamp', inplace=True)

    	print('Total number of CSO in domain = ',len(gdf))
    	gdf['dt'] = pd.to_datetime(gdf['timestamp'], format='%Y-%m-%dT%H:%M:%S').dt.date
    	gdf['Y'] = pd.DatetimeIndex(gdf['dt']).year
    	gdf['M'] = pd.DatetimeIndex(gdf['dt']).month
    	gdf['D'] = pd.DatetimeIndex(gdf['dt']).day
    
    return gdf  

#call the function to get dates
stdt, eddt = set_dates(st_dt,ed_dt,date_flag)
print(stdt, eddt)
#call the function to get cso data   
CSOgdf = get_cso(stdt, eddt)

#define name of output file. NOTE: I have decided to not create individual files for
#every day. Instead, we will have a single file, and simply append to this each day
path_to_file=OUTpath+'csodata.shp'
path=Path(path_to_file)

#below, we will only attempt to append data if the gdf is not empty (i.e., if new data exist)
if not CSOgdf.empty:
	#reduce it down to keep only desired variables
	cols=['geometry','depth','Y','M','D']
	CSOgdf_subset=CSOgdf[cols]
	print(CSOgdf_subset)
	print('Data were found; write to file')
	#next, we check to see if this file exists. The first time you call data, it won't.
	#if it does exist, use append mode
	if path.is_file():
		print('File exists; appending data')
		CSOgdf_subset.to_file(path_to_file, mode='a')	
	#if it does not exist use regular model
	else:
		print('File does not yet exist; write data to new file')
		CSOgdf_subset.to_file(path_to_file)
else:
	print('No measurements found')		
		
		
		