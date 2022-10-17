#!/usr/bin/env python
# coding: utf-8

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

#########################################################################
############################ USER INPUTS ################################
#########################################################################

#path for created shape files
OUTpath = '/scratch/ms_shapefiles/'

# TIME - note that fetched data will range from 12 am on st_dt to 12 am
# on ed_dt
# for now, I am only using manual option for a test.  
st_dt = '2022-01-01'
ed_dt = '2022-10-02'

stdt=st_dt
eddt=ed_dt

print(stdt, eddt)

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
      "limit": 5000,
    }

    csodata_resp = requests.get("https://api.communitysnowobs.org/observations", params=params)
    csodatajson = csodata_resp.json()
    
    #turn into geodataframe
    gdf = gpd.GeoDataFrame.from_features(csodatajson, crs=stn_proj)
    gdf['timestamp'] = pd.to_datetime(gdf.timestamp)
    gdf.sort_values('timestamp', inplace=True)

    print('Total number of CSO in domain = ',len(gdf))
    ingdf=gdf
    
    ingdf['dt'] = pd.to_datetime(ingdf['timestamp'], format='%Y-%m-%dT%H:%M:%S').dt.date
    ingdf['Y'] = pd.DatetimeIndex(ingdf['dt']).year
    ingdf['M'] = pd.DatetimeIndex(ingdf['dt']).month
    ingdf['D'] = pd.DatetimeIndex(ingdf['dt']).day
    #ingdf["longitude"] = ingdf.geometry.x
    #ingdf["latitude"] = ingdf.geometry.y
    
    return ingdf  
   
CSOgdf = get_cso(stdt, eddt)
print(CSOgdf)
cols=['geometry','depth','Y','M','D']
CSOgdf_subset=CSOgdf[cols]
print(CSOgdf_subset)  

CSOgdf_subset.to_file(OUTpath+stdt+'_'+eddt+'_cso.shp')	



