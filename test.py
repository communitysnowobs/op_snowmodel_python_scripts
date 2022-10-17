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

# TIME
# choose if want to set 'manual' or 'auto' date 
date_flag = 'auto'
# If you choose 'manual' set your dates below  
st_dt = '2019-01-01'
ed_dt = '2019-01-05'

# ASSIM OPTIONS
# select the data source to be assimilated
# can be set to 'none','cso', 'both' or 'snotel'
assim_mod = 'both'
print(assim_mod)
#########################################################################

# Date setup function
def set_dates(st_dt,ed_dt,date_flag):
    if date_flag == 'auto':
        # ###automatically select date based on today's date 
        hoy = date.today()
        antes = timedelta(days = 3)
        #end date 3 days before today's date
        fecha = hoy - antes
        eddt = fecha.strftime("%Y-%m-%d")  
        #whole water year
        if (hoy.month == 10) & (hoy.day == 3):
            eddt = fecha.strftime("%Y-%m-%d")
            stdt = '2021-10-01'
        #start dates
        elif fecha.month <10:
            stdt = '2021-10-01'
        else:
            stdt = '2021-10-01'
    elif date_flag == 'manual':
        stdt = st_dt
        eddt = ed_dt 
    return stdt, eddt

stdt, eddt = set_dates(st_dt,ed_dt,date_flag)
print(stdt, eddt)

# Function to build geodataframe of CSO point observations 

def get_cso(st, ed):
    ''' 
    st = start date 'yyyy-mm-dd'
    ed = end date 'yyyy-mm-dd'
    domain = string label of defined CSO domain
    '''
    
    #path to CSO domains
    domains_resp = requests.get("https://raw.githubusercontent.com/snowmodel-tools/preprocess_python/master/CSO_domains.json")
    domains = domains_resp.json()

    Bbox = {"latmax": 90, "latmin": -90, "lonmin": -180, "lonmax": 180}
    print(Bbox)
    stn_proj = "epsg:4326"
    print(stn_proj)
    
    #Issue CSO API observations request and load the results into a GeoDataFrame
    params = {
      "bbox": f"{Bbox['lonmin']},{Bbox['latmax']},{Bbox['lonmax']},{Bbox['latmin']}",
      "start_date": st,
      "end_date": ed,
      "format": "geojson",
      "limit": 10,
    }

    csodata_resp = requests.get("https://api.communitysnowobs.org/observations", params=params)
    csodatajson = csodata_resp.json()
    #turn into geodataframe
    gdf = gpd.GeoDataFrame.from_features(csodatajson, crs=stn_proj)
    
    #mask = (gdf['timestamp'] >= st) & (gdf['timestamp'] <= ed)
    #gdf = gdf.loc[mask]
    print(gdf)
    gdf=gdf.reset_index(drop=True)
    print('Total number of CSO in domain = ',len(gdf))
    print(gdf)
    ingdf=gdf

    #ingdf = extract_meta(gdf,domain,dem_path,lc_path)
    #ingdf = swe_calc(gdf)
    #ingdf_proj = ingdf.to_crs(mod_proj)
    
    ingdf['dt'] = pd.to_datetime(ingdf['timestamp'], format='%Y-%m-%dT%H:%M:%S').dt.date
    ingdf['Y'] = pd.DatetimeIndex(ingdf['dt']).year
    ingdf['M'] = pd.DatetimeIndex(ingdf['dt']).month
    ingdf['D'] = pd.DatetimeIndex(ingdf['dt']).day
    ingdf["longitude"] = ingdf.geometry.x
    ingdf["latitude"] = ingdf.geometry.y
    
    return ingdf
    
def df_to_geojson(df, properties, lat='latitude', lon='longitude'):
    # create a new python dict to contain our geojson data, using geojson format
    geojson = {'type':'FeatureCollection', 'features':[]}

    # loop through each row in the dataframe and convert each row to geojson format
    for _, row in df.iterrows():
        # create a feature template to fill in
        feature = {'type':'Feature',
                   'properties':{},
                   'geometry':{'type':'Point',
                               'coordinates':[]}}

        # fill in the coordinates
        feature['geometry']['coordinates'] = [row[lon],row[lat]]

        # for each column, get the value and add it as a new feature property
        for prop in properties:
            feature['properties'][prop] = row[prop]
        
        # add this feature (aka, converted dataframe row) to the list of features inside our dict
        geojson['features'].append(feature)
    
    return geojson    
    
    
    
CSOgdf = get_cso(stdt, eddt)
cols=['geometry','depth','Y','M','D','longitude','latitude']
#cols=['geometry','depth','Y','M','D']
CSOgdf_subset=CSOgdf[cols]
print(CSOgdf_subset)  
cols2=['depth','Y','M','D']
geojson=df_to_geojson(CSOgdf_subset,cols2) 
#geojson=json.dumps(geojson)
print(geojson)
#geojson.to_file('test.json',driver="GeoJSON")
with open('test.json','w') as outfile:
	dump(geojson, outfile)



