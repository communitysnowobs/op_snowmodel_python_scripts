#!/usr/bin/env python
# coding: utf-8
#this script is a 'one-time' script. It will simply pull the locations of snotel sites
#and then these should get written to a shapefile...

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


#########################################################################
# SNOTEL Functions
#########################################################################

# functions to get SNOTEL stations as geodataframe
def sites_asgdf(ulmo_getsites, stn_proj):
    """ Convert ulmo.cuahsi.wof.get_sites response into a point GeoDataframe
    """
    
    # Note: Found one SNOTEL site that was missing the location key
    sites_df = pd.DataFrame.from_records([
        OrderedDict(code=s['code'], 
        longitude=float(s['location']['longitude']), 
        latitude=float(s['location']['latitude']), 
        name=s['name'], 
        elevation_m=s['elevation_m'])
        for _,s in ulmo_getsites.items()
        if 'location' in s
    ])

    sites_gdf = gpd.GeoDataFrame(
        sites_df, 
        geometry=gpd.points_from_xy(sites_df['longitude'], sites_df['latitude']),
        crs=stn_proj
    )
    return sites_gdf
    
#function to get snotel stations.    
def get_snotel_stns():
    
    #USA
    Bbox = {"latmax": 72, "latmin": 18, "lonmin": -172, "lonmax": -66}
    #Bbox = {"latmax": 44, "latmin": 43.5, "lonmin": -127, "lonmax": -122}
    stn_proj = "epsg:4326"

    # Convert the bounding box dictionary to a shapely Polygon geometry using sgeom.box
    box_sgeom = sgeom.box(Bbox['lonmin'], Bbox['latmin'], Bbox['lonmax'], Bbox['latmax'])
    box_gdf = gpd.GeoDataFrame(geometry=[box_sgeom], crs=stn_proj)
    
    # WaterML/WOF WSDL endpoint url 
    wsdlurl = "https://hydroportal.cuahsi.org/Snotel/cuahsi_1_1.asmx?WSDL"
    
    # get dictionary of snotel sites 
    sites = ulmo.cuahsi.wof.get_sites(wsdlurl,user_cache=True)

    #turn sites to geodataframe 
    snotel_gdf = sites_asgdf(sites,stn_proj)
    
    #clip snotel sites to domain bounding box
    gdf = gpd.sjoin(snotel_gdf, box_gdf, how="inner")
    gdf.drop(columns='index_right', inplace=True)
    gdf.reset_index(drop=True, inplace=True)

    return gdf
    
#call the function to get snotel stations
snotel_gdf = get_snotel_stns()
#print('Retrieved station locations')
#print(snotel_gdf)

#keep only the geometry so that we can plot the points
cols=['geometry']
snotel_gdf_subset=snotel_gdf[cols]
#print(snotel_gdf_subset)   

#define name of output file. NOTE: I have decided to not create individual files for
#every day. Instead, we will have a single file, and simply append to this each day
path_to_file=OUTpath+'snotel_locations.shp'
path=Path(path_to_file) 

#next, we check to see if this file exists. The first time you call data, it won't.
#if it does exist, just exit
if path.is_file():
	print('File exists; quitting now')	
#if it does not exist use regular model
else:
	print('File does not yet exist; write data to new file')
	snotel_gdf_subset.to_file(path_to_file)