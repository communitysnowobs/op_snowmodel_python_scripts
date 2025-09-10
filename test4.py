#!/usr/bin/env python
# coding: utf-8

import geopandas as gpd
import pandas as pd
import numpy as np
from datetime import datetime, timedelta, date
import requests
import json
from rasterstats import point_query
from shapely import geometry as sgeom
import ulmo
from collections import OrderedDict

df = pd.DataFrame(
    {
        "City": ["Chicago", "Portland"],
        "Country": ["USA", "USA"],
        "Latitude": [40, 45],
        "Longitude": [-110, -120],
    }
)

gdf = gpd.GeoDataFrame(
    df, geometry=gpd.points_from_xy(df.Longitude, df.Latitude), crs="EPSG:4326"
)

print(gdf)

TD = np.array([point_query([val], '/nfs/attic/dfh/data/depth2swe/td_final.txt')[0] for val in gdf.geometry])
print(TD)