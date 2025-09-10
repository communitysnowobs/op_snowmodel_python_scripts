#this script will pull CFS and GFS data and concatenate them

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
ee.Initialize(project='sunny-emissary-318920')

#########################################################################
############################ USER INPUTS ################################
#########################################################################
# PATHS
# path to temporary folder to store tif files from gee
TIFpath = '/nfs/depot/cce_u1/hill/dfh/op_snowmodel/get_met_data/GEE_Downloads_wa/'
# path to temporary folder to store tif files from gee gfs
TIFpath2 = '/nfs/depot/cce_u1/hill/dfh/op_snowmodel/get_met_data/GEE_Downloads_wa_gfs/'
# path to where you want your output met .dat fime
OUTpath = '/nfs/depot/cce_u1/hill/dfh/op_snowmodel/wa_ne_snowmodel/met/mm_wa.dat'

# DOMAIN
# choose the modeling domain
domain = 'WA'

# TIME
# choose if want to set 'manual' or 'auto' date 
date_flag = 'auto'
# If you choose 'manual' set your dates below  
# This will start on the 'begin' date at 0:00 and the last iteration will 
# be on the day before the 'end' date below.
st_dt = '2020-10-01'
ed_dt = '2020-10-15'
#########################################################################

#########################################################################
# Date setup function
# sddt --> October 1 of current water year
# eddt --> the day before 'today'. Note that in the GEE calls, eddt is 'exclusive'
# 	so if EDDT is 2023-02-12, CFS data get pulled running through 2023-02-11 (18Hr)
def set_dates(st_dt,ed_dt,date_flag):
    if date_flag == 'auto':
        # automatically select date based on today's date 
        hoy = date.today()
        # back up one day (to yesterday)
        antes = timedelta(days = 1)
        fecha = hoy - antes
        eddt = fecha.strftime("%Y-%m-%d") 
        
        #need to take some care for runs near the start of Oct
        if (hoy.month == 10) & (hoy.day <= 2):
            eddt = fecha.strftime("%Y-%m-%d") 
            stdt = str(hoy.year - 1)+'-10-01'
        #start dates
        elif fecha.month <10:
            stdt = str(fecha.year - 1)+'-10-01'
            #stdt='2023-02-20' #temp addition (Feb 16 2023). Delete when done testing
        else:
            stdt = str(fecha.year)+'-10-01'
            #stdt='2023-02-20' #temp addition (Feb 16 2023). Delete when done testing
    elif date_flag == 'manual':
        stdt = st_dt
        # add one day to end date because GEE ends on date before last date
        eddt = (datetime.strptime(ed_dt, "%Y-%m-%d")+timedelta(days = 1)).strftime("%Y-%m-%d") 
    return stdt, eddt
#########################################################################

#########################################################################
# Download CFSv2 met data function
def get_cfsv2(domain, TIFpath, stdt, eddt):

    #create directory with initiation date for ensemble if it doesn't exist
    get_ipython().system('mkdir -p $TIFpath')

    #path to CSO domains
    domains_resp = requests.get("https://raw.githubusercontent.com/snowmodel-tools/preprocess_python/master/CSO_domains.json")
    domains = domains_resp.json()

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

    #Define the final output projection using EPSG codes
    epsg_code = domains[domain]['mod_proj']

    my_domain = ee.Geometry.Rectangle(**{'coords':[minLong,minLat,maxLong,maxLat],'proj': 'EPSG:4326','geodesic':True,});
    my_domain_met = ee.Geometry.Rectangle([minLongMET,minLatMET,maxLongMET,maxLatMET])
    
    # download reanalysis data
    cfsv2 = ee.ImageCollection('NOAA/CFSV2/FOR6H')         .filterBounds(my_domain_met)         .filter(ee.Filter.date(stdt,eddt))
    data = cfsv2.select('Temperature_height_above_ground',         'Geopotential_height_surface',         'u-component_of_wind_height_above_ground',         'v-component_of_wind_height_above_ground',         'Pressure_surface',         'Specific_humidity_height_above_ground',         'Precipitation_rate_surface_6_Hour_Average')
    with contextlib.redirect_stdout(None):
        geemap.ee_export_image_collection(data, out_dir=TIFpath,region=my_domain_met,scale=22200,crs=epsg_code)  
#########################################################################

#########################################################################
# Download GFS met data function
def get_gfs(domain, TIFpath2, stdt2, eddt2):
    
    #create directory with initiation date for ensemble if it doesn't exist
    get_ipython().system('mkdir -p $TIFpath2')

    #path to CSO domains
    domains_resp = requests.get("https://raw.githubusercontent.com/snowmodel-tools/preprocess_python/master/CSO_domains.json")
    domains = domains_resp.json()

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

    #Define the final output projection using EPSG codes
    epsg_code = domains[domain]['mod_proj']

    my_domain = ee.Geometry.Rectangle(**{'coords':[minLong,minLat,maxLong,maxLat],'proj': 'EPSG:4326','geodesic':True,});
    my_domain_met = ee.Geometry.Rectangle([minLongMET,minLatMET,maxLongMET,maxLatMET])
    
    # download reanalysis data
    gfs = ee.ImageCollection('NOAA/GFS0P25').filterBounds(my_domain_met).filter(ee.Filter.date(stdt2,eddt2)).filter(ee.Filter.Or(ee.Filter.eq('forecast_hours',6), ee.Filter.eq('forecast_hours',12),ee.Filter.eq('forecast_hours',18),ee.Filter.eq('forecast_hours',24), ee.Filter.eq('forecast_hours',30),ee.Filter.eq('forecast_hours',36),ee.Filter.eq('forecast_hours',42),ee.Filter.eq('forecast_hours',48),ee.Filter.eq('forecast_hours',54),ee.Filter.eq('forecast_hours',60),ee.Filter.eq('forecast_hours',66),ee.Filter.eq('forecast_hours',72)))
    data2 = gfs.select('temperature_2m_above_ground', 'u_component_of_wind_10m_above_ground', 'v_component_of_wind_10m_above_ground',     'relative_humidity_2m_above_ground',         'total_precipitation_surface')
    with contextlib.redirect_stdout(None):
        geemap.ee_export_image_collection(data2, out_dir=TIFpath2,region=my_domain_met,scale=22200,crs=epsg_code) 
#########################################################################  

#########################################################################     
# function to check for missing dates (for the CFS data)
# cfs seems to drop time slices now and then...so far, with gfs forecast, I am not seeing
# that as a problem. So the check for missing dates is presently cfs only.
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
#########################################################################   

#########################################################################
# Format gee cfs files for SnowModel function
# this function will take care of the CFS data...there will be a second function that 
#will try to append the GFS data...Ha. Good luck.
def MET2SM(TIFpath, OUTpath, stdt, eddt):
    
    # create a 6-hourly timeseries with no missing values from the start to end date
    timesin = pd.date_range(start=stdt, end=eddt, freq='6H')[:-1]
    
    #load first tif to get dimensions
    ar = xr.open_rasterio(TIFpath+timesin[0].strftime('%Y%m%d%H')+'.tif')
    #print(TIFpath+timesin[0].strftime('%Y%m%d%H')+'.tif')
    
    # empty arrays for each met variable
    T = np.empty((len(timesin),ar.shape[1],ar.shape[2]))
    Z = np.empty((len(timesin),ar.shape[1],ar.shape[2]))
    U = np.empty((len(timesin),ar.shape[1],ar.shape[2]))
    V = np.empty((len(timesin),ar.shape[1],ar.shape[2]))
    P = np.empty((len(timesin),ar.shape[1],ar.shape[2]))
    H = np.empty((len(timesin),ar.shape[1],ar.shape[2]))
    PR = np.empty((len(timesin),ar.shape[1],ar.shape[2]))

    # extract met data from tifs 
    for i in range(len(timesin)):

        #load tif file
        nam = TIFpath+timesin[i].strftime('%Y%m%d%H')+'.tif'
        ar = xr.open_rasterio(nam)
        T[i,:,:] = ar[0,:,:]
        Z[i,:,:] = ar[1,:,:]
        U[i,:,:] = ar[2,:,:]
        V[i,:,:] = ar[3,:,:]
        P[i,:,:] = ar[4,:,:]
        H[i,:,:] = ar[5,:,:]
        PR[i,:,:] = ar[6,:,:]

    #number of timesteps per dat 
    pointsperday = 4

    #compute number of grid points and time steps from size of 3d matrix
    t,y,x=PR.shape
    gridpts=x*y
    tsteps=t

    #create y m d h vectors
    year = timesin.year
    month = timesin.month
    day = timesin.day
    hour = timesin.hour

    #create ID numbers for the grid points
    ID=1000000+np.linspace(1,gridpts,gridpts)

    #create matrices of x and y values
    X, Y = np.meshgrid(ar.x.values, ar.y.values)
    X=X.flatten(order='F')
    Y=Y.flatten(order='F')

    #elevation is static (doesn't change with time)
    elev=Z[1,:,:].flatten(order='F')

    #find number of grid points with <0 elevation. Note: this is related to the
    #subroutine met_data_check in the preprocess_code.f. that subroutine seems
    #to suggest that negative elevations are ok (say, death valley). But, the
    #code itself checks for negative elevations and stops execution is any
    #negatives are found.
    I = np.where(elev>=0)
    validgridpts=np.shape(I)[1]

    #remove data at points with neg elevations
    ID=ID[I]
    X=X[I]
    Y=Y[I]
    elev=elev[I]

    #we are now ready to begin our main loop over the time steps. The w+ below
    #overwrites file, if it exists.
    fid= open(OUTpath,"w+")

    for j in range(tsteps):
        #first we write the number of grid points
        fid.write('{0:6d}\n'.format(validgridpts))

        #prep data matrix for this time step. First, grab the jth time slice
        Prtmp=PR[j,:,:].flatten(order='F')
        Htmp=H[j,:,:].flatten(order='F')
        Ptmp=P[j,:,:].flatten(order='F')
        Ttmp=T[j,:,:].flatten(order='F')
        Utmp=U[j,:,:].flatten(order='F')
        Vtmp=V[j,:,:].flatten(order='F')

        #remove data at points with neg elevations
        Prtmp=Prtmp[I]
        Htmp=Htmp[I]
        Ptmp=Ptmp[I]
        Ttmp=Ttmp[I]
        Utmp=Utmp[I]
        Vtmp=Vtmp[I]

        #convert precip rate to precip DEPTH (mm) during time interval
        Prtmp=Prtmp*24*3600/pointsperday

        #convert specific hum. to RH from Clausius-Clapeyron. T is still in K
        RHtmp=0.263*Ptmp*Htmp*(np.exp(17.67*(Ttmp-273.16)/(Ttmp-29.65)))**(-1)

        #compute wind speed
        SPDtmp=np.sqrt(Utmp**2+Vtmp**2)

        #compute wind direction. 0-360, with 0 being true north! 90 east, etc.
        DIRtmp=np.degrees(np.arctan2(Utmp,Vtmp))
        K=np.where(DIRtmp>=180)
        J=np.where(DIRtmp<180)
        DIRtmp[K]=DIRtmp[K]-180
        DIRtmp[J]=DIRtmp[J]+180

        #put T in C
        Ttmp=Ttmp-273.16

        for z in range(len(Prtmp)): 

            fid.write('{:5.0f}\t'.format(int(year[j]))+'{:5.0f}\t'.format(int(month[j]))+
                      '{:3.0f}\t'.format(int(day[j]))+'{:6.3f}\t'.format(hour[j])+
                      '{:9.0f}\t'.format(int(ID[z]))+'{:12.1f}\t'.format(X[z])+
                      '{:12.1f}\t'.format(Y[z])+'{:8.1f}\t'.format(elev[z])+
                      '{:9.2f}\t'.format(Ttmp[z])+'{:9.2f}\t'.format(RHtmp[z])+
                      '{:9.2f}\t'.format(SPDtmp[z])+'{:9.2f}\t'.format(DIRtmp[z])+
                      '{:9.2f}\n'.format(Prtmp[z]))
    fid.close()  
#########################################################################

#########################################################################
# Format gee files for SnowModel function
# this function will take care of the GFS data.
def MET2SM2(TIFpath, TIFpath2, OUTpath, stdt, eddt):
    
    #lots of confusing timing stuff going on. Bottom line...the goal of this
    #is to grab three days of gfs data to tack on to the cfs data    
    eddt_tmp=(datetime.strptime(eddt,'%Y-%m-%d')+timedelta(days=3)).strftime('%Y-%m-%d')
    
    # create a 6-hourly timeseries with no missing values from the start to end date of the GFS data
    timesin = pd.date_range(start=eddt, end=eddt_tmp, freq='6H')[:-1]
    
    #dump listing of tifs (GFS) into array. Critical there are no extra files in there! 
    #there is probably a better way, but this seems the easiest for now...The sort is 
    #key to have the files in the right order.
    x = sorted(os.listdir(TIFpath2))
    
    #load first tif to get dimensions
    ar = xr.open_rasterio(TIFpath2+x[0])
    
    # empty arrays for each met variable
    T = np.empty((len(timesin),ar.shape[1],ar.shape[2]))
    U = np.empty((len(timesin),ar.shape[1],ar.shape[2]))
    V = np.empty((len(timesin),ar.shape[1],ar.shape[2]))
    H = np.empty((len(timesin),ar.shape[1],ar.shape[2]))
    PR = np.empty((len(timesin),ar.shape[1],ar.shape[2]))
    
    #note...Z data not in GFS, but GFS and CFS are same size. We will have to read Z from one of the CFS files, but 
    #go ahead and initialize size here.
    Z = np.empty((len(timesin),ar.shape[1],ar.shape[2]))

    # extract met data from tifs 
    for i in range(len(timesin)):

        #load tif file
        nam = TIFpath2+x[i]
        ar = xr.open_rasterio(nam)
        T[i,:,:] = ar[0,:,:]
        U[i,:,:] = ar[1,:,:]
        V[i,:,:] = ar[2,:,:]
        H[i,:,:] = ar[3,:,:]
        PR[i,:,:] = ar[4,:,:]

    #number of timesteps per dat 
    pointsperday = 4

    #compute number of grid points and time steps from size of 3d matrix
    t,y,x=PR.shape
    gridpts=x*y
    tsteps=t

    #create y m d h vectors
    year = timesin.year
    month = timesin.month
    day = timesin.day
    hour = timesin.hour

    #create ID numbers for the grid points
    ID=1000000+np.linspace(1,gridpts,gridpts)

    #create matrices of x and y values
    X, Y = np.meshgrid(ar.x.values, ar.y.values)
    X=X.flatten(order='F')
    Y=Y.flatten(order='F')
    
    #ok...so now let us deal with the Z data from one of the CFS files.
    # create a 6-hourly timeseries with no missing values from the start to end date
    timesin_tmp = pd.date_range(start=stdt, end=eddt, freq='6H')[:-1]
    
    #load first tif to get dimensions
    ar2 = xr.open_rasterio(TIFpath+timesin_tmp[0].strftime('%Y%m%d%H')+'.tif')
    print(TIFpath+timesin_tmp[0].strftime('%Y%m%d%H')+'.tif')
    Z[0,:,:] = ar2[1,:,:]

    #elevation is static (doesn't change with time)
    elev=Z[0,:,:].flatten(order='F')

    #find number of grid points with <0 elevation. Note: this is related to the
    #subroutine met_data_check in the preprocess_code.f. that subroutine seems
    #to suggest that negative elevations are ok (say, death valley). But, the
    #code itself checks for negative elevations and stops execution is any
    #negatives are found.
    I = np.where(elev>=0)
    validgridpts=np.shape(I)[1]

    #remove data at points with neg elevations
    ID=ID[I]
    X=X[I]
    Y=Y[I]
    elev=elev[I]

    #we are now ready to begin our main loop over the time steps. Note use of 'a' option
    #since we are appending to the snowmodel met file generated for the CFS data.
    fid= open(OUTpath,"a")

    for j in range(tsteps):
        #first we write the number of grid points
        fid.write('{0:6d}\n'.format(validgridpts))

        #prep data matrix for this time step. First, grab the jth time slice
        #note that precip is super confusing. The GEE data catalog documentation is pretty
        #poor. A better reference: https://benny.istan.to/blog/20210416-how-to-get-daily-rainfall-forecast-data-from-gfs-part-2
        #from that...it SEEMS, that if I just take the precip grids every 6 hours (which I have done),
        #then the grids will be precip depth over the previous 6 hour interval...which is precisely what
        #we want. It may take further investigation to confirm this.
        Prtmp=PR[j,:,:].flatten(order='F')
        Htmp=H[j,:,:].flatten(order='F')
        Ttmp=T[j,:,:].flatten(order='F')
        Utmp=U[j,:,:].flatten(order='F')
        Vtmp=V[j,:,:].flatten(order='F')

        #remove data at points with neg elevations
        Prtmp=Prtmp[I] #gfs precip is already depth...no need to convert
        Htmp=Htmp[I] #gfs humidity is relative already...no need to convert
        Ttmp=Ttmp[I] #gfs temp is C already, not K, no need to convert...
        Utmp=Utmp[I]
        Vtmp=Vtmp[I]

        #compute wind speed
        SPDtmp=np.sqrt(Utmp**2+Vtmp**2)

        #compute wind direction. 0-360, with 0 being true north! 90 east, etc.
        DIRtmp=np.degrees(np.arctan2(Utmp,Vtmp))
        K=np.where(DIRtmp>=180)
        J=np.where(DIRtmp<180)
        DIRtmp[K]=DIRtmp[K]-180
        DIRtmp[J]=DIRtmp[J]+180

        for z in range(len(Prtmp)): 

            fid.write('{:5.0f}\t'.format(int(year[j]))+'{:5.0f}\t'.format(int(month[j]))+
                      '{:3.0f}\t'.format(int(day[j]))+'{:6.3f}\t'.format(hour[j])+
                      '{:9.0f}\t'.format(int(ID[z]))+'{:12.1f}\t'.format(X[z])+
                      '{:12.1f}\t'.format(Y[z])+'{:8.1f}\t'.format(elev[z])+
                      '{:9.2f}\t'.format(Ttmp[z])+'{:9.2f}\t'.format(Htmp[z])+
                      '{:9.2f}\t'.format(SPDtmp[z])+'{:9.2f}\t'.format(DIRtmp[z])+
                      '{:9.2f}\n'.format(Prtmp[z]))
    fid.close()
#########################################################################

#########################################################################
# ok, let's run things.       

# set time parameters for CFS2
stdt, eddt = set_dates(st_dt,ed_dt,date_flag)
print(stdt)
print(eddt)     

# set time parameters for GFS. This will grab three additional days of forecast. The timing is tricky
# the cfs call grabs data that will end at hour 18 on the day before eddt. So, we wish to pick up with hour 0
# on eddt itself. BUT, the GFS forecast does not issue a complete forecast at hour 0. So, we need to
# go back to hour 18 of the day before eddt and grab GFS starting then. Whew. In addition to this, GFS forecasts
# are every hour, and we only want things every six hours, to match CFS.
stdt_tmp=(datetime.strptime(eddt,'%Y-%m-%d')+timedelta(days=-1)).strftime('%Y-%m-%d')
eddt_tmp=(datetime.strptime(eddt,'%Y-%m-%d')+timedelta(days=0)).strftime('%Y-%m-%d')
# append string to help narrow down 
stdt2=stdt_tmp + 'T18:00:00Z'
eddt2=eddt_tmp + 'T00:00:00Z' 
print(stdt2)
print(eddt2)  

# download GEE data for CFS
get_cfsv2(domain, TIFpath, stdt, eddt) 

# download GEE data for GFS
get_gfs(domain, TIFpath2, stdt2, eddt2) 

# fill in missing time slices or throw error if missing >4 slices. Currently, this only 
# operates on the CFS reanalysis and not the GFS forecast.
missing_slice_check(stdt, eddt, TIFpath)

# build SnowModel met file
MET2SM(TIFpath, OUTpath, stdt, eddt)

# build SnowModel met file
MET2SM2(TIFpath, TIFpath2, OUTpath, stdt, eddt)

# delete directory with tif files 
get_ipython().system('rm -rf $TIFpath')
get_ipython().system('rm -rf $TIFpath2')
#########################################################################
