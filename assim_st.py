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


#########################################################################
############################ USER INPUTS ################################
#########################################################################
# NOTE: to runn assim, set irun_data_assim = 1 in .par file


# DOMAIN
# choose the modeling domain
domain = 'ST'

# PATHS
dataPath = '/nfs/attic/dfh/Aragon2/CSOdmn/'+domain+'/'
#path to dem .tif
dem_path = dataPath + 'DEM_'+domain+'.tif'
#path to landcover .tif
lc_path = dataPath + 'NLCD2016_'+domain+'.tif'
#path to SnowModel
SMpath = '/nfs/depot/cce_u1/hill/dfh/op_snowmodel/st_snowmodel/'

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
        antes = timedelta(days = 2)
        #end date 3 days before today's date
        fecha = hoy - antes
        eddt = fecha.strftime("%Y-%m-%d")  
        #whole water year
        if (hoy.month == 10) & (hoy.day == 2):
            eddt = fecha.strftime("%Y-%m-%d")
            stdt = str(hoy.year - 1)+'-10-01'
        #start dates
        elif fecha.month <10:
            stdt = str(fecha.year - 1)+'-10-01'
        else:
            stdt = str(fecha.year)+'-10-01'
    elif date_flag == 'manual':
        stdt = st_dt
        eddt = ed_dt 
    return stdt, eddt

stdt, eddt = set_dates(st_dt,ed_dt,date_flag)
print(stdt, eddt)
#########################################################################
# CSO Functions
#########################################################################
# Function to get SWE from CSO Hs

def swe_calc(gdf):
    #convert snow depth to mm to input into density function
    H = gdf.depth.values*10
    #Get temp info at each point
    TD = np.array([point_query([val], '/nfs/attic/dfh/data/depth2swe/td_final.txt')[0] for val in gdf.geometry])
    #Get pr info at each point
    PPTWT = np.array([point_query([val], '/nfs/attic/dfh/data/depth2swe/ppt_wt_final.txt')[0] for val in gdf.geometry])
    #Determine day of year
    dates = pd.to_datetime(gdf.timestamp, format='%Y-%m-%dT%H:%M:%S').dt.date.values
    DOY = [date.toordinal(date(dts.year,dts.month,dts.day))-date.toordinal(date(dts.year,9,30)) for dts in dates]
    DOY = np.array([doy + 365 if doy < 0 else doy for doy in DOY])
    #Apply regression equation 
    a = [0.0533,0.948,0.1701,-0.1314,0.2922] #accumulation phase
    b = [0.0481,1.0395,0.1699,-0.0461,0.1804]; #ablation phase
    SWE = a[0]*H**a[1]*PPTWT**a[2]*TD**a[3]*DOY**a[4]*(-np.tanh(.01*\
            (DOY-180))+1)/2 + b[0]*H**b[1]*PPTWT**b[2]*TD**b[3]*DOY**b[4]*\
            (np.tanh(.01*(DOY-180))+1)/2;
    #convert swe to m to input into SM
    gdf['swe'] = SWE/1000
    gdf['doy'] = DOY
    gdf['H'] = H
    return gdf


# Function to build geodataframe of CSO point observations 

def get_cso(st, ed, domain):
    ''' 
    st = start date 'yyyy-mm-dd'
    ed = end date 'yyyy-mm-dd'
    domain = string label of defined CSO domain
    '''
    
    #path to CSO domains
    domains_resp = requests.get("https://raw.githubusercontent.com/snowmodel-tools/preprocess_python/master/CSO_domains.json")
    domains = domains_resp.json()

    Bbox = domains[domain]['Bbox']
    stn_proj = domains[domain]['stn_proj']
    mod_proj = domains[domain]['mod_proj']
    
    #Issue CSO API observations request and load the results into a GeoDataFrame
    params = {
      "bbox": f"{Bbox['lonmin']},{Bbox['latmax']},{Bbox['lonmax']},{Bbox['latmin']}",
      "start_date": st,
      "end_date": ed,
      "format": "geojson",
      "limit": 5000,
    }

    csodata_resp = requests.get("https://api.communitysnowobs.org/observations", params=params)
    csodatajson = csodata_resp.json()
    #turn into geodataframe
    gdf = gpd.GeoDataFrame.from_features(csodatajson, crs=stn_proj)
    
    mask = (gdf['timestamp'] >= st) & (gdf['timestamp'] <= ed)
    gdf = gdf.loc[mask]
    gdf=gdf.reset_index(drop=True)
    print('Total number of CSO in domain = ',len(gdf))

    #ingdf = extract_meta(gdf,domain,dem_path,lc_path)
    ingdf = swe_calc(gdf)
    ingdf_proj = ingdf.to_crs(mod_proj)
    
    ingdf['dt'] = pd.to_datetime(ingdf['timestamp'], format='%Y-%m-%dT%H:%M:%S').dt.date
    ingdf['Y'] = pd.DatetimeIndex(ingdf['dt']).year
    ingdf['M'] = pd.DatetimeIndex(ingdf['dt']).month
    ingdf['D'] = pd.DatetimeIndex(ingdf['dt']).day
    ingdf["x"] = ingdf_proj.geometry.x
    ingdf["y"] = ingdf_proj.geometry.y
    
    return ingdf


# QA/QC function for CSO data
def qaqc_iqr(csodf):
    print('Performing qa/qc on CSO data using IQR method')
    clim_dir = '/nfs/attic/dfh/data/snodas/snodas_tif/clim/'
    iqr_flag = []
    for i in range(len(csodf)):
        # get cso snow depth
        csohs = csodf.H[i]
        # get date
        dates = pd.to_datetime(csodf.timestamp[i], format='%Y-%m-%dT%H:%M:%S')

        # define path names for 1st and 3rd doy quantiles
        q1_Fname = clim_dir+dates.strftime("%m")+dates.strftime("%d")+'1036q1.tif'
        q3_Fname = clim_dir+dates.strftime("%m")+dates.strftime("%d")+'1036q3.tif'

        q1 = point_query([csodf.geometry[i]], q1_Fname)[0]
        q3 = point_query([csodf.geometry[i]], q3_Fname)[0]
        IQR = q3-q1
        # False = outlier
        iqr_flag.append((csohs > (q1-1.5*IQR)) & (csohs < (q3+1.5*IQR)))

    csodf['iqr_flag'] = iqr_flag
    csodf_clean = csodf.loc[csodf['iqr_flag'] == True]
    csodf_clean = csodf_clean.reset_index(drop=True)
    return csodf_clean

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

def get_snotel_stns(domain):
    
    #path to CSO domains
    domains_resp = requests.get("https://raw.githubusercontent.com/snowmodel-tools/preprocess_python/master/CSO_domains.json")
    domains = domains_resp.json()

    #Snotel bounding box
    Bbox = domains[domain]['Bbox']

    # Snotel projection
    stn_proj = domains[domain]['stn_proj']
    # model projection
    mod_proj = domains[domain]['mod_proj']

    # Convert the bounding box dictionary to a shapely Polygon geometry using sgeom.box
    box_sgeom = sgeom.box(Bbox['lonmin'], Bbox['latmin'], Bbox['lonmax'], Bbox['latmax'])
    box_gdf = gpd.GeoDataFrame(geometry=[box_sgeom], crs=stn_proj)
    
    # WaterML/WOF WSDL endpoint url 
    if domain == 'NH':
        wsdlurl = "https://hydroportal.cuahsi.org/Scan/cuahsi_1_1.asmx?WSDL"
    else:
        wsdlurl = "https://hydroportal.cuahsi.org/Snotel/cuahsi_1_1.asmx?WSDL"

    # get dictionary of snotel sites 
    sites = ulmo.cuahsi.wof.get_sites(wsdlurl,user_cache=True)

    #turn sites to geodataframe 
    snotel_gdf = sites_asgdf(sites,stn_proj)
    
    #clip snotel sites to domain bounding box
    gdf = gpd.sjoin(snotel_gdf, box_gdf, how="inner")
    gdf.drop(columns='index_right', inplace=True)
    gdf.reset_index(drop=True, inplace=True)

    #add columns with projected coordinates 
    CSO_proj = gdf.to_crs(mod_proj)
    gdf['easting'] = CSO_proj.geometry.x
    gdf['northing'] = CSO_proj.geometry.y
    
    return gdf


def fetch(sitecode, variablecode, domain,start_date, end_date):
    print(sitecode, variablecode, domain,start_date, end_date)
    values_df = None
    # WaterML/WOF WSDL endpoint url 
    if domain == 'NH':
        wsdlurl = "https://hydroportal.cuahsi.org/Scan/cuahsi_1_1.asmx?WSDL"
        network = 'SCAN:'
    else:
        wsdlurl = "https://hydroportal.cuahsi.org/Snotel/cuahsi_1_1.asmx?WSDL"
        network = 'SNOTEL:'

    try:
        #Request data from the server
        site_values = ulmo.cuahsi.wof.get_values(
            wsdlurl, network+sitecode, variablecode, start=start_date, end=end_date
        )
        #Convert to a Pandas DataFrame   
        values_df = pd.DataFrame.from_dict(site_values['values'])
        #Parse the datetime values to Pandas Timestamp objects
        values_df['datetime'] = pd.to_datetime(values_df['datetime'])
        #Set the DataFrame index to the Timestamps
        values_df.set_index('datetime', inplace=True)
        #Convert values to float and replace -9999 nodata values with NaN
        values_df['value'] = pd.to_numeric(values_df['value']).replace(-9999, np.nan)
        #Remove any records flagged with lower quality
        values_df = values_df[values_df['quality_control_level_code'] == '1']
    except:
        print("Unable to fetch %s" % variablecode)
    
    return values_df


# returns daily timeseries of snotel variables 
# https://www.wcc.nrcs.usda.gov/web_service/AWDB_Web_Service_Reference.htm#commonlyUsedElementCodes
# 'WTEQ': swe [in]
# 'SNWD': snow depth [in]
# 'PRCP': precipitation increment [in]
# 'PREC': precipitation accumulation [in]
# 'TAVG': average air temp [F]
# 'TMIN': minimum air temp [F]
# 'TMAX': maximum air temp [F]
# 'TOBS': observered air temp [F]
def get_snotel_data(gdf,sddt, eddt,var,domain,units='metric'):
    '''
    gdf - pandas geodataframe of SNOTEL sites
    st_dt - start date string 'yyyy-mm-dd'
    ed_dt - end date string 'yyyy-mm-dd'
    var - snotel variable of interest 
    units - 'metric' (default) or 'imperial'
    '''
    stn_data = pd.DataFrame(index=pd.date_range(start=stdt, end=eddt))
    if domain == 'NH':
        network = 'SCAN:'
    else:
        network = 'SNOTEL:'    

    for sitecode in gdf.code:
        try:
            data = fetch(sitecode,network+var+'_D', domain, start_date=stdt, end_date=eddt)
            #check for nan values
            if len(data.value[np.isnan(data.value)]) > 0:
                #check if more than 10% of data is missing
                if len(data.value[np.isnan(data.value)])/len(data) > .02:
                    print('More than 2% of days missing')
                    gdf.drop(gdf.loc[gdf['code']==sitecode].index, inplace=True)
                    continue
            stn_data[sitecode] = data.value
        except:
            gdf.drop(gdf.loc[gdf['code']==sitecode].index, inplace=True)     
    
    gdf.reset_index(drop=True, inplace=True)
    if units == 'metric':
        if (var == 'WTEQ') |(var == 'SNWD') |(var == 'PRCP') |(var == 'PREC'):
            #convert SNOTEL units[in] to [m]
            for sitecode in gdf.code:
                stn_data[sitecode] = 0.0254 * stn_data[sitecode]
        elif (var == 'TAVG') |(var == 'TMIN') |(var == 'TMAX') |(var == 'TOBS'):
            #convert SNOTEL units[F] to [C]
            for sitecode in gdf.code:
                stn_data[sitecode] = (stn_data[sitecode] - 32) * 5/9
    return gdf, stn_data



#########################################################################
# Functions to format CSO & SNOTEL data for SM 
#########################################################################

def make_SMassim_file(CSOdata,outFpath):
    '''
    CSOdata = dataframe with CSO data 
    
    outFpath = output path to formated assim data for SM 
    '''
    print('Generating assim file')
    f= open(outFpath,"w+")

    tot_obs=len(CSOdata)
    uq_day = np.unique(CSOdata.dt)
    num_days = len(uq_day)
    f.write('{:02.0f}\n'.format(num_days))
    for j in range(len(uq_day)):
        obs = CSOdata[CSOdata['dt']==uq_day[j]]
        d=CSOdata.D[CSOdata['dt']==uq_day[j]].values
        m=CSOdata.M[CSOdata['dt']==uq_day[j]].values
        y=CSOdata.Y[CSOdata['dt']==uq_day[j]].values
        date = str(y[0])+' '+str(m[0])+' '+str(d[0])
        obs_count = str(len(obs))
        f.write(date+' \n')
        f.write(obs_count+' \n')
        for k in range(len(obs)):
            ids = 100+k
            x= obs.x[obs.index[k]]
            y=obs.y[obs.index[k]]
            swe=obs.swe[obs.index[k]]
            f.write('{:3.0f}\t'.format(ids)+'{:10.0f}\t'.format(x)+'{:10.0f}\t'.format(y)+'{:3.2f}\n'.format(swe))
    f.close() 
    

def make_SMassim_file_snotel(STswe,STmeta,outFpath):
    '''
    STmeta = dataframe with SNOTEL sites
    
    STswe = dataframe with SWE data 
    
    outFpath = output path to formated assim data for SM 
    '''
    print('Generating assim file')
    f= open(outFpath,"w+")

    tot_obs=np.shape(STswe)[0]*np.shape(STswe)[1]
    uq_day = np.shape(STswe)[0]
    stn = list(STswe.columns)
    f.write('{:02.0f}\n'.format(uq_day))
    for j in range(uq_day):
        d=STswe.index[j].day
        m=STswe.index[j].month
        y=STswe.index[j].year
        date = str(y)+' '+str(m)+' '+str(d)
        stn_count = np.shape(STswe)[1]
        f.write(date+' \n')
        f.write(str(stn_count)+' \n')
        ids = 100
        for k in stn:
            ids = ids + 1 
            x = STmeta.easting.values[STmeta.code.values == k][0]
            y = STmeta.northing.values[STmeta.code.values == k][0]
            swe = STswe[k][j]
            f.write('{:3.0f}\t'.format(ids)+'{:10.0f}\t'.format(x)+'{:10.0f}\t'.format(y)+'{:3.2f}\n'.format(swe))
    f.close() 


def make_SMassim_file_both(STswe,STmeta,CSOdata,outFpath):
    '''
    STmeta = dataframe with SNOTEL sites
    
    STswe = dataframe with SWE data 
    
    CSOdata = dataframe with CSO data
    
    outFpath = output path to formated assim data for SM 
    '''
    print('Generating assim file')

    f= open(outFpath,"w+")
    
    #determine number of days with observations to assimilate
    if STswe.shape[1]>0:
        uq_day = np.unique(np.concatenate((STswe.index.date,CSOdata.dt.values)))
        f.write('{:02.0f}\n'.format(len(uq_day)))
    else:
        uq_day = np.unique(CSOdata.dt.values)
        f.write('{:02.0f}\n'.format(len(uq_day)))
    
    # determine snotel stations 
    stn = list(STswe.columns)
    
    # ids for CSO observations - outside of loop so each observation is unique
    IDS = 500
    
    #add assimilation observations to output file
    for i in range(len(uq_day)):

        SThoy = STswe[STswe.index.date == uq_day[i]]
        CSOhoy = CSOdata[CSOdata.dt.values == uq_day[i]]

        d=uq_day[i].day
        m=uq_day[i].month
        y=uq_day[i].year

        date = str(y)+' '+str(m)+' '+str(d)

        if len(SThoy)>0:
            stn_count = len(stn) + len(CSOhoy)
        else:
            stn_count = len(CSOhoy)
        
        if stn_count > 0:
            f.write(date+' \n')
            f.write(str(stn_count)+' \n')

        #go through snotel stations for that day 
        ids = 100
        if len(SThoy) > 0:
            for k in stn:
                ids = ids + 1 
                x = STmeta.easting.values[STmeta.code.values == k][0]
                y = STmeta.northing.values[STmeta.code.values == k][0]
                swe = SThoy[k].values[0]
                f.write('{:3.0f}\t'.format(ids)+'{:10.0f}\t'.format(x)+'{:10.0f}\t'.format(y)+'{:3.2f}\n'.format(swe))    
        #go through cso obs for that day 
        if len(CSOhoy) > 0:
            for c in range(len(CSOhoy)):
                IDS = IDS + 1 
                x= CSOhoy.x[CSOhoy.index[c]]
                y=CSOhoy.y[CSOhoy.index[c]]
                swe=CSOhoy.swe[CSOhoy.index[c]]
                f.write('{:3.0f}\t'.format(IDS)+'{:10.0f}\t'.format(x)+'{:10.0f}\t'.format(y)+'{:3.2f}\n'.format(swe))
    f.close()
    return len(uq_day)

    
    
    
#########################################################################
# Functions to edit SM files 
#########################################################################
# function to edit SnowModel Files other than .par
# for assim - have to adjust .inc file to specify # of obs being assimilated
def replace_line(file_name, line_num, text):
    ''' 
    file_name = file to edit
    line_num = line number in file to edit
    text = nex text to put in
    '''
    lines = open(file_name, 'r').readlines()
    lines[line_num] = text
    out = open(file_name, 'w')
    out.writelines(lines)
    out.close()


#edit par file for correct number of timesteps 
parFile = SMpath + 'snowmodel.par'
replace_line(parFile,7,str(datetime.strptime(stdt,'%Y-%m-%d').year) +'                          !iyear_init - start year\n')
replace_line(parFile,8,str(datetime.strptime(stdt,'%Y-%m-%d').month) +'                        !imonth_init - start month\n')
replace_line(parFile,9,str(datetime.strptime(stdt,'%Y-%m-%d').day) +'                            !iday_init - start day\n')
value = str((datetime.strptime(eddt,'%Y-%m-%d')-datetime.strptime(stdt,'%Y-%m-%d')).days*4+4)
replace_line(parFile,11,value +'			!max_iter - number of model time steps\n')


# Run SM with CSO assim

outFpath = SMpath+'swe_assim/swe_obs_test.dat'
codepath = SMpath+'/code/'
incFile = SMpath+'code/snowmodel.inc'

if assim_mod == 'none':
    print('Executing SnowModel without assimilation')
    replace_line(parFile,35,'0			!irun_data_assim - 0 for straight run; 1 for assim run\n')      
    #compile SM
    get_ipython().run_line_magic('cd', '$codepath')
    get_ipython().system(' ./compile_snowmodel.script')
    #run snowmodel
    get_ipython().run_line_magic('cd', '$SMpath')
    get_ipython().system(' ./snowmodel')
    
elif assim_mod == 'cso':
    CSOgdf = get_cso(stdt, eddt, domain)
    if len(CSOgdf) < 1:
        print('Executing SnowModel without assimilation')
        replace_line(parFile,35,'0			!irun_data_assim - 0 for straight run; 1 for assim run\n')      
    else:    
        print('Creating assim input file using CSO observations')
        replace_line(parFile,35,'1			!irun_data_assim - 0 for straight run; 1 for assim run\n')
        CSOgdf_clean = qaqc_iqr(CSOgdf)
        make_SMassim_file(CSOgdf_clean,outFpath)
#     #edit .inc file
    replace_line(incFile, 30, '      parameter (max_obs_dates='+str(len(CSOgdf_clean)+1)+')\n')
    #compile SM
    get_ipython().run_line_magic('cd', '$codepath')
    get_ipython().system(' ./compile_snowmodel.script')
    #run snowmodel
    get_ipython().run_line_magic('cd', '$SMpath')
    get_ipython().system(' ./snowmodel')

elif assim_mod == 'snotel':
    print('Creating assim input file using SNOTEL observations')
    replace_line(parFile,35,'1			!irun_data_assim - 0 for straight run; 1 for assim run\n')
    snotel_gdf = get_snotel_stns(domain)
    SNOTELgdf, swe = get_snotel_data(snotel_gdf,stdt,eddt,'WTEQ',domain)
    delta = 5
    sample = swe.iloc[::delta,:]
    make_SMassim_file_snotel(sample,SNOTELgdf,outFpath)
    #edit .inc file
    replace_line(incFile, 30, '      parameter (max_obs_dates='+str(len(sample)+1)+')\n')
    #compile SM
    get_ipython().run_line_magic('cd', '$codepath')
    get_ipython().system(' ./compile_snowmodel.script')
     #run snowmodel
    get_ipython().run_line_magic('cd', '$SMpath')
    get_ipython().system(' ./snowmodel')

elif assim_mod == 'both':
    print('Creating assim input file using CSO & SNOTEL observations')
    replace_line(parFile,35,'1			!irun_data_assim - 0 for straight run; 1 for assim run\n')
    CSOgdf = get_cso(stdt, eddt, domain)
    # set delta time 
    delta = 5
    if len(CSOgdf)>=1:
        CSOgdf_clean = qaqc_iqr(CSOgdf)  
        CSOdata = CSOgdf_clean.sort_values(by='dt',ascending=True)
        CSOdata = CSOdata.reset_index(drop=True)
        newCSO = CSOdata
        snotel_gdf = get_snotel_stns(domain)
        SNOTELgdf, STswe = get_snotel_data(snotel_gdf,stdt,eddt,'WTEQ',domain)
        newST = SNOTELgdf
        newSTswe = STswe
        num_obs = make_SMassim_file_both(newSTswe,newST,newCSO,outFpath)
        #edit .inc file
        replace_line(incFile, 30, '      parameter (max_obs_dates='+str(num_obs+1)+')\n')
    else:
        print('No CSO observations. Creating assim input file using SNOTEL observations')
        snotel_gdf = get_snotel_stns(domain)
        SNOTELgdf, STswe = get_snotel_data(snotel_gdf,stdt,eddt,'WTEQ',domain)
        newST = SNOTELgdf
        newSTswe = STswe
        make_SMassim_file_snotel(newSTswe,newST,outFpath)
         #edit .inc file
        replace_line(incFile, 30, '      parameter (max_obs_dates='+str(len(newSTswe)+1)+')\n')
#     #compile SM
    get_ipython().run_line_magic('cd', '$codepath')
    get_ipython().system(' ./compile_snowmodel.script')
    #run snowmodel
    get_ipython().run_line_magic('cd', '$SMpath')
    get_ipython().system(' ./snowmodel')

else:
    print("No valid assim mode was entered. Select 'none','cso', 'snotel' or 'both'.")