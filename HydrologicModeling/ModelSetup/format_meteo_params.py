import os
import netCDF4
import numpy as np
import pandas as pd
from osgeo import gdal
from osgeo.gdalnumeric import *  
from osgeo.gdalconst import *
from datetime import datetime

t1 = datetime.now()

def find_nearest_idx(xx,yy,xval,yval):    
    xidx = (np.abs(xx-xval)).argmin()
    yidx = (np.abs(yy-yval)).argmin()
    
    ridx = yidx / xx.shape[1]
    cidx = xidx % xx.shape[1]
        
    return [ridx,cidx]
    
band = 1

infiles = [r'~\gis\GCSfiles\LMB_mask_GCS_01deg_full.tif']

data = np.zeros((138,98, len(infiles)))

for i in range(len(infiles)):
    
    ds = gdal.Open(infiles[i],GA_ReadOnly)
    b1 = ds.GetRasterBand(band)
    var = BandReadAsArray(b1)
    
    if i == 0:
        gt = ds.GetGeoTransform()
        lon0 = gt[0]
        lon1 = gt[0] + (98*0.1)
        lat0 = gt[3] - (138*0.1)
        lat1 = gt[3]
    
    data[:,:,i] = var[:,:]
        
    ds = None
    b1 = None
    
lons = np.linspace(lon0,lon1,98)
lats = np.linspace(lat0,lat1,138)

xx,yy = np.meshgrid(lons,lats)

yy = np.flipud(yy)

years = np.arange(1981,2011)

for yr in years:
    #for climate data
    
    if yr < 2006:
    
        ncdfs = [r'~\pr_historic_bcc-csm1-1_{0}.nc'.format(yr),
                r'~\tasmax\tasmax_historic_bcc-csm1-1_{0}.nc'.format(yr),
                r'~\tasmin\tasmin_historic_bcc-csm1-1_{0}.nc'.format(yr),
                r'~\models\input\forcing\ERAdata\WIND_E_{0}.nc'.format(yr),
                r'~\models\input\forcing\ERAdata\WIND_N_{0}.nc'.format(yr)]
    else:      
        ncdfs = [r'~\pr\pr_rcp45_bcc-csm1-1_{0}.nc'.format(yr),
                r'~\tasmax\tasmax_rcp45_bcc-csm1-1_{0}.nc'.format(yr),
                r'~\tasmin\tasmin_rcp45_bcc-csm1-1_{0}.nc'.format(yr),
                r'~\models\input\forcing\ERAdata\WIND_E_{0}.nc'.format(yr),
                r'~\models\input\forcing\ERAdata\WIND_N_{0}.nc'.format(yr)]
    
    #for ERA/CHIRPS data
    #ncdfs = [r'D:/Kel/UAH/classes/ESS/thesis/models/input/forcing/chirps_climo/chirps-v2.0.{0}.days_p05.nc'.format(yr),
    #        r'~\models\input\forcing\ERAdata\TMAX_{0}.nc'.format(yr),
    #        r'~\models\input\forcing\ERAdata\TMIN_{0}.nc'.format(yr),
    #        r'~\models\input\forcing\ERAdata\WIND_E_{0}.nc'.format(yr),
    #        r'~\models\input\forcing\ERAdata\WIND_N_{0}.nc'.format(yr)]
    
    #for climate data    
    prnc = netCDF4.Dataset(ncdfs[0])
    prv = prnc.variables['pr']
    platnc = prnc.variables['lat'];plat=platnc[:]
    plonnc = prnc.variables['lon'];plon=plonnc[:]
    
    tmaxnc = netCDF4.Dataset(ncdfs[1])
    tmaxv = tmaxnc.variables['tasmax']
    
    tminnc = netCDF4.Dataset(ncdfs[2])
    tminv = tminnc.variables['tasmin']
    
    #for ERA/CHIRPS data
    #prnc = netCDF4.Dataset(ncdfs[0])
    #prv = prnc.variables['precip']
    #platnc = prnc.variables['latitude'];plat=platnc[:]
    #plonnc = prnc.variables['longitude'];plon=plonnc[:]
    #
    #tmaxnc = netCDF4.Dataset(ncdfs[1])
    #tmaxv = tmaxnc.variables['MX2T_GDS4_SFC']
    #
    #tminnc = netCDF4.Dataset(ncdfs[2])
    #tminv = tminnc.variables['MN2T_GDS4_SFC']
    
    uwinnc = netCDF4.Dataset(ncdfs[3])
    uwinv = uwinnc.variables['10U_GDS4_SFC']
    
    vwinnc = netCDF4.Dataset(ncdfs[4])
    vwinv = vwinnc.variables['10V_GDS4_SFC']
    latncv = vwinnc.variables['g4_lat_1'];latnc=latncv[:]
    lonncv = vwinnc.variables['g4_lon_2'];lonnc=lonncv[:]    
    
    lonncs,latncs = np.meshgrid(lonnc,latnc)
    plons,plats = np.meshgrid(plon,plat)
    
    mask = data[:,:,0].astype(uint8)
    mask = np.ma.masked_where(mask!=1,mask)
    
    if yr%4 == 0:
        days = 366
    else:
        days = 365
        
    print "Processing year: {0}".format(yr)
    
    for i in range(yy.shape[0]):
        for j in range(yy.shape[1]):
            x = xx[i,j]
            y = yy[i,j]
            
            if mask.mask[i,j] == False: 
            
                #print x,y
                
                tidx = find_nearest_idx(lonncs,latncs,x,y)
                pidx = find_nearest_idx(plons,plats,x,y)
    
                cnt = 1
    
                meteofile = 'D:\\Kel\\UAH\\classes\\ESS\\thesis\\models\\input\\forcing\\historic_clim\\forcing_{0:.4f}_{1:.4f}'.format(y,x)
                #if os.path.exists(meteofile) == False:
                with open(meteofile, 'a') as f:
                    for t in range(days):
                        
                        range0 = t*4
                        range1 = range0+4
                        
                        range2 = t*2
                        range3 = range2+2
                            
                        #for cliamte data
                        if t < 365:
                            td = t
                        else:
                            td = t-1
                        
                        prval = prv[td,pidx[0],pidx[1]]*86400
                        tmaxval = tmaxv[td,pidx[0],pidx[1]]-273.15
                        tminval = tminv[td,pidx[0],pidx[1]]-273.15
                        
                        ##for ERA/CHIRPS data
                        #prval = np.percentile(prv[t,pidx[0]-1:pidx[0]+1,pidx[1]-1:pidx[1]+1],75)
                        #tmaxval = tmaxv[range2:range3,0,tidx[0],tidx[1]].max()-273.15
                        #tminval = tminv[range2:range3,0,tidx[0],tidx[1]].min()-273.15
                        
                        uwndval = uwinv[range0:range1,tidx[0],tidx[1]].mean()
                        vwndval = vwinv[range0:range1,tidx[0],tidx[1]].mean()
                        windval = np.sqrt(uwndval**2 + vwndval**2)
                        
                        if prval < 0.0:
                            prval = np.max(prv[t,pidx[0]-1:pidx[0]+1,pidx[1]-1:pidx[1]+1])
                                            
                        
                        #if (t > 365) & (t < 730):
                        #    calday = datetime.strptime(str(t+1-365),'%j')
                        #    mon = calday.month
                        #    windval = data[i,j,mon]
                        #elif t >= 730:
                        #    calday = datetime.strptime(str(t+1-730),'%j')
                        #    mon = calday.month
                        #    windval = data[i,j,mon]
                        #    
                        #else:
                        #    calday = datetime.strptime(str(t+1),'%j')
                        #    mon = calday.month
                        #    windval = data[i,j,mon]
            
                        f.write('{0} {1} {2} {3} {4} {5}\n'.format(prval,
                                tmaxval,tminval,windval,vwndval,uwndval))
                        
    prnc.close()
    tmaxnc.close()
    tminnc.close()
    vwinnc.close()
    uwinnc.close()

                    
dt = datetime.now() - t1

print "COMPLETE\nProcessing time: {0}".format(dt)
