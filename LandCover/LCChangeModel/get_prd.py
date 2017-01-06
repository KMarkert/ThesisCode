import numpy as np
from osgeo import gdal
from osgeo.gdalnumeric import *
from osgeo.gdalconst import *
from demand_handler import *

def get_cls(cls_probs,old):    
    if old == 0:
        poss = np.array((0,1,2,4))
    elif old == 1:
        poss = np.array((1,2,4))
    elif old == 2:
        poss = np.array((1,2,4))
    elif old == 3:
        poss = np.array((3))
    elif old == 4:
        poss = np.array((4))
    elif old == 5:
        poss = np.array((1,5))
    else:
        poss = np.arange(0,6)
    
    val = cls_probs[poss].max()
    
    idx = np.where(cls_probs==val)
    
    if len(idx[0]) < 1:
        cls = 0
    else:
        cls = idx[0][0]
    
    return cls

def get_data(files):    
    band = 1
    
    data = np.zeros((1807,1534, len(files)))
    
    for i in range(len(files)):
        
        ds = gdal.Open(files[i],GA_ReadOnly)
        b1 = ds.GetRasterBand(band)
        var = BandReadAsArray(b1)
        
        data[:,:,i] = var[:,:]
        
    return data
    
def get_prob(clf,data,itervar):
    
    elas = np.array((0.5,0.8,0.7,1.,0.9,0.3))
    
    probs = np.zeros((data.shape[0],data.shape[1],6))
    
    for i in range(6):
        coef = clf.coef_[i]
        intercept = clf.intercept_[i]
    
        prob = ((data[:,:,0]*coef[0]) + (data[:,:,1]*coef[1]) + (data[:,:,2]*coef[2]) + 
                (data[:,:,3]*coef[3]) + (data[:,:,4]*coef[4]) + (data[:,:,5]*coef[5]) + 
                (data[:,:,6]*coef[6]) + (data[:,:,7]*coef[7]) + (data[:,:,8]*coef[8]) + 
                (data[:,:,9]*coef[9]) + (data[:,:,10]*coef[10]) + (data[:,:,11]*coef[11]) + 
                (data[:,:,12]*coef[12]) + (data[:,:,13]*coef[13]) + (data[:,:,14]*coef[14]) + intercept)
            
        probs[:,:,i] = (np.exp(prob)/(1+np.exp(prob))) + elas[i] + itervar[i]
    
    return probs
    
def get_prd(clf,files):
    
    #read in water mask data for water class
    fWater = r'~\gis\MODIS\WaterMask\MOD44W_WATER_2000_BB.tif'
    dsWM = gdal.Open(fWater, GA_ReadOnly)
    b1 = dsWM.GetRasterBand(1)
    WM = BandReadAsArray(b1)
    
    #read in Protected area mask
    fPrtA = r'~\gis\misc\WDPA_Sept2015-shapefile\PrtArea_raster.tif'
    dsPA = gdal.Open(fPrtA, GA_ReadOnly)
    b2 = dsPA.GetRasterBand(1)
    PA = BandReadAsArray(b2)
    
    dsWM = None
    b1 = None
    dsPA = None
    b2 = None
    
    data = get_data(files)
    
    ctrl = 0
    
    itervar = np.zeros(6)
    
    demand = get_demand(data[:,:,1])
    print demand
    
    while ctrl <= 10:
    
        probs = get_prob(clf,data,itervar)
        
        lcOut = np.zeros((data.shape[0],data.shape[1]))
        lcOut[:,:] = -1
            
        for i in range(probs.shape[0]):
            for j in range(probs.shape[1]):
                if PA[i,j] == 0:
                    lcOut[i,j] = data[i,j,1]
                else:
                    cls = get_cls(probs[i,j,:],data[i,j,1])
                    lcOut[i,j] = cls

                
        idx = np.where(WM==1)
        
        lcOut[idx] = 3
        
        ndidx = np.where(data[:,:,1]==15)
        lcOut[ndidx] = 15
        
        rchk = rate_check(lcOut,demand)
        
        print itervar
        
        if np.any(rchk == False) == True:
            itervar = update_itervar(lcOut,demand,itervar,rchk)
        else:
            ctrl = 1
        ctrl+=1
    
    return lcOut
