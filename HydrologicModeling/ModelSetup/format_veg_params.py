import json
import numpy as np
from osgeo import gdal
from osgeo.gdalnumeric import *  
from osgeo.gdalconst import *

def format_veg_params(basinMask,lcData,,outVeg,scheme='IGBP'):
    
    if scheme == 'IGBP':
        attriFile = '~/core/veg_type_attributes_igbp.json'
    elif scheme == 'GLCC':
        attriFile = '~/core/veg_type_attributes_glcc.json'
    elif scheme == 'IPCC':
        attriFile = '~/core/veg_type_attributes_ipcc.json'
    else:
        raise NameError('Land cover classification scheme not supported')
    
    band = 1
    
    with open(attriFile) as data_file:    
        attriData = json.load(data_file)
    
    clsAttributes = attriData['classAttributes']

    infiles = [basinMask,lcData]
   
    ds = gdal.Open(infiles[0],GA_ReadOnly)
    b1 = ds.GetRasterBand(band)
    mask = BandReadAsArray(b1)
    maskRes = ds.GetGeoTransform()[1]
    ds = None
    b1 = None
    
    ds = gdal.Open(infiles[1],GA_ReadOnly)
    b1 = ds.GetRasterBand(band)
    lccls = BandReadAsArray(b1)  
    clsRes = ds.GetGeoTransform()[1]
    ds = None
    b1 = None
    
    ratio = maskRes/clsRes
    
    cells = 0
    
    vegfile = outVeg
    
    with open(vegfile, 'a') as f:
        
        cnt = 1
        
        for i in range(mask.shape[0]):
            y1 = i*ratio
            y2 = y1+ratio
            for j in range(mask.shape[1]):
        
                x1 = j*ratio
                x2 = x1+ratio
                
                tmp = lccls[y1:y2,x1:x2]
                                
                if mask[i,j] == 1:
                    if np.any(tmp<len(clsAttributes))==True:
                        negdx = np.where(tmp<len(clsAttributes))
                        tmp[negdx] = 0
                    uniqcnt = np.unique(tmp)
                    clscnt = np.bincount(tmp.ravel())
                    
                    if type(uniqcnt).__name__ == 'int':
                        Nveg = 1
                    else:
                        Nveg = len(uniqcnt)
                        
                    cells+=1
                    
                else:
                    Nveg=0
                    uniqcnt = 0
                    clscnt = 0
                
                f.write('{0} {1}\n'.format(cnt,Nveg))
                    
                if Nveg != 0:
                    for t in range(uniqcnt.size):
                        vegcls = int(uniqcnt[t])
                        Cv = np.float(clscnt[uniqcnt[t]])/np.float(clscnt.sum())
                        attributes = clsAttributes[vegcls]['properties']
                        rdepth1 = str(attributes['rootd1'])
                        rfrac1 = str(attributes['rootfr1'])
                        rdepth2 = str(attributes['rootd2'])
                        rfrac2 = str(attributes['rootfr2'])
                        rdepth3 = str(attributes['rootd3'])
                        rfrac3 = str(attributes['rootfr3'])
                    
                        f.write('\t{0} {1} {2} {3} {4} {5} {6} {7}\n'.format(vegcls,Cv,rdepth1,rfrac1,rdepth2,rfrac2,rdepth3,rfrac3))           
                
                cnt+=1
                
    return