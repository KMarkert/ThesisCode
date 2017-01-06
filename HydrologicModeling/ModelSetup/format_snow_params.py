import numpy as np
from osgeo import gdal
from osgeo.gdalnumeric import *  
from osgeo.gdalconst import *

band = 1

infiles = [r'~\gis\GCSfiles\LMB_mask_GCS_01deg_full.tif',
           r'~\gis\GCSfiles\BB_ELV_hires_FILL_GCS_hires_class.tif',
           r'~\gis\GCSfiles\BB_ELV_hires_FILL_GCS_hires.tif',
           r'~\gis\GCSfiles\BB_ELV_hires_FILL_GCS_hires_class_precip.tif'
          ]
   
ds = gdal.Open(infiles[0],GA_ReadOnly)
b1 = ds.GetRasterBand(band)
mask = BandReadAsArray(b1)  
ds = None
b1 = None

ds = gdal.Open(infiles[1],GA_ReadOnly)
b1 = ds.GetRasterBand(band)
elvcls = BandReadAsArray(b1)  
ds = None
b1 = None

ds = gdal.Open(infiles[2],GA_ReadOnly)
b1 = ds.GetRasterBand(band)
elv = BandReadAsArray(b1)  
ds = None
b1 = None

ds = gdal.Open(infiles[3],GA_ReadOnly)
b1 = ds.GetRasterBand(band)
pr = BandReadAsArray(b1)  
ds = None
b1 = None

elvcls = np.ma.masked_where(elvcls>10,elvcls)
elv = np.ma.masked_where(elv<0,elv)

cells = 0

snowfile = r'~\models\input\snow.params'

with open(snowfile, 'a') as f:
    
    cnt = 1
    
    for i in range(mask.shape[0]):
        y1 = i*20
        y2 = y1+20
        for j in range(mask.shape[1]):
    
            x1 = j*20
            x2 = x1+20
            
            tmp = elvcls[y1:y2,x1:x2]
            prtmp = pr[y1:y2,x1:x2]
                            
            if mask[i,j] > 0:
                
                uniqcnt = np.unique(tmp)
                clscnt = np.bincount(tmp.ravel())                
                
                f.write('{0} '.format(cnt))
                
                bands = []
                for c in range(8):
                    try:
                        idx = np.where(tmp==uniqcnt[c])
                        frac = np.float(idx[0].size) / 400.
                    except IndexError:
                        frac = 0
                    bands.append(frac)
                    
                bands = np.around(bands,4)
                
                bands = bands.astype(np.float)/np.float(bands.sum())
                        
                f.write('{0} {1} {2} {3} {4} {5} {6} {7} '.format(bands[0],bands[1],bands[2],bands[3],bands[4],bands[5],bands[6],bands[7]))
                
                for c in range(8):
                    try:
                        idx = np.where(tmp==uniqcnt[c])
                        muelv = np.mean(elv[y1:y2,x1:x2][idx])
                    except IndexError:
                        muelv = 0
                    f.write('{0} '.format(muelv))
                
                pfracs = []
                
                for c in range(8):
                    try:
                        idx = np.where(tmp==uniqcnt[c])
                        pfac = np.float(prtmp[idx].sum())/np.float(prtmp.sum())
                    except IndexError:
                        pfac = 0
                    pfracs.append(pfac)
                    
                pfracs = np.around(pfracs,4)
                
                pfracs = pfracs.astype(np.float)/np.float(pfracs.sum())
                    
                f.write('{0} {1} {2} {3} {4} {5} {6} {7}\n'.format(pfracs[0],pfracs[1],pfracs[2],pfracs[3],pfracs[4],pfracs[5],pfracs[6],pfracs[7]))
                    
            else:
                muelv = elv.mean()
                frac = 1
                f.write('{0} {1} {2} {3}\n'.format(cnt,frac,muelv,frac))
            
            cnt += 1
                
                
            
                
