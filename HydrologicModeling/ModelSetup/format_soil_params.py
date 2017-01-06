import numpy as np
import pandas as pd
from osgeo import gdal
from osgeo.gdalnumeric import *  
from osgeo.gdalconst import *

def get_soil_info(scls, sdata,subsoil):
    sdict = {1:[1490.,9218.4,11.2,0.19,0.07,0.92,0.43],2:[1520.,2608.8,10.98,0.36,0.14,0.82],
             3:[1570.,1257.6,12.68,0.53,0.22,0.60],4:[1420.,950.40,10.58,0.70,0.26,0.25],
             5:[1280.,2061.6,9.1,0.54,0.15,0.1],6:[1490.,472.80,13.6,0.67,0.33,0.40],
             7:[1600.,576.00,20.32,0.69,0.44,0.60],8:[1380.,1096.8,17.96,0.75,0.44,0.10],
             9:[1430.,424.80,19.04,0.74,0.46,0.35],10:[1570.,285.6,29.0,0.76,0.56,0.52],
             11:[1350.,708.0,22.52,0.76,0.51,0.1],12:[1390.,763.2,27.56,0.77,0.57,0.25]}

    usdaclss = sdata[1][np.where(sdata[0]==scls)]
    susdaclss = subsoil[1][np.where(subsoil[0]==scls)]
    try: usdacls = np.bincount(usdaclss).argmax()
    except ValueError: usdacls = 9   
    
    try: susdacls = np.bincount(susdaclss).argmax()
    except ValueError: susdacls = 9
      
    schar = sdict[usdacls]
    subchar = sdict[susdacls]
    
    tbden = np.nanmean(sdata[2][np.where(sdata[0]==scls)]) * 1000.
    sbden = np.nanmean(subsoil[2][np.where(subsoil[0]==scls)]) * 1000.
    
    if np.isnan(tbden)==True:
        tbden = 1381.
    if np.isnan(sbden)==True:
        sbden = 1419
    
    drn = sdata[3][np.where(sdata[0]==scls)]
    
    try: drncls = np.bincount(drn).argmax()
    except ValueError: drncls = 4
    
    drndict = {1:[0.23,0.001,0.80],2:[0.22,0.00105,0.81],3:[0.21,0.00115,0.84],4:[0.205,0.00125,0.85],5:[0.20,0.0013,0.86],6:[0.19,0.0014,0.88],7:[0.18,0.0015,0.9]}
    
    return schar,subchar,tbden,sbden,drndict[drncls]

band = 1

xlsxfile = r'~\gis\soils\HWSD_DATA.xlsx'

indata = pd.ExcelFile(xlsxfile)

s1 = indata.parse('HWSD_DATA')

soildata = [np.array(s1.MU_GLOBAL),np.array(s1.T_USDA_TEX_CLASS,dtype=np.int32),np.array(s1.T_BULK_DENSITY,dtype=np.float),np.array(s1.DRAINAGE,dtype=np.int32)]
subsoil = [np.array(s1.MU_GLOBAL),np.array(s1.S_USDA_TEX_CLASS,dtype=np.int32),np.array(s1.S_BULK_DENSITY,dtype=np.float)]

infiles = [r'~\gis\GCSfiles\LMB_mask_GCS_01deg_full.tif',
           r'~\gis\GCSfiles\HWSD_BB_GCS_01deg.tif',
           r'~\gis\GCSfiles\BB_ELV_hires_FILL_GCS_01deg_bands.TIF',
           r'~\gis\GCSfiles\annual_pr_GCS_01deg.tif',
           r'~\gis\GCSfiles\LMB__BasinSlope_GCS_01deg.tif']

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
    #
    data[:,:,i] = var[:,:]
        
    ds = None
    b1 = None
    
lons = np.linspace(lon0,lon1,98)
lats = np.linspace(lat0,lat1,138)

xx,yy = np.meshgrid(lons,lats)

yy = np.flipud(yy)

soilfile = r'~\models\input\soils.params'

cells = 0

with open(soilfile, 'a') as f:
    
    cnt = 1

    for i in range(138):
        for j in range(98):
            
            run = int(data[i,j,0])
            if run <= 0:
                run=0       
            if data[i,j,2] < 0:
                run = 0
            
            if run > 0:
                scls = data[i,j,1]
                #print scls
                tschar,sschar,tbden,sbden,drn = get_soil_info(scls,soildata,subsoil)
                cells+=1
            else:
                tschar = [-999,-999,-999,-999,-999]
                sschar = [-999,-999,-999,-999,-999]
                tbden = 1331.
                sbden = 1395.
                drn = [0.2,0.001,0.9]
            
            
            grdc = cnt
            lat = yy[i,j]
            lon = xx[i,j]
            infilt = drn[0]
            Ds = drn[1]
            Dsmax = data[i,j,4] * sschar[1]
            Ws = drn[2]
            c = 2
            expt = tschar[2]
            expt1 = sschar[2]
            tksat = tschar[1]
            sksat = sschar[1]
            phis = -999
            elev = data[i,j,2]
            depth = 0.10
            depth1 =0.35
            depth2 =1.05
            avg_t = 27
            dp = 4
            bub = 0.32*expt+0.43
            quartz = tschar[3]
            quartz1 = tschar[3]-0.05
            bulk_den = tbden
            bulk_den1 =sbden
            soil_den = 2650.
            soil_den1 = 2685
            off_gmt = lon * 24 / 360.
            wrc_frac = tschar[3]
            wrc_frac1 = sschar[3]
            wpwp_frac = tschar[4]
            wpwp_frac1 = sschar[4]
            rough = 0.001
            srough = 0.0005
            annprecip = data[i,j,3]
            resid = 0.065
            fs_act = 1
            init_moist = depth*wrc_frac*1000
            initmoist2 = depth1*wrc_frac1*1000
            initmoist3 = depth2*wrc_frac1*1000
            
            f.write('{0} {1} {2:.4f} {3:.4f} {4} {5} {6} {7} {8} {9} {10} {11} {12} {13} {14} {15} {16} {17} {18} {19} {20} {21} {22} {23} {24} {25} {26} {27} {28} {29} {30} {31} {32} {33} {34} {35} {36} {37} {38} {39} {40} {41} {42} {43} {44} {45} {46} {47} {48} {49} {50} {51} {52}\n'.format(run,
                                grdc,lat,lon,infilt,Ds,Dsmax,Ws,c,expt,expt1,expt1,tksat,sksat,sksat,
                                phis,phis,phis,init_moist,initmoist2,initmoist3,elev,depth,depth1,depth2,avg_t,dp,bub,bub,bub,quartz,quartz1,quartz1,
                                bulk_den,bulk_den1,bulk_den1,soil_den,soil_den1,soil_den1,off_gmt,wrc_frac,wrc_frac1,wrc_frac1,wpwp_frac,wpwp_frac1,wpwp_frac1,
                                rough,srough,annprecip,resid,resid,resid,fs_act))
                                                                
            cnt+=1
            
print "COMPLETE"
