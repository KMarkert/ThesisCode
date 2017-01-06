import os
import glob
import netCDF4
import numpy as np
from osgeo import gdal
from osgeo.gdalnumeric import *  
from osgeo.gdalconst import *
import scipy.stats as stats
import matplotlib.pyplot as plt

band =1

mons = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']

idx = {'Jan':['01',0,31],'Feb':['02',31,59],'Mar':['03',59,90],'Apr':['04',90,120],
       'May':['05',120,151],'Jun':['06',151,181],'Jul':['07',181,212],'Aug':['08',212,243],
       'Sep':['09',243,273],'Oct':['10',273,304],'Nov':['11',304,334],'Dec':['12',334,364]}
       
vicet = r'~\models\output\netCDFs\evap_2000.nc'

etnc = netCDF4.Dataset(vicet)
simet = etnc.variables['evap']

c = []
c1 = []

fobs = []
fsim = []
slps = []
intcs = []

dt = 731

for m in mons:
    
    midx = idx[m]
    
    t1 = midx[1] + dt
    t2 = midx[2] + dt

    jansim = np.sum(simet[t1:t2,:,:],axis=0)
    mask = jansim.mask
    
    modet = r'~\gis\MOD16\VICres\MOD16A2_ET_0.05deg_GEO_2002M{0}.tif_VICres.tif'.format(midx[0])
    
    ds = gdal.Open(modet,GA_ReadOnly)
    b1 = ds.GetRasterBand(band)
    janobs = BandReadAsArray(b1)
    ds = None
    b1 = None
    
    janobs = np.ma.masked_where(mask==True,janobs)
    janobs = np.ma.masked_where(jansim<0,janobs)
    jansim = np.ma.masked_where(jansim<0,jansim)
    janobs = np.ma.masked_where(janobs<0,janobs)
    jansim = np.ma.masked_where(janobs<0,jansim)
    
    janobs = janobs*0.1
    
    n = 0
    pr = 0
    
    c.append(jansim.mean())
    c1.append(janobs.mean())
    
    while (n < 300) and (pr < 0.2):
    
        ysamp = np.random.choice(jansim.shape[0],25,replace=False)
        xsamp = np.random.choice(jansim.shape[1],25,replace=False)
        
        samps = [[],[]]
    
        for i in range(25):
            for j in range(25):
                samps[0].append(ysamp[i])
                samps[1].append(xsamp[j])
        
        obvsamp = janobs[samps].ravel()
        simsamp = jansim[samps].ravel()
        
        obs = obvsamp[np.where(obvsamp.mask==False)]
        sim = simsamp[np.where(obvsamp.mask==False)]
        
        n = obs.size
    
    r = stats.pearsonr(obs,sim)
    bias = np.mean(obs-sim)
    pbias = np.sum(obs-sim)*100/np.sum(obs)
    rmse = np.mean(np.abs(obs-sim))
    mre = np.mean(np.abs(obs-sim)/obs)*100
    
    fobs.append(obs)
    fsim.append(sim)
    
    pr=r[0]
    
    slope, intercept, r_value, p_value, std_err = stats.linregress(obs,sim)
    
    slps.append(slope)
    intcs.append(intercept)
    
    diff = janobs - jansim
    
    wbias = diff.mean()
    wrmse = np.mean(np.abs(diff))
    
    print r[0], bias, pbias, rmse, mre
    #print wbias, wrmse

etnc.close()

fig,ax = plt.subplots(nrows=4,ncols=3,sharex=True,sharey=True,figsize=(6.5,8.5))

cnt = 0

for i in range(4):
    for j in range(3):
        
        x = np.arange(np.min(fobs[cnt]),np.max(fobs[cnt]+1))
        y = slps[cnt]*x+intcs[cnt]
        
        ax[i,j].plot(fobs[cnt],fsim[cnt],'ko',markersize=0.75)
        ax[i,j].plot([0,200],[0,200],'k--',label='1:1')
        ax[i,j].plot(x,y,'r',label='Best Fit')
        ax[i,j].set_xlim(0,200)
        ax[i,j].set_ylim(0,200)
        ax[i,j].set_title('{0} ET'.format(mons[cnt]),fontsize=10)
        
        cnt+=1
        
fig.text(0.5, 0.009, 'Observed ET [mm]', ha='center')
fig.text(0.009, 0.5, 'Simulated ET [mm]', va='center', rotation='vertical')

ax[0,-1].legend(loc='upper right',frameon=True,fontsize=8)

plt.xticks([0,50,100,150,200])
plt.yticks([0,50,100,150,200])

fig.tight_layout()
        
plt.show()

#fig.savefig(r'~\ET_validation.jpg',dpi=300)
