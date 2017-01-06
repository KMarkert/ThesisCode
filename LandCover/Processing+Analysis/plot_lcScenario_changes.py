import numpy as np
from osgeo import gdal
from osgeo.gdalnumeric import *
from osgeo.gdalconst import *
import matplotlib.pyplot as plt

drv = gdal.GetDriverByName('GTiff')

band = 1

inmask = '~/LMB_raster_GSC_mask_001.tif'

ds = gdal.Open(inmask,GA_ReadOnly)
b1 = ds.GetRasterBand(band)
mask = BandReadAsArray(b1)

ds = None
b1 = None

sims = sims = ['base','For05','For10','Agr05','Agr10']

yrs = np.arange(2015,2055,5)

inpath = r'~\gis\LCSim\GSCLCSim\\'

curves = np.zeros([len(sims),len(yrs),6])

for i in range(len(sims)):
    for j in range(len(yrs)):
        
        if sims[i] == 'base':
            filename = 'LandCover_prediction_{0}_GSC.tif'.format(yrs[j])
        else:
            filename = 'LandCover_prediction_{0}_{1}per_GSC.tif'.format(yrs[j],sims[i])
        
        infile = inpath + filename
        
        ds = gdal.Open(infile,GA_ReadOnly)
        b1 = ds.GetRasterBand(band)
        var = BandReadAsArray(b1)
        
        for k in range(6):
            curves[i,j,k] = (len(np.where(var[np.where(mask==0)]==k)[0]) / float(len(np.where(mask==0)[0])))*100
        
        ds = None
        b1 = None       
        
fig,ax = plt.subplots(ncols=3,nrows=2,sharex=True,sharey=True)
cnt = 0

labels = ['Forest','Grassland','Agriculture','Water','Urban','Other']
colors = ['darkgreen','limegreen','goldenrod','blue','red','gray']
snames = ['(a) Baseline LCC', '(b) Forest 5%', '(c) Forest 10%', '(d) Agr 5%', '(e) Agr 10%']

for i in range(2):
    for j in range(3):
        if cnt < 5:
            for k in range(6):
                colorVal = colors[k]
                ax[i,j].plot(yrs,curves[cnt,:,k],color=colorVal,linewidth=2,label=labels[k])
            
            ax[i,j].set_ylim(0,60)
            ax[i,j].set_xlim(np.min(yrs),np.max(yrs))
            ax[i,j].set_xticklabels(yrs,rotation=45,fontsize=10)
            ax[i,j].set_title((snames[cnt]),fontsize=11)
            
        else:
            for k in range(6):
                colorVal = colors[k]
                ax[i,j].plot([0,0],[0,0],color=colorVal,linewidth=2,label=labels[k])
         
        cnt+=1
        
ax[1,2].axis('off')
ax[1,2].legend(loc='center',fontsize=10)


fig.text(0.04, 0.5, 'Percent Area [%]', va='center', rotation='vertical')

plt.show()

plt.savefig(r'~\LCScenario_changes.jpg',dpi=500)
