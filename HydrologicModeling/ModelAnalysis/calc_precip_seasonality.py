import netCDF4
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
from osgeo import gdal
from osgeo.gdalnumeric import *  
from osgeo.gdalconst import *

def calc_seasonality(arr):
    angles = np.deg2rad(np.array([15.8,44.9,74.0,104.1,134.1,164.2,194.3,224.9,255.,285.,315.1,345.2]))
    
    S = np.sum(arr*np.sin(angles))
    C = np.sum(arr*np.cos(angles))
    
    intensity = np.sqrt(S**2 + C**2)
    direction = np.rad2deg(np.arctan(S/C))
    
    seasonality = intensity / np.sum(arr)
    
    #print S,C
    
    if C < 0:
        direction = direction + 180
    elif (S<0) and (C>0):
        direction = direction + 360
    else:
        pass
        
    return direction,seasonality
    
def plot_seasonality(direc,seasn):

    values = np.array([2010,2015,2020,2025,2030,2035,2040,2045,2050])    
    jet = plt.get_cmap('inferno')
    jet.set_under(color='black')
    cNorm  = colors.Normalize(vmin=min(values), vmax=max(values))
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)
    
    Z = [[0,0],[0,0]]
    levels = np.arange(2010,2060,5)
    cb = plt.contourf(Z, levels, cmap=jet)
    plt.close()
    
    snames = ['Vientiane','Mukdahan','Pakse','Stung Treng']
    
    fig,ax = plt.subplots(nrows=2,ncols=2,subplot_kw=dict(projection='polar'))#figsize=(6.5,6.5),
    
    cnt = 0
    
    width = np.pi/30
    
    for i in range(2):
        for j in range(2):
            #print 'yay'
            for r in range(len(direc[cnt])):
                
                colorVal = scalarMap.to_rgba(values[r])
            
                #ax[i,j].bar(np.deg2rad(direc[cnt][r]), seasn[cnt][r], width=width, bottom=0.0,color=colorVal)
                ax[i,j].arrow(0,0,np.deg2rad(direc[cnt][r]),seasn[cnt][r],head_width=0.01,lw=2.5,color=colorVal)
            
            ax[i,j].set_xticks(np.deg2rad([1.0,31.6,59.2,89.8,119.3,149.9,179.5,210.1,240.7,270.2,300.8,330.4]))
            ax[i,j].set_rgrids([0.5,1.0], angle=75,fontsize=9)
            ax[i,j].set_title((snames[cnt]+'                                             '),fontsize=12)
    
            ax[i,j].set_xticklabels(['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'],fontsize=10)
            
            cnt+=1
            
    fig.subplots_adjust(bottom=0.15)
    
    #ax3 = fig.add_axes([0.2, 0.2125, 0.6, 0.01])
    
    ax3 = fig.add_axes([0.1, 0.075, 0.8, 0.02])
    
    tcb = fig.colorbar(cb,cax=ax3,orientation='horizontal')
    tcb.set_ticks(values+2.5)
    tcb.ax.set_xticklabels(['Historic\n  Mean','2015','2020','2025','2030','2035','2040','2045','2050',''])
    #fig.tight_layout()
    
    fig.savefig(r'~\precip_seasonality_org.jpg',dpi=500)
    plt.show()
    
#arr = np.array([4.01,3.48,2.69,1.30,0.48,0.11,0.01,0.02,0.19,0.74,1.57,4.09])
#
#seasons = calc_seasonality(arr)

mons = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']

nidx = {'Jan':['01',0,31],'Feb':['02',31,59],'Mar':['03',59,90],'Apr':['04',90,120],
       'May':['05',120,151],'Jun':['06',151,181],'Jul':['07',181,212],'Aug':['08',212,243],
       'Sep':['09',243,273],'Oct':['10',273,304],'Nov':['11',304,334],'Dec':['12',334,364]}
       
lidx = {'Jan':['01',0,31],'Feb':['02',31,60],'Mar':['03',60,91],'Apr':['04',91,121],
       'May':['05',121,152],'Jun':['06',152,182],'Jul':['07',182,213],'Aug':['08',213,244],
       'Sep':['09',244,274],'Oct':['10',274,305],'Nov':['11',305,335],'Dec':['12',335,365]}

yrs = np.arange(2015,2055,5)

#stns = [70103,11201,11903,11901,13101,13402,13801,13901,14501,19802]
stns = [11901,13402,13901,14501]

areaprecip = np.zeros([len(yrs),len(stns),12])

clmprecip = np.zeros([len(stns),12])

cncfile = netCDF4.Dataset(r'~\models\output\netCDFs\ppt_climatology_clim.nc')

cmn = cncfile.variables['MonClimo']

catchfile = r'~\gis\GCSfiles\StnCatchment\XXXX_Catchment_GCS_01deg.tif'

direc = []
seasn = []

for i in range(len(stns)):
    
    hydros = []
    
    stnID = 'Stn{0}'.format(str(stns[i]))
    
    ds = gdal.Open(catchfile.replace('XXXX',stnID),GA_ReadOnly)
    b1 = ds.GetRasterBand(1)
    mask = BandReadAsArray(b1)
    
    idx = np.where(mask>0)
    
    for j in range(len(mons)):
        tmp = cmn[j,:,:]
        clmprecip[i,j] = np.mean(tmp[idx])
        
    hydros.append(clmprecip[i,:])
    
    for yr in range(len(yrs)):
        ncfile = r'~\models\output\netCDFs\ppt_{0}_org.nc'.format(str(yrs[yr]))
        
        nc = netCDF4.Dataset(ncfile)
        
        data = nc.variables['ppt'][:,:,:]
        
        if (yrs[yr]-1)%4 == 0:
            dt=366
        else:
            dt=365
            
        if yrs[yr]%4 == 0:
            yidx = lidx
        else:
            yidx = nidx
        
        for m in range(len(mons)):
            
            midx = yidx[mons[m]]
        
            t1 = midx[1] + dt
            t2 = midx[2] + dt
            
            areaprecip[yr,i,m] = np.mean(np.sum(data[t1:t2,:,:],axis=0))
            
        hydros.append(areaprecip[yr,i,:])
        
    stndirec = []
    stnseasn = []
    
    for hy in hydros:
        d,s = calc_seasonality(hy)
        stndirec.append(d)
        stnseasn.append(s)
    
    direc.append(stndirec)
    seasn.append(stnseasn)
    
plot_seasonality(direc,seasn)
