import netCDF4
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from scipy import stats
import matplotlib.colors as colors

sims = ['climate','base','For05','For10','Agr05','Agr10']

sim = sims[0]

yrs = np.arange(2015,2055,5)

clmdata = np.zeros([4,138,98])

var = ['ppt','runoff','base','evap']

for i in range(len(var)):
    if sim == 'climate':
        ncfile = r'~/models/output/netCDFs/{0}_climatology_clim.nc'.format(var[i])
    else:
        ncfile = r'~/models/output/netCDFs/{0}_climatology.nc'.format(var[i])
        
    nc = netCDF4.Dataset(ncfile)
    
    clmdata[i,:,:] = nc.variables['YrClimo'][:,:]
    
    lat = nc.variables['lat'][:]
    lon = nc.variables['lon'][:]
    
    nc.close()
    
lons,lats = np.meshgrid(lon,lat)

clmdata = np.ma.masked_where(clmdata==-9999.,clmdata)

hydroclm = np.zeros([3,138,98])

hydroclm[0,:,:] = (clmdata[0,:,:] - (clmdata[1,:,:] + clmdata[2,:,:] + clmdata[3,:,:]))

hydroclm[1,:,:] = clmdata[1,:,:] / clmdata[0,:,:]

hydroclm[2,:,:] = clmdata[3,:,:]

prjdata = np.zeros([4,138,98,yrs.size])

for y in range(yrs.size):
    if (yrs[y]-1)%4 == 0:
        idx = 366
    else:
        idx = 365
    
    for i in range(len(var)):
        if sim == 'climate':
            ncfile = r'~/models/output/netCDFs/{0}_{1}_org.nc'.format(var[i],yrs[y])
        else:
            ncfile = r'~/models/output/netCDFs/LCSim/{0}/{1}_{0}_{2}.nc'.format(sim,var[i],yrs[y])
        
        nc = netCDF4.Dataset(ncfile)
        
        data = nc.variables[var[i]][:,:]
    
        prjdata[i,:,:,y] = np.sum(data[idx:,:,:],axis=0) 
    
        nc.close()
        
hydrovars = np.zeros([3,138,98,yrs.size])

for i in range(hydrovars.shape[0]):
    for j in range(yrs.size):
        if i == 0:
            hydrovars[i,:,:,j] = (prjdata[0,:,:,j] - (prjdata[1,:,:,j] + prjdata[2,:,:,j] + prjdata[3,:,:,j]))
        elif i == 1:
            hydrovars[i,:,:,j] = prjdata[1,:,:,j] / prjdata[0,:,:,j]
        else:
            hydrovars[i,:,:,j] = prjdata[3,:,:,j]
            
trends = np.zeros([3,138,98])

for v in range(hydrovars.shape[0]):
    for i in range(hydrovars.shape[1]):
        for j in range(hydrovars.shape[2]):
            trends[v,i,j] = (stats.theilslopes(hydrovars[v,i,j,:])[0]*40)/hydroclm[v,i,j]
                
trends = np.ma.masked_where(clmdata[:-1,:,:].mask==True,trends)

labels = ['(a) $\mathrm{\Delta}$S', '(b) R/P Ratio', '(c) ET']

zdata = np.zeros(trends.shape)

clevels = [-100,-75,-50,-25,-10,-5,-4,-3,-2,-1,0,1,2,3,4,5,10,25,50,75,100]

fig,ax = plt.subplots(ncols=3,sharex=True,sharey=True)

for i in range(3):
    m = Basemap(projection='aea', llcrnrlat=lat.min()-0.5,urcrnrlat=lat.max()+0.5,
            llcrnrlon=lon.min()-0.5,urcrnrlon=lon.max()+0.5,lon_0=lon.mean(),lat_0=lat.mean(),
            resolution='i',ax=ax[i])
            
    xx,yy = m(lons,lats)
        
    xlabels = [0,0,0,1]        
 
    if i == 0:
        ylabels = [1,0,0,0]
    elif i == 2:
        ylabels = [0,1,0,0]
    else:
        ylabels = [0,0,0,0]
    
    m.drawparallels(np.arange(5,30,5),labels=ylabels,fontsize=8)
    m.drawmeridians(np.arange(90,115,5),labels=xlabels,fontsize=8)
    
    m.drawcountries(linewidth=0.5,zorder=11)
    m.drawcoastlines(linewidth=0.5,zorder=12)
    m.fillcontinents(color='silver')
    
    #if cnt == 0:
    #
    #    cb = m.pcolormesh(xx,yy,dSclm,cmap='seismic_r',vmin=-1000,vmax=1000)
    #    ax[i,j].set_title('Baseline',fontsize=10)
    #else:
    cb = m.pcolormesh(xx,yy,trends[i,:,:],cmap='bwr_r',norm=colors.BoundaryNorm(boundaries=clevels, ncolors=256),zorder=10)

    ax[i].text(0.5, 1700000,labels[i],fontsize=11)
    ax[i].axis('off')

ax3 = fig.add_axes([0.1, 0.2125, 0.8, 0.01])
cbar = fig.colorbar(cb,cax=ax3,orientation='horizontal',ticks=clevels)
cbar.ax.set_xticklabels(clevels,fontsize=11)
cbar.set_label('% $\mathrm{\Delta}$')
    
plt.show()
plt.savefig(r'~/{0}_hydrotrends_org.jpg'.format(sim),dpi=500)
