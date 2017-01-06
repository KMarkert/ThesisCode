import netCDF4
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from matplotlib.patches import Rectangle

yrs = np.arange(2015,2055,5)

clmdata = np.zeros([3,138,98])

var = ['ppt','runoff','base']

for i in range(len(var)):

    ncfile = r'~\models\output\netCDFs\{0}_climatology_clim.nc'.format(var[i])
    
    nc = netCDF4.Dataset(ncfile)
    
    clmdata[i,:,:] = nc.variables['YrClimo'][:,:]
    
    lat = nc.variables['lat'][:]
    lon = nc.variables['lon'][:]
    
    nc.close()
    
clmdata = np.ma.masked_where(clmdata==-9999.,clmdata)
    
lons,lats = np.meshgrid(lon,lat)

elas = np.zeros([yrs.size,138,98])

prjdata = np.zeros([3,138,98,yrs.size])

for y in range(yrs.size):
    
    if (yrs[y])%4 == 0:
        days = 366
    else:
        days = 365
    
    yelas = np.zeros([days,138,98])
    
    if (yrs[y]-1)%4 == 0:
        idx = 366
    else:
        idx = 365
    
    for i in range(len(var)):
        ncfile = r'~\models\output\netCDFs\{0}_{1}_org.nc'.format(var[i],yrs[y])
        
        nc = netCDF4.Dataset(ncfile)
        
        data = nc.variables[var[i]][:,:]
    
        prjdata[i,:,:,y] = np.sum(data[idx:,:,:],axis=0) 
    
        nc.close()
    
    elas[y,:,:] = ((prjdata[1,:,:,y] - clmdata[1,:,:])/(prjdata[0,:,:,y] - clmdata[0,:,:])) * (clmdata[0,:,:]/clmdata[1,:,:])
    
    print elas[y,:,:].max(),elas[y,:,:].min()
    
    #dSpct[y,:,:] = np.ma.masked_where(dSclm<-100,dSpct[y,:,:])
    elas[y,:,:] = np.ma.masked_where(clmdata[0,:,:].mask==True,elas[y,:,:])
    

elasmed = np.zeros(yrs.size)

for i in range(yrs.size):
    tmp = elas[i,:,:]
    elasmed[i] = np.median(tmp[np.where(clmdata[0,:,:].mask==False)])
    
Z = [[0,0],[0,0]]
levels = np.arange(-5,8)
cb = plt.contourf(Z, levels, cmap='bwr_r')
plt.close()

fig, ax = plt.subplots()

cnt = 0

m = Basemap(projection='aea', llcrnrlat=lat.min()-0.5,urcrnrlat=lat.max()+0.5,
    llcrnrlon=lon.min()-0.5,urcrnrlon=lon.max()+0.5,lon_0=lon.mean(),lat_0=lat.mean(),
    resolution='i',ax=ax)

xx,yy = m(lons,lats)

xlabels = [1,1,1,1]

ylabels = [1,1,1,1]

m.drawparallels(np.arange(0,90,5),labels=ylabels,fontsize=8)
m.drawmeridians(np.arange(0,180,5),labels=xlabels,fontsize=8)

m.drawcountries(linewidth=0.5,zorder=11)
m.drawcoastlines(linewidth=0.5,zorder=12)
m.fillcontinents(color='silver')

emed = np.median(elas,axis=0)

emed = np.ma.masked_where(clmdata[0,:,:].mask==True,emed)

cb1 = m.pcolormesh(xx,yy,emed,cmap='bwr_r',vmin=-5,vmax=7,zorder=10)

#ax.set_title(yrs[cnt],fontsize=10)
ax.axis('off')

cnt+=1
        
plt.subplots_adjust(bottom=0.255)

ax1 = fig.add_axes([0.1, 0.04, 0.8, 0.1])
    
ax1.plot(yrs,elasmed,'b',linewidth=2)
#ax1.plot(yrs[np.where(dSmed>=0)],dSmed[np.where(dSmed>=0)],'b',linewidth=2)
#ax1.plot([2032.1,2035],[0,dSmed[4]],'b',linewidth=2)
ax1.plot([yrs[0],yrs[-1]],[1,1],'silver')
#ax1.set_ylim(-30,30)
ax1.set_ylabel('$\mathrm{\\varepsilon}_\mathit{P}$',fontsize=10)
ax1.set_yticks([1,2,3])
ax1.set_yticklabels([1,2,3],fontsize=10)
ax1.set_xticklabels(np.arange(2015,2055,5),fontsize=10)

ax3 = fig.add_axes([0.2, 0.2125, 0.6, 0.01])
cbar = fig.colorbar(cb,cax=ax3,orientation='horizontal',ticks=np.arange(-2,8))
cbar.set_label('$\mathrm{\\varepsilon}_\mathit{P}$',labelpad=-2)
cbar.ax.set_xlim([cbar.norm(-2), cbar.norm(7)])
cbar.outline.set_clip_box([0, 0, 1, 1])

plt.show()

plt.savefig(r'~\runoff_elasticity_changes_org.jpg',dpi=500)
