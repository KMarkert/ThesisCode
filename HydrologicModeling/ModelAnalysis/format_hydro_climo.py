import netCDF4
import numpy as np
import datetime

mons = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']

nidx = {'Jan':['01',0,31],'Feb':['02',31,59],'Mar':['03',59,90],'Apr':['04',90,120],
       'May':['05',120,151],'Jun':['06',151,181],'Jul':['07',181,212],'Aug':['08',212,243],
       'Sep':['09',243,273],'Oct':['10',273,304],'Nov':['11',304,334],'Dec':['12',334,364]}
       
lidx = {'Jan':['01',0,31],'Feb':['02',31,60],'Mar':['03',60,91],'Apr':['04',91,121],
       'May':['05',121,152],'Jun':['06',152,182],'Jul':['07',182,213],'Aug':['08',213,244],
       'Sep':['09',244,274],'Oct':['10',274,305],'Nov':['11',305,335],'Dec':['12',335,365]}
       
#var = ['runoff','swe','ppt','evap']
var = ['ppt']

for i in range(len(var)):

    filet1 = r'~/output/netCDFs/{0}_1980.nc'.format(var[i])
    
    nc = netCDF4.Dataset(filet1)
    draw = nc.variables[var[i]]
    
    clon = nc.variables['longitude'][:]
    clat = nc.variables['latitude'][:]
        
    yrs = np.arange(1980,2011)
        
    dt = -365
    vi = 0
    
    outdata = np.zeros([yrs.size,len(mons),clat.size,clon.size])
    outmon = np.zeros([len(mons),clat.size,clon.size])
    
    yrvars = np.zeros([yrs.size,clat.size,clon.size])
        
    for y in range(yrs.size):
        if yrs[y]%4 == 0:
            idx = lidx
            dt+=366
        else:
            idx = nidx
            dt+=365
            
        if dt >= 3652:
            dt = 0
            vi += 1
        
        for m in range(len(mons)):
            
            midx = idx[mons[m]]
        
            t1 = midx[1] + dt
            t2 = midx[2] + dt
            
            outdata[y,m,:,:] = np.sum(draw[t1:t2,:,:],axis=0)
            
    outdata = np.ma.masked_where(outdata<=0,outdata)
            
    outyr = np.sum(outdata,axis=1)
    
    outmon = np.mean(outdata,axis=0)
    
    avgyr = np.mean(outyr,axis=0)
    
    ncfile = netCDF4.Dataset(r'~/output/netCDFs/{0}_climatology.nc'.format(var[i]), "w")
    
    ncfile.Conventions = "CF-1.6"
    ncfile.title = "VIC hydrologic flux climatology"
    ncfile.source = 'VIC hydrologic model 4.2.d'
    ncfile.history = "Created using the NASA SERVIR VICUtils package. " + datetime.date.today().isoformat()
    ncfile.date_created = str(datetime.datetime.now())
    ncfile.references = "N/A"
    ncfile.comment = "N/A"
    
    #create dimensions
    ncfile.createDimension("longitude", clon[:].size)
    ncfile.createDimension("latitude", clat[:].size)
    ncfile.createDimension("Year",yrs.size)
    ncfile.createDimension("Month",len(mons))
    
    #create variables
    latvar = ncfile.createVariable("latitude", float, ("latitude",))
    latvar.long_name = "Latitude"
    latvar.units = "degrees_north"
    latvar[:] = clat[:]
    
    lonvar = ncfile.createVariable("longitude", float, ("longitude",))
    lonvar.long_name = "Longitude"
    lonvar.units = "degrees_east"
    lonvar[:] = clon[:]
    
    timevar = ncfile.createVariable("mon", float, ("Month",))
    timevar.long_name = "Month"
    timevar.units = "Calendar Month"
    timevar[:] = range(len(mons))
    
    yrvar = ncfile.createVariable("yr", float, ("Year",))
    yrvar.long_name = "Year"
    yrvar.units = "Year since 0000"
    yrvar[:] = yrs[:]
    
    data_rvar = ncfile.createVariable('monraw', float, ("Year","Month","latitude","longitude"))
    data_rvar.long_name = 'Monthly accumulated {0} for each year'.format(var[i])
    data_rvar.missing_value = -9999.
    data_rvar.units = "mm/mon"
    data_rvar[:] = outdata[:,:,:,:]
    
    data_yvar = ncfile.createVariable('yrraw', float, ("Year","latitude","longitude"))
    data_yvar.long_name = 'Yearly accumulated {0} for each year'.format(var[i])
    data_yvar.missing_value = -9999.
    data_yvar.units = "mm/yr"
    data_yvar[:] = outyr[:,:,:]
    
    data_mvar = ncfile.createVariable('monclimo', float, ("Month","latitude","longitude"))
    data_mvar.long_name = 'Monthly accumulated {0} climatology'.format(var[i])
    data_mvar.missing_value = -9999.
    data_mvar.units = "mm/mon"
    data_mvar[:] = outmon[:,:,:]
    
    data_avar = ncfile.createVariable('yrclimo', float, ("latitude","longitude"))
    data_avar.long_name = 'Yearly accumulated {0} climatology'.format(var[i])
    data_avar.missing_value = -9999.
    data_avar.units = "mm/yr"
    data_avar[:] = avgyr[:,:]
    
    ncfile.close()
