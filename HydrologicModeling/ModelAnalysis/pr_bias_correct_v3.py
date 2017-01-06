import netCDF4
import numpy as np
import matplotlib.pyplot as plt

days = 365
mons = np.array(['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'])
yrs = np.arange(1981,2011)

nidx = {'Jan':['01',0,31],'Feb':['02',31,59],'Mar':['03',59,90],'Apr':['04',90,120],
    'May':['05',120,151],'Jun':['06',151,181],'Jul':['07',181,212],'Aug':['08',212,243],
    'Sep':['09',243,273],'Oct':['10',273,304],'Nov':['11',304,334],'Dec':['12',334,364]}
    
lidx = {'Jan':['01',0,31],'Feb':['02',31,60],'Mar':['03',60,91],'Apr':['04',91,121],
    'May':['05',121,152],'Jun':['06',152,182],'Jul':['07',182,213],'Aug':['08',213,244],
    'Sep':['09',244,274],'Oct':['10',274,305],'Nov':['11',305,335],'Dec':['12',335,365]}
    
b = [0.437765858722,0.829201263724,1.08704406569,1.31989494531,1.69123948463,1.15165720864,1.65220033753,1.34956520929,1.58507304895,0.868883042746,0.760661902419,0.42592117102]

dscl = np.zeros([days,138,98])
mscl = np.zeros([mons.size,138,98])

oprfilet1 = r'~\models\output\netCDFs\ppt_1981.nc'
oprfilet2 = r'~\models\output\netCDFs\ppt_1991.nc'
oprfilet3 = r'~\models\output\netCDFs\ppt_2001.nc'

oprnct1 = netCDF4.Dataset(oprfilet1)
oprvart1 = oprnct1.variables['ppt']
oprnct2 = netCDF4.Dataset(oprfilet2)
oprvart2 = oprnct2.variables['ppt']
oprnct3 = netCDF4.Dataset(oprfilet3)
oprvart3 = oprnct3.variables['ppt']

cprfilet1 = r'~\models\output\netCDFs\ppt_1981_clim.nc'
cprfilet2 = r'~\models\output\netCDFs\ppt_1991_clim.nc'
cprfilet3 = r'~\models\output\netCDFs\ppt_2001_clim.nc'

cprnct1 = netCDF4.Dataset(cprfilet1)
cprvart1 = cprnct1.variables['ppt']
cprnct2 = netCDF4.Dataset(cprfilet2)
cprvart2 = cprnct2.variables['ppt']
cprnct3 = netCDF4.Dataset(cprfilet3)
cprvart3 = cprnct3.variables['ppt']
clon = cprnct3.variables['X']
clat = cprnct3.variables['Y']

s98 = np.zeros([mons.size,138,98])
s03 = np.zeros([mons.size,138,98])

pthresh = 1.0

for i in range(dscl.shape[1]):
    print i
    for j in range(dscl.shape[2]):
        
        try:
            if oprvart1[0,i,j].mask == True:
                dscl[:,i,j] = -9999.0
                mscl[:,i,j] = -9999.0
                continue
        except AttributeError:
        
            dtmp = np.zeros([yrs.size,days])
            mtmp = np.zeros([yrs.size,mons.size])
            
            oidx = np.concatenate([oprvart1[:,i,j],oprvart2[:,i,j],oprvart3[:,i,j]])
            cidx = np.concatenate([cprvart1[:,i,j],cprvart2[:,i,j],cprvart3[:,i,j]])
            
            #oidx[np.where(oidx<pthresh)] = 0.0
            #cidx[np.where(cidx<pthresh)] = 0.0
            
            t = -1
            dt = -365
            
            for y in range(yrs.size):
                if yrs[y]%4 == 0:
                    julian = 366
                    dt+=366
                    idx = lidx
                else:
                    julian = 365
                    dt+=365
                    idx = nidx
                    
                for d in range(julian):
                    t += 1
                    
                    if d == 365:
                        continue
                    else:
                        
                        t1 = t-15
                        t2 = t+15
                    
                        if t1 < 0:
                            pobs = oidx[0:t2]==0,oidx[0:t2]
                            pscn = cidx[0:t2]<pthresh,cidx[0:t2]
                            
                            #dtmp[y,d] = np.mean(oidx[0:t2])/(np.mean(cidx[0:t2]).astype(np.float))
                        elif t2 > oidx.size-1:
                            pobs = oidx[t1:oidx.size-1]==0,oidx[t1:oidx.size-1]
                            pscn = cidx[t1:oidx.size-1]<pthresh,cidx[t1:oidx.size-1]
                            
                            #dtmp[y,d] = np.mean(oidx[t1:oidx.size-1])/(np.mean(cidx[t1:oidx.size-1]).astype(np.float))
                        else:
                            pobs = oidx[t1:t2]==0,oidx[t1:t2]
                            pscn = cidx[t1:t2]<pthresh,cidx[t1:t2]

                        
                        if d < 31:
                            bv = b[0]
                        if (d>=31) & (d<59):
                            bv = b[1]
                        if (d>=59) & (d<90):
                            bv = b[2]
                        if (d>=90) & (d<120):
                            bv= b[3]
                        if (d>=120) & (d<151):
                            bv = b[4]
                        if (d>=151) & (d<181):
                            bv = b[5]
                        if (d>=181) & (d<212):
                            bv = b[6]
                        if (d>=212) & (d<243):
                            bv = b[7]
                        if (d>=243) & (d<273):
                            bv = b[8]
                        if (d>=273) & (d<304):
                            bv = b[9]
                        if (d>=304) & (d<334):
                            bv = b[10]
                        if (d>=334):
                            bv = b[11]
                            
                        pscninv = pscn[1] ** bv
                            
                        try:                            
                            dtmp[y,d] = np.mean(pobs)/np.float(np.mean(pscninv))
                        except AttributeError:
                            dtmp[y,d] = 0.0
                            
            
                for m in range(mons.size):
                    midx = idx[mons[m]]
                    
                    t1 = midx[1] + dt
                    t2 = midx[2] + dt
                    
                    #pobs = np.sum(np.ma.masked_where(oidx[t1:t2]<0,oidx[t1:t2]))
                    #pscn = np.sum(np.ma.masked_where(cidx[t1:t2]<pthresh,cidx[t1:t2]))
                    
                    pobs = oidx[t1:t2]==0,oidx[t1:t2]
                    pscn = cidx[t1:t2]<pthresh,cidx[t1:t2]
                    
                    #sigobs = np.std(pobs)
                    #sigscn = np.std(pscn)
                    #muobs = np.mean(pobs)
                    #muscn = np.mean(pscn)
                    #
                    #b = (sigobs / muobs) - (sigscn / muscn)
                    #
                    pscninv = pscn[1] ** b[m]
                    
                    
                    try:                            
                        mtmp[y,m] = np.mean(pobs)/(np.float(np.mean(pscninv)))
                    except AttributeError:
                        mtmp[y,m] = 0.0
                        
                    if yrs[y] == 1998:
                        s98[m,i,j] = np.mean(pobs)/(np.float(np.mean(pscninv)))
                        
                    if yrs[y] == 2003:
                        s03[m,i,j] = np.mean(pobs)/(np.float(np.mean(pscninv)))
                    
                    #print mons[m],np.mean(pobs), np.mean(pscninv) , b, mtmp[y,m]
                    
                    #mtmp[y,m] = np.mean(oidx[t1:t2])/(np.mean(cidx[t1:t2]).astype(np.float))
                    
                        
            dscl[:,i,j] = np.mean(dtmp,axis=0)
            mscl[:,i,j] = np.mean(mtmp,axis=0)
            
            
            
#dscl[np.where(dscl<0)] = -9999.0
#mscl[np.where(mscl<0)] = -9999.0

ncfile = netCDF4.Dataset(r'~\models\input\forcing\chirps_climo\chips_nex_scaling.nc', "w")

ncfile.Conventions = "COARDS"
ncfile.history = "Created by Kel Markert (UAH,ATS)."
ncfile.production = "NEX Precip Linear Bias Correction"

#create dimensions
ncfile.createDimension("Longitude", clon[:].size)
ncfile.createDimension("Latitude", clat[:].size)
ncfile.createDimension("Time", days)
ncfile.createDimension("Month",mons.size)


#create variables
latvar = ncfile.createVariable("Y", float, ("Latitude",))
latvar.long_name = "Latitude"
latvar.units = "degrees_north"
latvar[:] = clat[:]

lonvar = ncfile.createVariable("X", float, ("Longitude",))
lonvar.long_name = "Longitude"
lonvar.units = "degrees_east"
lonvar[:] = clon[:]

timevar = ncfile.createVariable("T", float, ("Time",))
timevar.long_name = "Time"
timevar.units = "Julian Date "
timevar[:] = range(0, days)

dbdt_var = ncfile.createVariable('Monthly B Var', float, ("Month",))
dbdt_var.long_name = 'Monthly exponential scaling factor for NEX Precip'
dbdt_var.missing_value = -9999.0
dbdt_var.units = "unitless"
dbdt_var[:] = b

data_var = ncfile.createVariable('Daily Scaling Factor', float, ("Time","Latitude","Longitude"))
data_var.long_name = 'Daily smoothed scaling factor for NEX Precip'
data_var.missing_value = -9999.0
data_var.units = "unitless"
data_var[:] = dscl[:,:,:]

mdata_var = ncfile.createVariable('Monthly Scaling Factor', float, ("Month","Latitude","Longitude"))
mdata_var.long_name = 'Monthly scaling factor for NEX Precip'
mdata_var.missing_value = -9999.0
mdata_var.units = "unitless"
mdata_var[:] = mscl[:,:,:]

ncfile.close()
