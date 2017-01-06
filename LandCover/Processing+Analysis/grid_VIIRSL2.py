import numpy as np
from scipy.interpolate import griddata
from scipy.spatial import ConvexHull
from pyhdf.SD import SD, SDC
from osgeo import gdal
from osgeo.gdalnumeric import *  
from osgeo.gdalconst import *
from pyproj import Proj, transform
from PIL import Image, ImageDraw

def update_vertices(xx,yy,hull):
    verts = hull.vertices
    
    vertsOut = []
    
    for i in range(verts.size):
        xval = hull.points[verts[i],0]
        yval = hull.points[verts[i],1]
        
        xidx = (np.abs(xx-xval)).argmin()
        yidx = (np.abs(yy-yval)).argmin()
        
        ridx = yidx / xx.shape[1]
        cidx = xidx % xx.shape[1]
        
        vertsOut.append((cidx,ridx))
        
    return vertsOut

fname = r'NPP_VSUT_L2.A2015241.0550.AGG_03000.2015242110711.hdf'
gname = r'NPP_VMAE_L1.A2015241.0550.AGG_03000.2015242104324.hdf'

HDFdata = SD(fname, SDC.READ)
HDFgeo = SD(gname, SDC.READ)

ST = HDFdata.select('SurfaceType')
Lat = HDFgeo.select('Latitude')
Lon = HDFgeo.select('Longitude')

st = ST[:,:]
lats = Lat[:,:]
lons = Lon[:,:]

HDFdata.end()
HDFgeo.end()

inProj = Proj(init='epsg:4326')
outProj = Proj(init='epsg:32648')

lonutm,latutm = transform(inProj,outProj,lons,lats)

latgrd = np.arange(latutm.min(),latutm.max(),1000)
longrd = np.arange(lonutm.min(),lonutm.max(),1000)

xx,yy = np.meshgrid(longrd,latgrd)

pts = np.zeros((lons.size,2))

refpts = np.vstack((xx.ravel(),yy.ravel())).T

pts[:,0] = lonutm.ravel()
pts[:,1] = latutm.ravel()

hull = ConvexHull(pts)

verts = update_vertices(xx,yy,hull)

coors = []

img = Image.new('L', (xx.shape[1], xx.shape[0]), 0)
ImageDraw.Draw(img).polygon(verts, fill=1)
mask = numpy.array(img)
mask = np.flipud(mask)

stgrd = griddata(pts, st.ravel(), (xx,yy), method='nearest')

stgrd = np.flipud(stgrd)

ndidx = np.where((stgrd>20) | (mask==0))
stgrd[ndidx] = 253

drv = gdal.GetDriverByName("GTiff")

outfile = 'gridded/NPP_VSUT_L2.A2015241.0550.AGG.TIF'

dsOut = drv.Create(outfile, xx.shape[1],  xx.shape[0], 1, gdal.GDT_Byte)

gt = [xx.min(), 1000, 0, yy.max(), 0, -1000]

dsOut.SetGeoTransform(gt)

#srs = osr.SpatialReference()
#srs.SetWellKnownGeogCS('WGS84')
#dest_wkt = srs.ExportToWkt()

dest_wkt = '''PROJCS["WGS 84 / UTM zone 48N",
                GEOGCS["WGS 84",
                DATUM["WGS_1984",
                    SPHEROID["WGS 84",6378137,298.257223563,
                        AUTHORITY["EPSG","7030"]],
                    AUTHORITY["EPSG","6326"]],
                PRIMEM["Greenwich",0,
                    AUTHORITY["EPSG","8901"]],
                UNIT["degree",0.01745329251994328,
                    AUTHORITY["EPSG","9122"]],
                AUTHORITY["EPSG","4326"]],
            UNIT["metre",1,
                AUTHORITY["EPSG","9001"]],
            PROJECTION["Transverse_Mercator"],
            PARAMETER["latitude_of_origin",0],
            PARAMETER["central_meridian",105],
            PARAMETER["scale_factor",0.9996],
            PARAMETER["false_easting",500000],
            PARAMETER["false_northing",0],
            AUTHORITY["EPSG","32648"],
            AXIS["Easting",EAST],
            AXIS["Northing",NORTH]]'''

dsOut.SetProjection(dest_wkt)

bandOut=dsOut.GetRasterBand(1)
bandOut.SetNoDataValue(253)
BandWriteArray(bandOut, stgrd)

dsOut = None
bandOut = None 
