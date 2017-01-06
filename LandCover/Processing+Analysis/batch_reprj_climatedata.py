import arcpy
import glob

arcpy.CheckOutExtension("Spatial")
arcpy.env.overwriteOutput= True

template = r'D:\Kel\UAH\classes\ESS\thesis\gis\LandScan\LandScan2014_BB.tif'
tmpras = arcpy.Raster(template)
arcpy.env.snapRaster = tmpras
arcpy.env.extent = tmpras
arcpy.env.outputCoordinateSystem = tmpras
arcpy.env.cellSize = tmpras

outcoor = template
resamp = 'NEAREST'
csize = 1000

outprjfolder = r'D:\Kel\UAH\classes\ESS\thesis\gis\climate\NEX\prj\\'

inpath = r'D:\Kel\UAH\classes\ESS\thesis\gis\climate\NEX\org\*.tif'

files = glob.glob(inpath)

for f in files:
    
    prjout = outprjfolder+f.split('\\')[-1][:-4] + '_prj.tif'
    rasout = outprjfolder[:-5]+f.split('\\')[-1][:-4] + '_prj_BB.tif'
    
    arcpy.ProjectRaster_management(f,prjout,out_coor_system=outcoor,resampling_type=resamp,cell_size=csize)
    arcpy.ProjectRaster_management(f,prjout,out_coor_system=outcoor)
    
    outras = arcpy.sa.ExtractByMask(f,template)
    
    outras.save(rasout)
    
    del outras