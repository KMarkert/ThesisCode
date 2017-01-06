import arcpy
import glob

arcpy.env.cellSize = 0.01
arcpy.env.overwriteOutput = True
arcpy.CheckOutExtension("Spatial")

inpath= '~/gis/LCSim/*.tif'

outfolder = '~gis/LCSim/GSCLCSim/'

template = r'~\gis\GCSfiles\MCD12Q1.A2011001.MOSAIC.LCReclass.1KM_BB_GCS_gridalign.tif'

arcpy.env.outputCoordinateSystem = template
arcpy.env.extent = template
arcpy.env.snapRaster = template

flist = glob.glob(inpath)

for f in flist:
    fname = f.split('\\')[-1]
    
    tmp = outfolder + 'tmp.tif'
    
    outpath = outfolder + fname[:-4] +'_GSC.tif'
    
    arcpy.Resample_management(f, tmp, "0.01", "NEAREST")
    
    outras = arcpy.sa.ExtractByMask(tmp, template)
    
    outras.save(outpath)
    
    del outras
