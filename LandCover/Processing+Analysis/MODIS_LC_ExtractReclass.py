import arcpy

arcpy.CheckOutExtension("Spatial")
arcpy.env.overwriteOutput = True

sr = arcpy.SpatialReference()
sr.factoryCode = 32648
sr.create()
arcpy.env.outputCoodinateSystem = sr

inpath = r'~\gis\MODIS\tiles\org\\'
outpath = r'~\gis\MODIS\tiles\org\\'

#rvals = []
#
#for i in range(0,256):
#    tmp = '{0:08b}'.format(i)
#    if tmp[-2:] == '00':
#        rvals.append([i,1])
#    elif tmp[-2:] == '01':
#        if tmp[-4:-2] == '00':
#            rvals.append([i,1])

rvals = [[0,3],[1,0],[2,0],[3,0],[4,0],[5,0],[6,1],[7,1],[8,1],[9,1],[10,1],[11,3],[12,2],[13,4],[14,2],[15,5],[16,5]]

remap = arcpy.sa.RemapRange(rvals)

arcpy.env.workspace = inpath

flist = arcpy.ListFiles('*.hdf')

nfiles = len(flist)

cnt = 0

for i in flist:
    cnt+=1
    
    scene = i[:24]
    
    outfile = outpath+scene+'LCReclass.TIF'
    
    extract = arcpy.ExtractSubDataset_management(i,subdataset_index='0')
    
    reclass = arcpy.sa.Reclassify(extract,'VALUE',remap,"NODATA")
    reclass.save(outfile)
    
    del extract
    del reclass
    
    print "Done processing file {0} of {1}".format(cnt,nfiles)
