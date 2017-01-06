import datetime
import pandas as pd
import numpy as np
from osgeo import gdal
from osgeo.gdalnumeric import *
from osgeo.gdalconst import *
from sklearn.linear_model import LogisticRegression
from sklearn.cross_validation import cross_val_score
#from get_prd import *

def get_model():
    #begin fitting logistic regression model
    infile = r'~/validation/random_samples_training_stratified_extract.xlsx'
    
    indata = pd.ExcelFile(infile)
    
    s1 = indata.parse('random_samples_training_stratif')
    
    target = np.array(s1.LC_11)
    preLC2 = np.array(s1.LC_06)
    preLC1 = np.array(s1.LC_01)
    R_dist = np.array(s1.RoadDis)
    Slope = np.array(s1.Slope)
    U_dist = np.array(s1.UrbDis)
    Prt_dis = np.array(s1.PrtaDis)
    Max_temp = np.array(s1.MaxTemp)
    Min_temp = np.array(s1.MinTemp)
    Ann_precip = np.array(s1.AnnPreci)
    strdis = np.array(s1.StrDis)
    Hand = np.array(s1.HAND)
    PopDen = np.array(s1.PopDen)
    Soil = np.array(s1.SoilCls)
    Aspect = np.array(s1.Aspect)
    Curva = np.array(s1.Curva)
    Eleva = np.array(s1.ELV)
    
    clf = LogisticRegression()
    
    X = np.zeros((target.size,13))
    
    X[:,0] = PopDen[:]
    X[:,1] = preLC2[:]
    X[:,2] = Aspect[:]
    X[:,3] = Eleva[:]
    X[:,4] = Hand[:]
    X[:,5] = Curva[:]
    X[:,6] = U_dist[:]
    X[:,7] = R_dist[:]
    X[:,8] = Prt_dis[:]
    X[:,9] = Max_temp[:]
    X[:,10] = Min_temp[:]
    X[:,11] = Slope[:]
    X[:,12] = Ann_precip[:]
    #X[:,13] = strdis[:]
    #X[:,14] = Soil[:]
    #X[:,15] = preLC1[:]
    
    #Xtrain = X[rsel,:]
    #Ttrain = target[rsel]
    #
    #target[rsel] = -1
    #
    #tsel = np.where(target>=0)
    #
    #Xtest = X[tsel[0],:]
    #Ttest = target[tsel[0]]
    
    clf.fit(X,target)
        
    scores = cross_val_score(clf, X, target)
    print scores.mean()
    #print clf.coef_
    
    #pred = clf.predict(Xtest)
    #
    #preds = np.zeros((250,5))
    #trues = np.zeros((250,5))
    #
    #for j in range(5):
    #    for i in range(250):
    #        cls = j
    #        if pred[i] == cls:
    #            preds[i,j] = 1
    #            
    #        if Ttest[i] == cls:
    #            trues[i,j] = 1
    #
    #aucs = []
    #
    #for i in range(5):
    #    fpr,tpr, _ = metrics.roc_curve(trues[:,i],preds[:,i])
    #    auc = metrics.auc(fpr,tpr)
    #    aucs.append(auc)
    #    
    #print aucs
    #print np.mean(aucs)
    
    del X
    del target

    return clf

#******* begin main program ********
drv = gdal.GetDriverByName('GTiff')

files = [r'~\gis\MODIS\MCD12Q1.A2006001.MOSAIC.LCReclass.1KM_BB.TIF',
         r'~\gis\MODIS\MCD12Q1.A2011001.MOSAIC.LCReclass.1KM_BB.TIF',
         r'~\gis\ELV\BB_GMTED_MEANELV300.TIF',
         r'~\gis\ELV\BB_GMTED_SLOPE.TIF',
         r'~\gis\ELV\hydro\BB_GMTED_HAND.TIF',
         r'~\gis\ELV\BB_GMTED_CURVATURE.TIF',
         r'~\gis\distance\Urban_distance_BB.tif',
         r'~\gis\distance\Road_distance_BB.tif',
         r'~\gis\distance\stream_distance_BB.tif',
         r'~\gis\distance\PrtArea_distance_BB.tif',
         r'~\gis\climate\NEX\Nor\annual_tasmin_rcp45_NorESM1-M_2015_prj_BB.tif',
         r'~\gis\climate\NEX\Nor\annual_tasmax_rcp45_NorESM1-M_2015_prj_BB.tif',
         r'~\gis\climate\NEX\Nor\annual_pr_rcp45_NorESM1-M_2015_prj_BB.tif',
         r'~\gis\LandScan\LandScan2014_BB.tif',
         r'~\gis\soils\HWSD_BB.tif',]
         
years = [2015,2020,2025,2030,2035,2040,2045,2050]
step = [2015,2020,2025,2030,2035,2040,2045,2050,2055]

clf = get_model()


for i in range(len(years)):

    t1 = datetime.datetime.now()

    print "Modeling for {0}".format(years[i])

    lcOut = get_prd(clf,files)

    fWater = r'~\gis\MODIS\WaterMask\MOD44W_WATER_2000_BB.tif'
    ds = gdal.Open(fWater, GA_ReadOnly)

    outfile = r'~\gis\LCSim\LandCover_prediction_{0}_baseline.tif'.format(years[i])

    dsOut = drv.Create(outfile, ds.RasterXSize, ds.RasterYSize, 1, gdal.GDT_Byte)
    CopyDatasetInfo(ds,dsOut)
    bandOut = dsOut.GetRasterBand(1)
    BandWriteArray(bandOut, lcOut)

    bandOut=None
    dsOut = None
    ds = None

    files[0] = files[1]
    files[1] = outfile
    files[6] = r'~\gis\climate\NEX\Nor\annual_tasmin_rcp45_NorESM1-M_{0}_prj_BB.tif'.format(step[i])
    files[7] = r'~\gis\climate\NEX\Nor\annual_tasmax_rcp45_NorESM1-M_{0}_prj_BB.tif'.format(step[i])
    files[8] = r'~\gis\climate\NEX\Nor\annual_pr_rcp45_NorESM1-M_{0}_prj_BB.tif'.format(step[i])

    dt = datetime.datetime.now()-t1

    print "Run time for {0}: {1}".format(i,dt)

print 'Complete'
