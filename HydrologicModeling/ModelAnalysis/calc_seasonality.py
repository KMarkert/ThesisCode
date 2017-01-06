import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx

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
    
def plot_seasonality(direc,seasn,sim):

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
    
    fig.savefig(r'~\{0}_seasonality.jpg'.format(sim),dpi=500)
    plt.show()
    
#arr = np.array([4.01,3.48,2.69,1.30,0.48,0.11,0.01,0.02,0.19,0.74,1.57,4.09])
#
#seasons = calc_seasonality(arr)

yrs = [2015,2020,2025,2030,2035,2040,2045,2050]

#stns = [70103,11201,11903,11901,13101,13402,13801,13901,14501,19802]
stns = [11901,13402,13901,14501]

sims = ['climate','base','For05','For10','Agr05','Agr10']

sim = sims[5]

if sim == 'climate':
    Qsim = 'NEXsimscl'
else:
    Qsim = 'basesim'

xlsxfiles = [r'~\models\output\scenario_hydrographs\{0}_yr2015_hydrographs.xlsx'.format(sim),
             r'~\models\output\scenario_hydrographs\{0}_yr2020_hydrographs.xlsx'.format(sim),
             r'~\models\output\scenario_hydrographs\{0}_yr2025_hydrographs.xlsx'.format(sim),
             r'~\models\output\scenario_hydrographs\{0}_yr2030_hydrographs.xlsx'.format(sim),
             r'~\models\output\scenario_hydrographs\{0}_yr2035_hydrographs.xlsx'.format(sim),
             r'~\models\output\scenario_hydrographs\{0}_yr2040_hydrographs.xlsx'.format(sim),
             r'~\models\output\scenario_hydrographs\{0}_yr2045_hydrographs.xlsx'.format(sim),
             r'~\models\output\scenario_hydrographs\{0}_yr2050_hydrographs.xlsx'.format(sim)]

xls = []

for i in range(len(xlsxfiles)):
    xl = pd.ExcelFile(xlsxfiles[i])
    xls.append(xl)
    
direc = []
seasn = []

for stn in stns:
    
    hydros = []
        
    stnID = 'Stn{0}'.format(str(stn))
    
    Qfile = r'~\models\output\discharge_raw\{0}_{1}_Q.npy'.format(stnID,Qsim)
    Q = np.load(Qfile)
        
    monQs = np.array(np.split(Q,(Q.size/12.)))
        
    hydros.append(np.mean(monQs,axis=0))
    
    for yr in range(len(yrs)):
        sheet = xls[yr].parse(stnID)
        hydros.append(np.array(sheet.Simulated))
        
    stndirec = []
    stnseasn = []
    
    for hy in hydros:
        d,s = calc_seasonality(hy)
        stndirec.append(d)
        stnseasn.append(s)
    
    direc.append(stndirec)
    seasn.append(stnseasn)
    
plot_seasonality(direc,seasn,sim)
