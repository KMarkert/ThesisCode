import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from sklearn.metrics import confusion_matrix

from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})

infile = '~/gis/validation/MODIS_VIIRS_Comparison_true.csv'
data = pd.read_csv(infile)
siml = np.array(data.MCD12Q1_A2)
true = np.array(data.VIIRS_2013)

#shts = ['Sim_2011_05y', 'Sim_2015_05y', 'Sim_2015_10y']
#titles = ['2011 5yr Prediction','2015 5yr Prediction','2015 10yr Prediction']
#
#infile = r'D:\Kel\UAH\classes\ESS\thesis\gis\validation\LCSim_confusion_data.xlsx'
#
#xl = pd.ExcelFile(infile)

fig,ax = plt.subplots(ncols=3,nrows=1,sharey=True,sharex=True,figsize=(10, 3.8))

fig = plt.subplot()

k = []
o = []
p = []
u = []

i = 0

#for i in range(3):
#    
#    s1 = xl.parse(shts[i])
#    
#    siml = np.array(s1.SimVal)
#    true = np.array(s1.TrueVal)
    
## Compute confusion matrix
cm = confusion_matrix(true,siml)

# Normalize the confusion matrix by row (i.e by the number of samples
# in each class)
cm_normalized = cm.astype('float') / cm.sum(axis=1)[:, np.newaxis]
cm_normalized1 = cm.astype('float') / cm.sum(axis=0)[:, np.newaxis]
#print('Normalized confusion matrix')
#print(cm)

cls = []
cls1 = []
cls2 = []
pes = []

for j in range(6):
    cls.append(cm[j,j])
    cls1.append(cm_normalized[j,j])
    cls2.append(cm_normalized1[j,j])
    pei = (np.sum(cm[j,:])*np.sum(cm[:,j]))/float(np.sum(cm))
    pes.append(pei)
    
acc = (float(np.sum(cls))/float(np.sum(cm)))*100

pa = acc/100.
pe = np.sum(pes) / float(np.sum(cm))
kap = (pa-pe)/(1-pe)

o.append(acc)
k.append(kap)
p.append(np.mean(cls1)*100)
u.append(np.mean(cls2)*100)

    
tlabels = ['Forest','Grassland','Agriculture','Water', 'Urban','Other']
im = ax[i].imshow(cm_normalized, interpolation='nearest', cmap='Blues',aspect='auto')
im.set_clim(vmin=0,vmax=1)
#ax[i].set_title(titles[i],fontsize=11)

tick_marks = np.arange(6)
ax[i].set_xlim(-0.5,5.5)
ax[i].set_ylim(5.5,-0.5)
ax[i].set_xticks(tick_marks)
ax[i].set_yticks(tick_marks)
ax[i].set_xticklabels(tlabels,fontsize=10,rotation=45)
ax[i].set_yticklabels(tlabels,fontsize=10)
if i == 0:
    ax[i].set_ylabel('Predicted Class',fontsize=11)
if i == 1:
    ax[i].set_xlabel('Observed Class',fontsize=11)


print "YEAR\t\t\tOverall\t\tProducer's\tUser's\t\tKappa"
#for i in range(3):
#    print titles[i]+':\t',o[i],'\t',p[i],'\t',u[i],'\t',k[i]
print 'MODISvVIIRS:',o[i],'\t',p[i],'\t',u[i],'\t',k[i]

#plt.tight_layout()

#fig.subplots_adjust(right=0.8)
#fig.subplots_adjust(left=0.15)
cbar_ax = fig.add_axes([0.92, 0.225, 0.015, 0.7])
cb = fig.colorbar(im, cax=cbar_ax)
cb.set_label('Normalized by Observation')
cb.set_ticks([0.0,0.5,1])

fig.subplots_adjust(bottom=0.25)

#plt.tight_layout()
plt.show()
#plt.savefig(r'D:\Kel\UAH\classes\ESS\thesis\images\graphs\LCSim_confusion_matrix.jpg',dpi=500)
