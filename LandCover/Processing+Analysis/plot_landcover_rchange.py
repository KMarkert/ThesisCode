import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.stats as stats

infile = r'~\gis\validation\rate_of_change.xlsx'

indata = pd.ExcelFile(infile)

insht = indata.parse('rchange')

yrs = np.array(insht.Year)
frst = np.array(insht.Frst)*100
grss = np.array(insht.Grss)*100
agri = np.array(insht.Agri)*100
watr = np.array(insht.Watr)*100
urbn = np.array(insht.Urbn)*100
othr = np.array(insht.Othr)*100

fslope, fintercept, fr_value, fp_value, fstd_err = stats.linregress(yrs,frst)
gslope, gintercept, gr_value, gp_value, gstd_err = stats.linregress(yrs,grss)
aslope, aintercept, ar_value, ap_value, astd_err = stats.linregress(yrs,agri)
wslope, wintercept, wr_value, wp_value, wstd_err = stats.linregress(yrs,watr)
uslope, uintercept, ur_value, up_value, ustd_err = stats.linregress(yrs,urbn)
oslope, ointercept, or_value, op_value, ostd_err = stats.linregress(yrs,othr)

plt.plot(yrs,frst,'darkgreen',linewidth=2,label='Forest')
plt.plot(yrs,grss,'limegreen',linewidth=2,label='Grassland')
plt.plot(yrs,agri,'goldenrod',linewidth=2,label='Agriculture')
plt.plot(yrs,watr,'blue',linewidth=2,label='Water')
plt.plot(yrs,urbn,'red',linewidth=2,label='Urban')
plt.plot(yrs,othr,'gray',linewidth=2,label='Other')

plt.plot(yrs,(yrs*fslope+fintercept),'darkgreen',linestyle='--',linewidth=1)
plt.plot(yrs,(yrs*gslope+gintercept),'limegreen',linestyle='--',linewidth=1)
plt.plot(yrs,(yrs*aslope+aintercept),'goldenrod',linestyle='--',linewidth=1)
plt.plot(yrs,(yrs*wslope+wintercept),'blue',linestyle='--',linewidth=1)
plt.plot(yrs,(yrs*uslope+uintercept),'red',linestyle='--',linewidth=1)
plt.plot(yrs,(yrs*oslope+ointercept),'gray',linestyle='--',linewidth=1)

plt.ylim(0,60)
plt.xlim(2000,2016)

plt.ylabel('Percent Area [%]')

plt.legend(loc='upper center',fontsize=12,frameon=False,ncol=3)

plt.show()

plt.savefig(r'D:\Kel\UAH\classes\ESS\thesis\images\graphs\lc_rchange.jpg',dpi=500)

print fslope*5, gslope*5, aslope*5, wslope*5, uslope*5, oslope*5
