import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import itertools
import pickle

def exponential_growth(t,r,x0):
    out = x0*np.exp(r*t)
    return out


rows = ['A','B','C','D','E','F','G','H']
cols = np.array(np.arange(1,13),dtype=str)

tags= np.array(list(itertools.product(cols,rows)),dtype=str)
for i in range(len(tags)):
    tag = tags[i]
    tags[i] = tag[1]+tag[0]

tags = tags[:,0]

print(tags)


day = '04'
month = '07'
year = '25'
medium = 'glu0.2-r3-sp-vpy-kinetic24h' #it is the only recording of dynamics of F2 in glucose.
##medium = 'cit0.2-r3-sp-vpy-kinetik24h'

filename = 'R3/{}-{}-{}-{}.xlsx'.format(day,month,year,medium)

data = pd.read_excel(filename) #data organized in rows, then cols: A1,A2,...,A12,B1,...


''' Convert data to array and reshape with plate shape '''
t_vec = data.iloc[:87,0]
t_vect = [t[:-1] for t in t_vec]
t_vec = np.array(t_vect,dtype=float)
t_saturation = 8 #in hours, at which we will stop fitting exponential behavior
saturation_index = np.sum(t_vec/3600<t_saturation)
print('Saturation index',saturation_index)


data = data.iloc[:,2:]
data = np.array(data,dtype=float)
data_F2 = data[:87,61]


from scipy.optimize import curve_fit
popt,pcov,infodict,mesg,ier = curve_fit(exponential_growth,
                                        t_vec[:saturation_index]/3600,
                                        data_F2[:saturation_index],
                                        full_output=True)
r,x0 = popt
std = pcov[0,0]**0.5
print('r',r,'+-',std)

derivative = (np.roll(data_F2,-1)-data_F2)[:-1]

n0 = data_F2[0]
nt = n0*np.exp(r*(t_vec-t_vec[0])/3600)

fig,ax = plt.subplots(figsize=(4,3))
ax.plot(t_vec/3600,
        data_F2,
        c='tab:red',
        label='Acinetobacter guillouiae')
ax.plot(t_vec/3600,
        nt,'k--',
        label='Exponential Fit')
ax.text(7,0.12,
        'g = {:.3f} $\pm$ {:.3f} (1/h)'.format(r,std))
ax.set_yscale('log')
ax.set_xlabel('Time (h)')
ax.set_ylabel('OD')
ax.legend(loc='best',frameon=False)
plt.tight_layout()
plt.savefig('f2_growth_rate_in_glucose.png',dpi=300,
            format='png',transparent = True)
plt.savefig('f2_growth_rate_in_glucose.svg',
            format='svg',transparent = True)

with open('Acid_specialist_glucose_growth_rate.dat','wb') as f:
    pickle.dump((r,std),
                f)
    f.close()


plt.show()

