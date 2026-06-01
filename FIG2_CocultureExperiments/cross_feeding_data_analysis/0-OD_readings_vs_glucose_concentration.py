import pickle
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

cmap = plt.get_cmap('RdBu')
output_filename = 'OD_readings'

logarithmic = True
if logarithmic:
    output_filename = 'logarithmic-'+output_filename

fig,ax = plt.subplots(figsize=(5,3))
for ii,plate in enumerate(range(1,3)):
    filename = 'OD_readings_plate_{}.xlsx'.format(plate)
    plate = pd.read_excel(filename)
    ODs = plate.iloc[0,1:].to_numpy().reshape((12,8)).T
    well_names = plate.columns[1:].to_numpy().reshape((12,8)).T
##    print(well_names) #just for troubleshooting
    control = np.mean(ODs[6,:]) #controls in the current plate
    
    if ii==0:
        concentrations = 1/1.5**np.arange(12)
    elif ii==1:
        concentrations = 1/1.5**np.arange(12,24)

    for col in range(12):
        color = cmap(col/24) if ii==0 else cmap((col+12)/24) #red/blue for high/low concnetraionts
##        print(ODs[:6,col])
        ax.scatter(
            [concentrations[col]]*6,
            ODs[:6,col]-control,
            color=color,
            s=50,edgecolor='k',
            alpha=0.5)
        mean = np.mean(ODs[:6,col]-control)
        std = np.std(ODs[:6,col]-control)
        ax.errorbar(concentrations[col],
                    mean,yerr=std,
                    color=color)
        ax.scatter(concentrations[col],
                   mean,
                   s=50,marker='s',
                   edgecolor='k',color=color)
                    

ax.set_xlabel(r'$C/C_0$')
ax.set_ylabel('OD Final')
##ax.plot(np.linspace(0,1,100),
##        np.linspace(0,0.5,100),
##        'k--')
if logarithmic:
    ax.set_xscale('log')
    ax.set_yscale('log')
plt.tight_layout()
plt.savefig('{}.svg'.format(output_filename),format='svg',transparent=True)
plt.savefig('{}.png'.format(output_filename),format='png',transparent=True)
plt.show()
