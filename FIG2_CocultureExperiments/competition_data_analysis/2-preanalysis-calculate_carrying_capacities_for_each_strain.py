import numpy as np
import matplotlib.pyplot as plt
import pickle
import pandas as pd
from scipy.ndimage import gaussian_filter1d
from scipy.optimize import curve_fit
import matplotlib.cm as cm

'''
Load Data
'''
rows = 8
cols = 12

filename = '29.01.2025.yfp-k12.victor.xlsx'

data = pd.read_excel(filename)

time = list(data.iloc[3:86,0].to_numpy())
time_new = [float(t[:-1]) for t in time] #eliminate the 's' of secs and turn to float
time = np.round(np.array(time_new)/3600,2) #in hourse
t0 = 0 #index from which to start plotting.
time = time[t0:]
ODs = data.iloc[3:86,2:]
ODs = ODs.to_numpy() #Rows are timepoints, cols are positions in the well
ODs = ODs[t0:,:] #from t0 onwards

which_controls_are_bad = [0,48,59,71,95]
which_controls_are_good = [11,12,
                           23,24,
                           35,36,
                           47,
                           60,
                           73,
                           83,84]


controls = np.vstack((ODs[:,which_controls_are_good]))
control_means_in_time = np.mean(controls,axis=1) #avg of both controls at each point in time.
fluo = data.iloc[86:169,2:]
fluo = fluo.to_numpy(dtype=float) #Rows are timepoints, cols are positions in the well
fluo = fluo[t0:,:]


'''
Proceed
'''

conditions = 6
cmap = cm.get_cmap('tab10')

label_dict = {} #for each condition, we'll label it with the corresponding gluc. conc and initial OD
YFP_percentages = [100,90,60,40,10,0]
for condition in range(conditions):
    which_percentage = YFP_percentages[condition]

    label_dict[condition] = which_percentage

label_dict[-1] = 'Control'

'''
Consider different dynamics, both exponential growth and with different
exponent values
'''

def exponential_growth(t,r,x0):
    out = x0*np.exp(r*t)
    return out

def alfa_growth(t,r,x0,alfa):
    num = x0
    den = ((1-alfa)*r*t*x0**(alfa-1)+1)**(1/(alfa-1))
    out = num/den
    return out
    


def which_condition(i):
    '''
    Excluding controls, returns which condition the plot belongs to.
    We will color according to the %YFP in the co-culture, from 100% to 0%.
    Lighter color will mean less initial OD (set through alpha when plotting).
    '''
    #do the row, modulo 2, times 5, then add the column -1 modulo 2
    which_row = i//cols
    which_col = i%cols
    if which_col==0 or which_col==11:
        return -1
    if which_col == 1: condition=0
    elif which_col == 10: condition=5
    else:
        condition = which_col//2
    return condition



'''
Use the kneepoint method to find K_yfp and K_12.
Do not count bubbly samples.
Or samples that do not saturate.
'''
def find_kneepoint(OD):
    '''
    Returns the index corresponding to the guessed point of saturation.
    Input: a smoothed out version of the OD curve
    Output: the index corresponding to the point of minimum 2nd derivative, which in principle corresponds to the saturation
    '''

    '''
    Step 2: check the second derivative of the OD and find the minimum of the second derivative,
    which (in principle) belongs to the kneepoint of the OD curve.
    '''

    grad = np.gradient(OD,time)
    gradgrad = np.gradient(grad,time) #second derivative

    where_kneepoint = np.where(gradgrad==np.min(gradgrad))[0][0]

    return where_kneepoint+1 #the +1 seems to be better at giving the saturation point.

K_yfp = [] #Should be in units of Delta_OD
K_k12 = []

indices_yfp = [#1, #Remove A2, bubbly
           #13, #Remove B2, bubbly
           #25, #Remove C2, bubbly
           37,
           49,
           61,
           73,
           85]
indices_k12 = [10,
               22,
               34,
               46,
               58,
               70,
               #82, #not saturated
               #94  #not saturated
               ]

from scipy.ndimage import gaussian_filter1d        
for i in indices_yfp: #however many controls we have for YFP. THe first 3 controls are not used (bubbly)
    index = i
    which_row = index//cols
    which_col = index%cols
    
    OD = ODs[:,index]-control_means_in_time[:] #noralize by controls
    OD = np.array(OD,dtype=float)
    OD = OD-OD[0]#consider Delta_OD
    OD_filtered = gaussian_filter1d(OD,sigma=4) #this seems to be enough to filter out sudden jumps/bumps and still get the correct saturation point


    where_kneepoint = find_kneepoint(OD_filtered)
    '''
    Find K_saturation from the kneepoint.
    '''

    saturation_time = time[where_kneepoint]
    saturation_OD = OD[where_kneepoint] #its already Delta_OD
    K_yfp.append(saturation_OD)
    
for i in indices_k12: #however many controls we have for YFP. THe first 3 controls are not used (bubbly)
    index = i
    which_row = index//cols
    which_col = index%cols
    
    OD = ODs[:,index]-control_means_in_time[:] #noralize by controls
    OD = np.array(OD,dtype=float)
    OD = OD-OD[0]#consider Delta_OD
    OD_filtered = gaussian_filter1d(OD,sigma=4) #this seems to be enough to filter out sudden jumps/bumps and still get the correct saturation point


    where_kneepoint = find_kneepoint(OD_filtered)
    '''
    Step 3: find K_saturation from the kneepoint.
    '''

    saturation_time = time[where_kneepoint]
    saturation_OD = OD[where_kneepoint]
    K_k12.append(saturation_OD)


'''
Plot average and std of K_yfp and K_12.
'''

print('Carrying capacity for YFP  \n',
      np.mean(K_yfp),'+-',
      np.std(K_yfp,ddof=1))
print('Carrying capacity for K12  \n',
      np.mean(K_k12),'+-',
      np.std(K_k12,ddof=1))

with open('monoculture_yields.dat','wb') as f: #Values of yields.
    pickle.dump((np.mean(K_k12),np.std(K_k12,ddof=1),
                 np.mean(K_yfp),np.std(K_yfp,ddof=1)),f)
    f.close()



fig,ax = plt.subplots(figsize=(4,3))
plotted_labels = []
for K in K_yfp:
    label='YFP'
    ax.scatter(0,K,
               c='tab:blue',edgecolor='k',
               s=50,
               label=label if label not in plotted_labels else None)
    if label not in plotted_labels:plotted_labels.append(label)
ax.errorbar(0.01,np.mean(K_yfp),
            yerr=np.std(K_yfp,ddof=1),
            elinewidth=2,
            color='tab:blue',
            capsize=5)
ax.scatter(0.01,np.mean(K_yfp),c='tab:blue',s=50)
for K in K_k12:
    label='K12'
    ax.scatter(0.1,K,
               c='tab:brown',edgecolor='k',
               s=50,
               label=label if label not in plotted_labels else None)
    if label not in plotted_labels:plotted_labels.append(label)

ax.errorbar(0.11,np.mean(K_k12),
            yerr=np.std(K_k12,ddof=1),
            elinewidth=2,
            color='tab:brown',
            capsize=5)
ax.scatter(0.11,np.mean(K_k12),c='tab:brown',s=50)

ax.set_xticks([])
ax.set_xlim(-0.05,0.2)
ax.set_ylabel(r'$\Delta OD$ upon saturation')
ax.legend(loc='best',frameon=False,
          title='Strain')
plt.close('all')

fig,ax = plt.subplots(figsize=(1.5,1.5))
for i in range(2):
    if i==0: data = K_yfp; facecolor = 'khaki';
    elif i==1: data = K_k12; facecolor='peru';
    ax.bar(i,
           np.mean(data),
           facecolor=facecolor,edgecolor='k',
           width=1,linewidth=1,
        )
    ax.errorbar(i,
           np.mean(data),
           yerr=np.std(data),
           ecolor='k',elinewidth=2,
        )

ax.set_xticks([0,1])
ax.set_xticklabels(['YFP','K12'])
ax.set_ylabel(r'Yield $(OD)$')
plt.tight_layout()
plt.savefig('barplot-yields.svg',format='svg',transparent=True)
plt.savefig('barplot-yields.png',format='png',dpi=300,transparent=True)
plt.show()
 


