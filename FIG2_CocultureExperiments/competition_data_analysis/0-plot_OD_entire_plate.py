import numpy as np
import matplotlib.pyplot as plt
import pickle
import pandas as pd
from scipy.ndimage import gaussian_filter1d

import matplotlib.cm as cm

plot_log_OD = False

conditions = 6
cmap = cm.get_cmap('viridis',conditions) #creates a cmap with as many colors as conditions
cmap = cm.get_cmap('tab10')

label_dict = {} #for each condition, we'll label it with the corresponding gluc. conc and initial OD
YFP_percentages = [100,90,60,40,10,0]
for condition in range(conditions):
    which_percentage = YFP_percentages[condition]

    label_dict[condition] = which_percentage

label_dict[-1] = 'Control'
    
def find_kneepoint(OD,time):
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

    return where_kneepoint+1

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

figures_dir = 'figures'
rows = 8
cols = 12

filename = '29.01.2025.yfp-k12.victor.xlsx'

data = pd.read_excel(filename)

time = list(data.iloc[3:86,0].to_numpy())
time_new = [float(t[:-1]) for t in time] #eliminate the 's' of secs and turn to float
time = np.round(np.array(time_new)/3600,2) #in hourse
t0 = 0 #index from which to start plotting.
ODs = data.iloc[3:86,2:]
ODs = ODs.to_numpy() #Rows are timepoints, cols are positions in the well, ordered by column first
### The controls are columns 1 and 12, minus the wells that have bubbles
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


OD_0 = 0.15
OD_0_theoretical = np.array([1,1/3,1/9,1/27])*OD_0

t_sat_vector = []
plotted_conditions = []
fig,ax = plt.subplots(8,12,figsize=(8,4))
for i in range(np.shape(ODs)[1]):
    plot_row = i//cols
    plot_col = i%cols


    condition = which_condition(i)
    alpha = 1-(plot_row//2)/5

    OD_rescaled_by_control = (ODs[t0:,i] #OD as is
                                -control_means_in_time[t0:]) #remove control

    OD_rescaled_by_control_and_normalized_to_Delta_OD = OD_rescaled_by_control
    OD_rescaled_by_control_and_normalized_to_Delta_OD = np.array(OD_rescaled_by_control_and_normalized_to_Delta_OD,
                                                                 dtype=float)
    if plot_log_OD:
        ax[plot_row,plot_col].plot(time[t0:],
                                   np.log(OD_rescaled_by_control_and_normalized_to_Delta_OD)
                                    ,
                                   linewidth=3,color=cmap(condition) if condition !=-1 else 'k',
                                    label = label_dict[condition] if condition not in plotted_conditions else None,
                                   alpha = alpha)
    else:
        ax[plot_row,plot_col].plot(time[t0:],
                                   OD_rescaled_by_control_and_normalized_to_Delta_OD
                                    ,
                                   linewidth=3,color=cmap(condition) if condition !=-1 else 'k',
                                    label = label_dict[condition] if condition not in plotted_conditions else None,
                                   alpha = alpha)
    if condition not in plotted_conditions: plotted_conditions.append(condition)
    ax[plot_row,plot_col].set_xticks([])
    ax[plot_row,plot_col].set_yticks([])
    if plot_log_OD==False:  
        ax[plot_row,plot_col].set_ylim(0,1)
    
    fig.legend(
               bbox_to_anchor=(0.98,1),frameon=True,
              title='%YFP',
               fontsize='small')

    ''' Find the kneepoint and scatter-plot it in red '''
    OD_filtered = gaussian_filter1d(np.array(OD_rescaled_by_control_and_normalized_to_Delta_OD,dtype=float),sigma=4)
    t_sat = find_kneepoint(OD_filtered,time)
    t_sat_vector.append(t_sat)
    if (plot_row,plot_col) == (2,2):
        '''
        This sample includes a large popping bubble and the algorithm mistakenly takes the popping time as
        as the saturation.
        Apply the algorithm from t>t_bubble_pop to ensure it gets the correct saturation point
        '''
##        print('Handling bubbly sample with popping')
        OD_filtered = gaussian_filter1d(np.array(OD_rescaled_by_control_and_normalized_to_Delta_OD,dtype=float),sigma=4)
        initial_cutoff = len(OD_filtered)//3
        OD_filtered = OD_filtered[initial_cutoff:] #chosen by inspection
        time_filtered = time[initial_cutoff:]
        t_sat = initial_cutoff+find_kneepoint(OD_filtered,time_filtered)

    if plot_col!=0 and plot_col!=11: #do not plot for empty controls
        if plot_log_OD:
            ax[plot_row,plot_col].scatter(time[t_sat],np.log(OD_rescaled_by_control_and_normalized_to_Delta_OD[t_sat]),
                                  c='red')
        else:
            ax[plot_row,plot_col].scatter(time[t_sat],OD_rescaled_by_control_and_normalized_to_Delta_OD[t_sat],
                                  c='red')
   
ax[0,0].set_yticks([0,1])
ax[0,0].set_ylabel(r'$ OD$')
ax[-1,0].set_xticks([time[0],time[-1]])
ax[-1,0].set_xlabel('Time (h)')
plt.suptitle(filename)
plt.savefig('OD.png')


read_fluorescence=True
if read_fluorescence==True:
    fluo = data.iloc[86:169,2:]
    
    fluo = fluo.to_numpy(dtype=float) #Rows are timepoints, cols are positions in the well

    controls = np.hstack((fluo[:,:8],fluo[:,-8:]))
    control_means_by_row = np.zeros((len(time),8))
    for i in range(8):
        control_means_by_row[:,i] = 0.5*(controls[:,i]+controls[:,i+8])
    control_means_in_time = np.mean(controls,axis=1) #avg of all controls at each point in time.
    control_means_total = np.mean(controls)

    
    plotted_conditions = []
    fig,ax = plt.subplots(8,12,figsize=(8,4))
    for i in range(np.shape(fluo)[1]):
        plot_row = (i)//cols
        plot_col = (i)%cols

        ### Remove previously calculated control average.
##        fluo_to_plot = np.array(fluo[:,i]-control_means_by_row[:,plot_row],dtype=float)
        ### Or calculate the average now.
        fluo_to_plot = (gaussian_filter1d(fluo[:,i],4)
                        -0.5*gaussian_filter1d(fluo[:,plot_row*cols],4)#subtract avg of left control
                        -0.5*gaussian_filter1d(fluo[:,(plot_row+1)*cols-1],4)) #and right control
        condition = which_condition(i)
        alpha = 1-(plot_row//2)/5
        ax[plot_row,plot_col].plot(time[t0:],
                                   fluo_to_plot[t0:],
                                   linewidth=3,
                                   color=cmap(condition) if condition !=-1 else 'k',
                                   label = label_dict[condition] if condition not in plotted_conditions else None,
                                   alpha=alpha)
        if condition not in plotted_conditions: plotted_conditions.append(condition)
        ax[plot_row,plot_col].set_xticks([])
        ax[plot_row,plot_col].set_yticks([])
        ax[plot_row,plot_col].set_ylim(-10,45)

        ''' Find the kneepoint and scatter-plot it in red '''
        ### Use the saturation point found in the OD curve.
##        if plot_col!=0 and plot_col!=11: #do not plot for empty controls
##            ax[plot_row,plot_col].scatter(time[t_sat_vector[i]],fluo_to_plot[t_sat_vector[i]],
##                                          c='red')
        
    ax[0,0].set_yticks([-10,45])
    ax[0,0].set_ylabel('YF')
    ax[-1,0].set_xticks([time[0],time[-1]])
    ax[-1,0].set_xlabel('Time (h)')
    fig.legend(
               bbox_to_anchor=(0.98,1),frameon=True,
              title='%YFP',
               fontsize='small')
    plt.suptitle(filename)
    
plt.savefig('YF.png')
plt.show()
