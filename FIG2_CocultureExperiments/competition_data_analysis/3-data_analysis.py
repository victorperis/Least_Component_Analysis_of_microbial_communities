import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.ndimage import gaussian_filter1d
from scipy.optimize import curve_fit
from scipy.integrate import cumulative_trapezoid
import pickle

'''
Function Used to lighten a given color in plt.
Taken from StackOverflow:
https://stackoverflow.com/questions/37765197/darken-or-lighten-a-color-in-matplotlib
'''
def adjust_lightness(color, amount=0.5):
    import matplotlib.colors as mc
    import colorsys
    try:
        c = mc.cnames[color]
    except:
        c = color
    c = colorsys.rgb_to_hls(*mc.to_rgb(c))
    return colorsys.hls_to_rgb(c[0], max(0, min(1, amount * c[1])), c[2])

'''
Do stuff needed to plot each experimental condition in one color
'''
conditions=6
label_dict = {} #for each condition, we'll label it with the corresponding gluc. conc and initial OD
YFP_percentages = [100,90,60,40,10,0]
for condition in range(conditions):
    which_percentage = YFP_percentages[condition]

    label_dict[condition] = which_percentage

label_dict[-1] = 'Control'

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
Theoretical Values of the OD at t=0
'''
OD_0 = 0.15
OD_0_theoretical = np.round(np.array([1,1/3,1/9,1/27])*OD_0,3)


'''

'''
bubbly_samples = [(1,2), #samples that showed signs of bubbling in the plate, disrupting OD measurements. These wre not considered for the analysis.
                  (1,3),
                  (2,2),
                  (2,3),
                  (3,9),
                  (5,8),
                  (6,7)]
other_good_replicates = [((0,2),(0,3)), #FOR EACH bubbly datapoint, look at the replicates that are not bubbly.
                         ((0,2),(0,3)),
                         ((3,2),(3,3)),
                         ((3,2),(3,3)),
                         ((2,8),(2,9),(3,8)),
                         ((4,8),(4,9),(5,9)),
                         ((6,6),(7,6),(7,7))]



figures_dir = 'figures'
### Geometry of the 96-well plate.
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
ODs = ODs[t0:,:] #from t0 onwards #is shape (timepoints,rows x columns)

which_controls_are_bad = [0,48,59,71,95] #due to bubbling, so OD measurements are not reliable. Overall index position inthe plate.
which_controls_are_good = [11,12,
                           23,24,
                           35,36,
                           47,
                           60,
                           72,
                           83,84]



controls = np.vstack((ODs[:,which_controls_are_good]))
control_means_in_time = np.mean(controls,axis=1) #avg of all controls at each point in time.
fluo = data.iloc[86:169,2:] #fluorescence measurements
fluo = fluo.to_numpy(dtype=float) #Rows are timepoints, cols are positions in the well
fluo = fluo[t0:,:]

fluo_controls = fluo[t0:,which_controls_are_good]
filtered_fluo_controls = np.zeros(np.shape(fluo_controls))
for c in range(len(which_controls_are_good)):
    fluo_control = fluo_controls[t0:,c]
    filtered_fluo_controls[:,c] = fluo_control
fluo_control_mean = np.mean(filtered_fluo_controls,axis=1) #mean fluorescence in time, of the controls


'''
Choose one color for each condition to check
reproducibility for the same condition.
'''
import matplotlib.cm as cm

conditions = 20
cmap = cm.get_cmap('tab10')

def which_condition(i):
    '''
    Excluding controls, returns which condition the plot belongs to.
    We will color according to the %YFP in the co-culture, from 100% to 0%.
    Lighter color will mean less initial OD (set through alpha when plotting).
    Condition:
    0 - YFP 100%
    1 - YFP 90%
    2 - YFP 60&
    3 - YFP 40%
    4 - YFP 10%
    5 - YFP 0%
    '''
    #do the row, modulo 2, times 5, then add the column -1 modulo 2
    which_row = i//cols
    which_col = i%cols
    if which_col==0 or which_col==11:
        return -1 #these are controls
    if which_col == 1: condition=0 #monoculture
    elif which_col == 10: condition=5 #monoculture
    else:
        condition = which_col//2
    return condition


################# GET YFP controls from which to map fluorescence to OD #################
indices = [#1,  #Remove A2, bubbly
           #13, #Remove B2, bubbly
           #25, #Remove C2, bubbly
           37,
           49,
           61,
           73,
           85]
t_threshold = 0

k_values = np.zeros((len(time)-1-t_threshold,4))
k_fits = []
intercept_fits = []
k1_fits = []
k2_fits = []
counter=0


'''
Find K_yfp as the kneepoint of the OD curves. The OD curves are rescaled by the controls
'''
def find_kneepoint(OD,time):
    '''
    Returns the index corresponding to the guessed point of saturation.
    Input: a smoothed out version of the OD curve
    Output: the index corresponding to the point of minimum 2nd derivative, which in principle corresponds to the saturation
    '''

    grad = np.gradient(OD,time)
    gradgrad = np.gradient(grad,time) #second derivative

    where_kneepoint = np.where(gradgrad==np.min(gradgrad))[0][0]

    return where_kneepoint+1 #The +1 correction gets a better estimate of the saturation point. (see by plotting Figure A below).

K_yfp = [] #Should be in units of Delta_OD
for i in indices: #however many controls we have for YFP. The first 3 controls are not used (bubbly)
    index = i
    which_row = index//cols
    which_col = index%cols
    
    OD = ODs[:,index]-control_means_in_time[:] #noralize by controls
    OD = np.array(OD,dtype=float)
    OD = OD-OD[0]#consider Delta_OD
    OD_filtered = gaussian_filter1d(OD,sigma=4) #this seems to be enough to filter out sudden jumps/bumps and still get the correct saturation point
    where_kneepoint = find_kneepoint(OD_filtered,time)
    
    '''
    Step 3: find K_saturation from the kneepoint.
    '''
    saturation_time = time[where_kneepoint]
    saturation_OD = OD[where_kneepoint]
    K_yfp.append(saturation_OD)
##    print(saturation_time,saturation_OD)

    ### Figure A
##    fig,ax = plt.subplots()
##    ax.plot(time,OD)
##    ax.scatter(saturation_time,saturation_OD,color='red')
##    plt.show()

'''
For the YFP moncultures,
Create a list of YF-vs-DeltaOD curves so that later we can find the mean curve and interpolate.
'''
curves = []
fluo_mean_error = []
for index in indices:
    which_row = index//cols
    which_col = index%cols


    OD = ODs[t_threshold:,index]-control_means_in_time[t_threshold:]
    '''
    Remove initial value of the dataset
    '''
    OD = OD-OD[0] #remove initial condition to have Delta_OD

    ### Rescale fluorescence by controls of the same row (which has the same conc. of glucose)
    YF = np.array(gaussian_filter1d(fluo[t_threshold:,index],4)
                  -gaussian_filter1d(fluo_control_mean[t_threshold:],4))

    fluo_error = YF-(fluo[t_threshold:,index]-fluo_control_mean[t_threshold:]) #difference between true measurement and low-pass-filtered measurement

    fluo_mean_error.append(fluo_error) #append all discrepancies

    curves.append((YF,OD))

fluo_mean_error = np.std(fluo_mean_error)

### Option 2. Interpolation.
'''
Try Interpolation in 1D:

For every YF-vs-OD curve, use interpolation in 1D to get DeltaOD(YF).
Then create a curve that is the average of all those curves, at the chosen OD intervals,
and create a new function that is the interpolation of that mean curve.

We can then use this function to get DeltaOD(YF) for any measured YF.
'''
from scipy.interpolate import interp1d

#Create a common x-vector
x_min = min(min(x) for x, y in curves)
x_max = max(max(x) for x, y in curves)
x_common = np.linspace(x_min, x_max, 100)  # Fine grid
#Interpolate each curve on the common x-vector
y_interp = [] #interpolated curve for every real curve
for x, y in curves:
    f = interp1d(x, y, kind='linear',
                 fill_value=(y[0],y[-1]), #extrapolation is done by taking known extreme points.
                 bounds_error=False)
    y_interp.append(f(x_common)) #get interpolated curve at common x-vector.
    

# Step 3: Compute the average y-values
y_mean = np.mean(y_interp, axis=0) #average OD(YF) curve.
y_std = np.std(y_interp,axis=0,ddof=1) #so we have an idea of the error in the estimation of OD from YF.

# Plot all curves and the average curve
plt.figure(figsize=(4,3))
ii = 0
for x, y in curves:
    plt.plot(x, y,
             '-', alpha=0.5,
             label='Measured Curve' if ii==0 else None)
    ii +=1

plt.plot(x_common, y_mean, 'k-', linewidth=2,
         label='Average Curve')
plt.fill_between(x_common,y_mean-y_std,y_mean+y_std,
                 facecolor='lightgray',zorder=0)
plt.legend(loc='best',frameon=False)
plt.ylabel(r'$\Delta OD$')
plt.xlabel('YF (a.u.)')
plt.title(r'Average $\Delta$OD(YF) Curve')
plt.tight_layout()
plt.savefig('OD_from_YF_curve.png',format='png',transparent=True)
plt.savefig('OD_from_YF_curve.svg',format='svg',transparent=True)
plt.tight_layout()

### Get function so we can interpolate OD(YF) at any value of YF.
average_curve_function = interp1d(x_common,y_mean,kind='linear',
                                  fill_value='extrapolate')
### Error in OD comming from the mapping from YF to OD.
OD_YF_mapping_error = interp1d(x_common,y_std,kind='linear', 
                                fill_value = 'extrapolate')
### Error in OD coming from the uncertainty in the original measurement of YF, in the low-pass filtering
OD_error_from_YF_error = interp1d(x_common,np.gradient(y_mean,x_common)*fluo_mean_error,kind='linear',
                                fill_value = 'extrapolate')
                                  

'''
Check predictions vs real OD for the in-sample indices that we used to make the inference,
just to check results.
'''
fig,ax = plt.subplots()
ii=0
for index in indices:
    which_row = index//cols
    which_col = index%cols
    
    OD_true = ODs[t_threshold:,index]-control_means_in_time[t_threshold:]
    '''
    Remove OD[0] to get delta_OD from dataset
    '''
    OD_true = OD_true-OD_true[0]

    ### Rescale fluorescence by control of the same row (which has the same conc. of glucose)
    YF = np.array(gaussian_filter1d(fluo[t_threshold:,index],4)
                  -gaussian_filter1d(fluo_control_mean[t_threshold:],4))
##                            -0.5*gaussian_filter1d(fluo[t_threshold:,which_row*cols],4)#subtract avg of left control
##                            -0.5*gaussian_filter1d(fluo[t_threshold:,(which_row+1)*cols-1],4)) #and right control

    
    t_vect = time[t_threshold:]
    t_vect = (t_vect-t_vect[0])


    condition = which_condition(index)
    alpha = 1-(which_row//2)/5
    
    color = cmap(condition)


    '''
    Option 1: Use prediction through interpolation
    '''
    OD_interpolation = average_curve_function(YF)
    ax.plot(OD_true,OD_interpolation,
               color=color,
            alpha=alpha,
            linestyle='--',
            label = 'Interpolated' if ii==0 else None,
            )
    ii += 1

x_vect  = np.linspace(0,0.44,1000)
ax.plot(x_vect,x_vect,'k--')
ax.legend(loc='best',frameon=False)
ax.set_xlabel(r'$\Delta OD$ True')
ax.set_ylabel(r'$\Delta OD$ Predicted')

    



'''
Do True OUT-OF-SAMPLE predictions on DeltaOD(YF) in the co-cultures!
'''
'''
Out of sample indices INCLUDE bubbly samples, where the bubble popped during the 24h of growth.
They are not considered for the analysis of the inferred low variance mode, and are not included in the
final figure.
But the code here still considers them and tries to correct for the initial OD (not measurable because of the bubble)
by averaging over OD_0 of other replicates in the same conditions.
'''

out_of_sample_indices = [2,3,4,5,6,7,8,9,
                         14,15,16,17,18,19,20,21,
                         26,27,28,29,30,31,32,33,
                         38,39,40,41,42,43,44,45,
                         50,51,52,53,54,55,56,57,
                         62,63,64,65,66,67,68,69,
                         74,75,76,77,78,79,#80,81, #if not saturated, do not consider.
                         86,87,88,89,90,91,#92,93
    ]


saturation_times = []
plotted_conditions = []
plotted_dummies = []
fig,ax = plt.subplots()
YFP_data = []
K12_data = []
K12_error = []
data_conditions = [] #we will output also the condition of the data (ratio,OD_0) by index
''' Create a fake axis to plot a single color in various luminosities to plot initial OD '''
ax2 = ax.twinx()
color_list = [] #dummy list to be printed
label_value_list = [] #dummy list to be printed
alpha_list = [] #dummy list to be printed
for i,index in enumerate(out_of_sample_indices):
    condition = which_condition(i) #Percentage of YFP
    which_row = index//cols
    which_col = index%cols

    OD_total = ODs[t_threshold:,index]-control_means_in_time[t_threshold:]
    '''
    Consider DeltaOD.
    '''
    if (which_row,which_col) not in bubbly_samples:
        delta_OD_total = OD_total - OD_total[0] #OD curves are NOT low-pass filtered. The low-pass filter is only done to find the kneepoint time.
    else:
        '''
        Take the mean OD_0 of the replicates that are not bubbly and use that.
        '''
        where = [i  for i in range(len(bubbly_samples)) if bubbly_samples[i]==(which_row,which_col)][0]
        
        other_samples = other_good_replicates[where]
        OD_0 = 0
        for sample in other_samples:
            new_index = sample[0]*12+sample[1]
            OD_other = ODs[t_threshold:,new_index]-control_means_in_time[t_threshold:]
            OD_0 += OD_other[0]
        OD_0 /= len(other_samples)
        delta_OD_total = OD_total - OD_0
    
    ### Rescale fluorescence by control of the same row (which has the same conc. of glucose)
    YF = np.array(gaussian_filter1d(fluo[t_threshold:,index],4)
                            -0.5*gaussian_filter1d(fluo[t_threshold:,which_row*cols],4)#subtract avg of left control
                            -0.5*gaussian_filter1d(fluo[t_threshold:,(which_row)*cols-1],4)) #and right control
    
    t_vect = time[t_threshold:]
    t_vect = (t_vect-t_vect[0])

    #### Find saturation time and predict only until then
    OD_filtered = gaussian_filter1d(np.array(delta_OD_total,dtype=float),sigma=4)
    t_sat = find_kneepoint(OD_filtered,time)
##    print('tf',time[t_threshold+t_sat])
    if (which_row,which_col) == (2,2):
        '''
        This sample includes a large popping bubble and the algorithm mistakenly takes the popping time as
        as the saturation.
        Apply the algorithm from t>t_bubble_pop to ensure it gets the correct saturation point.
        This has no effect since the datapoint is not considered for analysis, but is necessary if plotting all datapoints including bubbly ones.
        '''
        OD_filtered = gaussian_filter1d(np.array(delta_OD_total,dtype=float),sigma=4)
        initial_cutoff = len(OD_filtered)//3
        OD_filtered = OD_filtered[initial_cutoff:] #chosen by inspection
        time_filtered = t_vect[initial_cutoff:]
        t_sat = initial_cutoff+find_kneepoint(OD_filtered,time_filtered)

    
    YF = YF[:t_sat]
    t_vect = t_vect[:t_sat]
    delta_OD_total = delta_OD_total[:t_sat]
    

    condition = which_condition(index)
    alpha = 1+(which_row//2)/5 #used to lighten the color later.

    '''
    Option 1: Use prediction of DeltaOD(YF), for the YFP strain,  through interpolation
    '''
    Delta_OD_interpolation = (average_curve_function(YF)) #OD(YFP) from interpolation
    error = np.sqrt(OD_error_from_YF_error(YF[-1])**2+OD_YF_mapping_error(YF[-1])**2) #err^2 = err1^2+err2^2
    Delta_OD_interpolation = Delta_OD_interpolation 

    Delta_OD_YFP = Delta_OD_interpolation[-1] #total change is the final value of DeltaOD

    '''
    Now we have predicted the OD of YFP.
    Find also the OD of K12 as
    OD_K12(t_sat) = OD_total(t_sat) -OD_YFP(t_sat)
    '''

    Delta_OD_K12 = delta_OD_total[-1]-Delta_OD_YFP

    marker = 'x' if (which_row,which_col) in bubbly_samples else 'o'
    color = cmap(condition)
    color = adjust_lightness(color,alpha)

    if (which_row,which_col) not in bubbly_samples:
        ### Append to output data only if not in a bubbly sample.
        YFP_data.append(Delta_OD_YFP)
        K12_data.append(Delta_OD_K12)
        K12_error.append(error) #which is the same for YFP.
        data_conditions.append((condition,which_row//2))
        saturation_times.append(time[t_sat])

    label_value = label_dict[condition] if condition not in plotted_conditions else None
    if marker == 'o':
        '''Append if not a bubbly sample '''
        color_list.append(color)
        label_value_list.append(label_value)
        if alpha not in alpha_list:
            alpha_list.append(alpha)
    ax.scatter(Delta_OD_YFP,Delta_OD_K12,
               s=50,edgecolor='k',
               marker=marker,
               color=color,
               zorder=2,
               label = label_value)
    ax2.scatter(-1e6,-1e6, #Dummy
               s=50,edgecolor='k',
               marker='o',
               color=adjust_lightness('gray',alpha),
               label = OD_0_theoretical[which_row//2] if alpha not in plotted_dummies else None)
    if alpha not in plotted_dummies: plotted_dummies.append(alpha)
    if condition not in plotted_conditions: plotted_conditions.append(condition)

    ax.errorbar(Delta_OD_YFP,Delta_OD_K12,
                xerr=error,yerr=error,
                ecolor='gray',zorder=1,
                elinewidth=0.5)

''' Output data to directory so that it can be preloaded for other things '''
with open('YFP_data.dat','wb') as f:
    pickle.dump(YFP_data,f)
    f.close()
with open('K12_data.dat','wb') as f:
    pickle.dump(K12_data,f)
    f.close()
with open('K12_error.dat','wb') as f:
    pickle.dump(K12_error,f)
    f.close()
with open('data_conditions.dat','wb') as f:
    pickle.dump(data_conditions,f)
    f.close()
with open('color_list.dat','wb') as f:
    pickle.dump(color_list,f)
    f.close()
with open('label_value_list.dat','wb') as f:
    pickle.dump(label_value_list,f)
    f.close()
with open('alpha_list.dat','wb') as f:
    pickle.dump(alpha_list,f)
    f.close()
with open('saturation_times.dat','wb') as f:
    pickle.dump(saturation_times,f)
    f.close()


####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################

#### Data analysis ends here. The plotting is done in a separate file for simplicity.
#### What follows is a preliminary version of the final figure that includes the bubbly samples.


with open('monoculture_yields.dat','rb') as f: #Monoculture yields.
    K_k12,K_k12_std,K_yfp,K_yfp_std = pickle.load(f)
    f.close()
'''
From monoculture measurements of the YFP yield and K12 yield,
plot the prediction for the co-culture experiment.
'''

##K_yfp = 0.4426232727272727
##K_yfp_std = 0.015
##K_k12 = 0.4170263333333333
##K_k12_std = 0.00629

#Plot predicted curve between k_yfp and K_k12
ax.plot(x_vect,K_k12-K_k12/K_yfp*x_vect,'k--',zorder=1,
        label='Prediction from \n monoculture')
polygon_vertices = np.array([[K_yfp-K_yfp_std,0],
                    [0,K_k12-K_k12_std],
                    [0,K_k12+K_k12_std],
                    [K_yfp+K_yfp_std,0]])
poly_x = polygon_vertices[:,0]
poly_y = polygon_vertices[:,1]
ax.fill(poly_x,poly_y,
        'lightgray',zorder=0)




'''Plot axes and other decorators'''
ax.axhline(y=0, color='gray')
ax.axvline(x=0, color='gray')
ax.set_xlim(-0.1,max(K_yfp,K_k12)*1.2)
ax.set_ylim(-0.1,max(K_yfp,K_k12)*1.2)
ax2.set_yticks([])
ax2.legend(loc='lower left',frameon=True,
           title='Initial \nBiomass', #style 1
           )
ax.legend(loc='upper right',frameon=False,
          title='Initial %YFP', #Style 1
          )
#### Axlabels 1
ax.set_xlabel(r'$\Delta OD$ YFP')
ax.set_ylabel(r'$\Delta OD$ K12')
ax.set_title('Preliminary Figure \n (Includes bubbly replicates)')
format_type = 'png'
plt.savefig('preliminary_final_figure.{}'.format(format_type),
            format=format_type,
            transparent=True)
plt.show()

















    
