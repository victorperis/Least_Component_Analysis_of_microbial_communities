import numpy as np
import matplotlib.pyplot as plt
import pickle
import pandas as pd
from scipy.ndimage import gaussian_filter1d
from scipy.optimize import curve_fit
import matplotlib.cm as cm

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



conditions = 6
cmap = cm.get_cmap('viridis',conditions) #creates a cmap with as many colors as conditions
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
t0 = 7 #index from which to start plotting.
ODs = data.iloc[3:86,2:]
ODs = ODs.to_numpy() #Rows are timepoints, cols are positions in the well

which_controls_are_bad = [0,48,59,71,95]
which_controls_are_good = [11,12,
                           23,24,
                           35,36,
                           47,
                           60,
                           72,
                           83,84]


controls = np.vstack((ODs[:,which_controls_are_good]))
control_means_in_time = np.mean(controls,axis=1) #avg of both controls at each point in time.

plotted_conditions = []
growth_rates = []
fig,ax = plt.subplots()
fig2,ax2 = plt.subplots(figsize=(4,3))
chosen_conditions = np.array([0,5]) #100% YFP and 100% K12
strain_names = ['YFP','K12']
growth_rates_by_strain = [[],[]]
for i in range(np.shape(ODs)[1]):
    plot_row = i//cols
    plot_col = i%cols

    '''
    Skip controls
    '''
    if i%cols ==0 or i%cols==11:
        continue

    

    OD = ODs[:,i]-control_means_in_time[:]
    OD = OD-OD[0] #remove initial to show \Delta OD

    ''' Find kneepoint '''
    OD = np.array(OD,dtype=float)
    OD_filtered = gaussian_filter1d(OD,sigma=4) #used only to get the kneepoint.
    where_kneepoint = find_kneepoint(OD_filtered,time)

    condition = which_condition(i)
    alpha = 1+(plot_row//2)/5
    color = cmap(condition)
    color = adjust_lightness(color,alpha)
    
    

    if condition not in chosen_conditions:
        continue
    if condition==0 and plot_row<3:
        ''' We have a bubbly sample of YFP, so we do not consider it for analysis. '''
        continue
    print(plot_row,plot_col,condition)
    which_chosen_condition = np.where(condition==chosen_conditions)[0][0]
    which_strain = strain_names[which_chosen_condition]

    '''
    Choose some OD interval [OD_min,OD_max], both of which are low, in
    which exponential growth can be assumed and find the exponential growth
    rate from those curves.
    '''
    if True:

        if which_strain == 'K12':
            OD_min = (0.1)# if np.min(OD)<0.01 else np.min(OD))
            OD_max = (0.3 if np.max(OD)>0.3 else np.max(OD))
        elif which_strain=='YFP':
            OD_min = (0.12)
            OD_max = (0.3 if np.max(OD)>0.3 else np.max(OD))

        '''
        For K12, skip the last 2 rows where it does not saturate.
        '''
        if which_strain=='K12' and plot_row in [6,7]: continue

        t_OD_min = np.where(OD<OD_min)[0][-1] #choose last time it passe a given OD, which will be good even for bubbbly samples.
        t_OD_max = np.where(OD>0.99*OD_max)[0][0]
        t_OD_max = where_kneepoint

        OD_difference_in_interval = np.log(OD[t_OD_max])-np.log(OD[t_OD_min])
        time_difference_in_interval = time[t_OD_max]-time[t_OD_min]

        OD_interval = OD[t_OD_min:t_OD_max]
        time_interval = time[t_OD_min:t_OD_max]

        '''
        Fit a growth rate through scipy curve fit
        '''
        popt,pcov,infodict,mesg,ier = curve_fit(exponential_growth,time_interval,OD_interval,
                                       full_output=True)
        r,x0 = popt
        print('r',r)
        growth_rates.append(r)
        growth_rates_by_strain[which_chosen_condition].append(r)
        exponential_error = np.mean(infodict['fvec'])

        exponential_curve = exponential_growth(time_interval,r,x0)


        ax2.scatter(which_chosen_condition, #plot each column in a different x-value
                    r,
                    s=50,edgecolor='k',
                    color=cmap(condition), #color according to strain only
                    label = which_strain if condition not in plotted_conditions else None)
        
    
        ax.plot(
##            time_interval,OD_interval,
                time[t_OD_min:],OD[t_OD_min:],
                                   linewidth=3,color=color if condition !=-1 else 'k',
                                label = ('YFP' if condition==0 else 'K12') if condition not in plotted_conditions else None,
##                label = ('Sp. 1' if condition==0 else 'Sp. 2') if condition not in plotted_conditions else None
                )
        

        ax.plot(time_interval,
                                   exponential_curve,
                                   '--',linewidth=2,
                color = cmap(condition))

        ax.set_yscale('log')
        if condition not in plotted_conditions: plotted_conditions.append(condition)

    else:
        pass
    
for i,which_strain in enumerate(strain_names):
    mean = np.mean(growth_rates_by_strain[i])
    std = np.std(growth_rates_by_strain[i])

    ax2.errorbar(i+0.1,mean,
            yerr=std,
            elinewidth=2,
            color=cmap(i if i==0 else 5),
            capsize=5)
    ax2.scatter(i+0.1,mean,c=cmap(i if i==0 else 5),s=50)

    print('{} has growth rate: \n'.format(which_strain),mean,'+-',std)
        
##ax[0,0].set_yticks([0,1])
##ax[0,0].set_ylabel('OD')
##ax[-1,0].set_xticks([time[0],time[-1]])
##ax[-1,0].set_xlabel('Time (h)')
ax.set_xlabel('Time (h)')
##ax.set_ylabel('OD(t)')
ax.set_ylabel('Biomass(t)')
ax.legend(title='Strain',frameon=False)
fig.tight_layout()
fig.savefig('exponential_growth_rates.png',format='png',dpi=600)
##plt.suptitle(filename)




##ax2.boxplot(growth_rates,positions=[-1],zorder=0)
ax2.legend(bbox_to_anchor=(1.01,1),frameon=False,title='Strain')
ax2.set_xticks([])
ax2.set_xlim(-0.2,1.5+0.2)
##ax2.set_xlabel(r'$OD_0$')
ax2.set_ylabel('Growth Rate '+r'$h^{-1}$')
fig2.tight_layout()
plt.savefig('monoculture-strain-growth-rates.svg',format='svg',transparent=True)
plt.savefig('monoculture-strain-growth-rates.png',format='png',transparent=True)



with open('monoculture_growth_rates.dat','wb') as f:
    to_save = (np.mean(growth_rates_by_strain[0]), #yfp
               np.std(growth_rates_by_strain[0]), #yfp
               np.mean(growth_rates_by_strain[1]), #k12
               np.std(growth_rates_by_strain[1]) #k12
               )
    pickle.dump(to_save,f)
    f.close()


fig,ax = plt.subplots(figsize=(1.5,1.5))
for i in range(2):
    if i==0: facecolor = 'khaki';
    elif i==1: facecolor='peru';
    ax.bar(i,
           np.mean(growth_rates_by_strain[i]),
           facecolor=facecolor,edgecolor='k',
           width=1,linewidth=1,
        )
    ax.errorbar(i,
           np.mean(growth_rates_by_strain[i]),
           yerr=np.std(growth_rates_by_strain[i]),
           ecolor='k',elinewidth=2,
        )

ax.set_xticks([0,1])
ax.set_xticklabels(['YFP','K12'])
ax.set_ylabel(r'Growth Rate $(h^{-1})$')
plt.tight_layout()
plt.savefig('barplot-growth-rates.svg',format='svg',transparent=True)
plt.savefig('barplot-growth-rates.png',format='png',dpi=300,transparent=True)
plt.show()


