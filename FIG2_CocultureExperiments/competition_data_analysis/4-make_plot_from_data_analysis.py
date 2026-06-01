import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.ndimage import gaussian_filter1d
from scipy.optimize import curve_fit
from scipy.integrate import cumulative_trapezoid
import pickle

figsize_factor = 0.67

OD_0 = 0.15
OD_0_theoretical = np.round(np.array([1,1/3,1/9,1/27])*OD_0,3)

def adjust_lightness(color, amount=0.5):
    import matplotlib.colors as mc
    import colorsys
    try:
        c = mc.cnames[color]
    except:
        c = color
    c = colorsys.rgb_to_hls(*mc.to_rgb(c))
    return colorsys.hls_to_rgb(c[0], max(0, min(1, amount * c[1])), c[2])


''' Output data to directory so that it can be preloaded for other things '''
with open('YFP_data.dat','rb') as f:
    YFP_data = pickle.load(f)
    f.close()
with open('K12_data.dat','rb') as f:
    K12_data= pickle.load(f)
    f.close()
with open('K12_error.dat','rb') as f:
    K12_error = pickle.load(f)
    f.close()
with open('data_conditions.dat','rb') as f:
    data_conditions = pickle.load(f)
    f.close()
with open('color_list.dat','rb') as f:
    color_list = pickle.load(f)
    f.close()
with open('label_value_list.dat','rb') as f:
    label_value_list = pickle.load(f)
    f.close()
with open('alpha_list.dat','rb') as f:
    alpha_list = pickle.load(f)
    f.close()
with open('monoculture_yields.dat','rb') as f: #Monoculture yields.
    K_k12,K_k12_std,K_yfp,K_yfp_std = pickle.load(f)
    f.close()

figsize = np.array([4,3],dtype=float)
figsize *= figsize_factor
fig, ax = plt.subplots(figsize=figsize)



# Plot experimental points from co-culture
ax2 = ax.twinx()
plotted_dummies = []
for i in range(len(K12_data)):

    condition,which_row_half = data_conditions[i]
    which_row = which_row_half*2
    
    
    ax.scatter(YFP_data[i],K12_data[i],
               c = color_list[i])

    ax.scatter(YFP_data[i],K12_data[i],
               s=50*figsize_factor,
               edgecolor='k', #add for old version
               color=color_list[i],
               zorder=2,
               label = label_value_list[i])
    ax.errorbar(YFP_data[i],K12_data[i],
                xerr=K12_error[i],yerr=K12_error[i],
                ecolor='gray',zorder=1, #old version
                elinewidth=0.5,
##                ecolor='lightgray',zorder=1, #new version, if plotting avg across replicates
##                elinewidth=2
                )
    
    alpha = 1+(which_row//2)/5

    label = OD_0_theoretical[which_row//2] if OD_0_theoretical[which_row//2] not in plotted_dummies else None
    ax2.scatter(-1e6,-1e6, #Dummy
               s=50*figsize_factor,edgecolor='k',
               marker='o',
               color=adjust_lightness('gray',alpha),
               label = label)

    if label not in plotted_dummies and label!= None: plotted_dummies.append(label)


'''
Plot averages per condition, as in the cross-feeding experiment.
'''
##K12_data = np.array(K12_data)
##YFP_data = np.array(YFP_data)
##plotted_conditions = []
##for i in range(len(K12_data)):
##
##    condition,which_row_half = data_conditions[i]
##    which_row = which_row_half*2
##
##
##    if (condition,which_row_half) not in plotted_conditions:
##        plotted_conditions.append((condition,which_row_half))
##
##        where = [(condition,which_row_half)==i for i in data_conditions]
##
##        avg_K12 = np.mean(K12_data[where])
##        avg_YFP = np.mean(YFP_data[where])
##        std_K12 = np.std(K12_data[where],ddof=1)
##        std_YFP = np.std(YFP_data[where],ddof=1)
##        
##
##        ax.scatter(avg_YFP,avg_K12,
##               s=50,edgecolor='k',marker='s',
##               color=color_list[i],
##               zorder=2,
##               )
##        ax.errorbar(avg_YFP,avg_K12,
##                    xerr=std_YFP,yerr=std_K12,
##                    ecolor='gray',zorder=1,
##                    elinewidth=2)

        
        

'''
Perform analysis of the correlation matrix and plot vector arrows
'''
YFP_data = np.array(YFP_data)
K12_data = np.array(K12_data)
complete_data = np.array((YFP_data,K12_data))
data_std = np.std(complete_data,axis=1,ddof=1)
print(data_std)
z_scored = True
print('Using z_scored abundances in PCA:',z_scored)

if z_scored:
    corr = np.corrcoef(complete_data)
else:
    corr = np.cov(complete_data)
eigvals,eigvects = np.linalg.eig(corr)
sorted_indices = np.argsort(eigvals)
eigvals = eigvals[sorted_indices]
eigvects = eigvects[:,sorted_indices]
print('Eigvals',eigvals)


for e,eigvect in enumerate(eigvects.T):
    if z_scored:
        eigvect = eigvect/data_std #if using the correlation matrix, rescale eigenvector to work in n, not z-scores
        eigvect_norm = np.sqrt(np.sum(eigvect**2))
        eigvect = eigvect/eigvect_norm #normalize to norm 1 for the projection.
    else:
        eigvect *= 10
        if e==0:
            eigvect *= -1 #switch direction of lowvar vector if using covariance, to match
                            #direction of inverse yields
        pass
        
    if e==1:
        color='tab:green'
        x0 = 0.12; y0 = 0.38
    else:
        color='tab:purple'
##        eigvect = -eigvect
        x0 = 0.425; y0 = 0.07

        factor = 120 #reduce arrow size so that it fits
        ax.arrow(0.4,0.09,
                 1/K_yfp/factor,1/K_k12/factor,
                 facecolor='k',
                 width = 0.008,zorder=0
                 )

    ### Calculate dx,dy for the arrows from the data itself.
    length = np.sqrt(np.var(eigvect@complete_data))*5
    
    dx,dy = eigvect*length#*2.75 #if *2.75, make it the same length as the black arrow just for visuals.


    if e==0:
        ''' Print ratio K_YFP/K_K12 from the eigenvector '''
        estimated_ratio = (eigvect)
        estimated_ratio = eigvect[1]/eigvect[0] #results in \hat{K}_yfp/\hat{K}_K12
        print('Estimated ratio: {}'.format(estimated_ratio))

    if e==0: #plot only Least Component
        ax.arrow(x0,y0,dx,dy,
                 facecolor=color,
                 length_includes_head = True, #else it's misleadingly large
                 head_length = 0.8*length,
                 width=0.008,
                 edgecolor='k')



    

##Plot predicted curve between k_yfp and K_k12
x_vect  =np.linspace(0,max(K_yfp,K_k12),1000)
ax.plot(x_vect,K_k12-K_k12/K_yfp*x_vect,'gray',linestyle='--',zorder=1,
        label='Prediction from \n monoculture \n yields',
        linewidth=3)



'''Plot axes and other decorators'''
ax.set_aspect(data_std[0]/data_std[1])
ax.set_xlim(-0.1,max(K_yfp,K_k12)*1.2)
ax.set_ylim(-0.3,max(K_yfp,K_k12)*1.4)
ax2.set_yticks([])
legend = ax2.legend(loc='lower left',frameon=True,
           title='Initial \nBiomass',fontsize=8,#style 1
           )
plt.setp(legend.get_title(),fontsize=8)
legend= ax.legend(loc='upper right',frameon=False,
          title='Initial %YFP Strain',fontsize=8,#Style 1
          )
plt.setp(legend.get_title(),fontsize=8)

ax.axhline(y=0,zorder=0,color='lightgray')#,alpha=0.5)
ax.axvline(x=0,zorder=0,color='lightgray')#,alpha=0.5)
#### Axlabels 2
ax.set_xlabel(r'YFP Strain Abundance ($OD_f-OD_i$)')
ax.set_ylabel(r'K12 Strain Abundance ($OD_f-OD_i$)')
ax.set_xticks([0,0.2,0.4])
ax.set_xticklabels(ax.get_xticks(),fontsize=8)
ax.set_yticks([0,0.2])
ax.set_yticklabels(ax.get_yticks(),fontsize=8)
format_type = 'svg'

plt.tight_layout()
plt.savefig('final_figure.{}'.format(format_type),
            format=format_type,
##            dpi=300,
            transparent=True)
plt.savefig('final_figure.png',
            format='png',
##            dpi=300,
            transparent=True)
##plt.show()


''' Find the estimations for K_k12 and K_yfp from co-culture '''
print('K_monoculture',
      'YFP',
      K_yfp,'+-',K_yfp_std, '\n',
      'K12',
      K_k12,'+-',K_k12_std,'\n',
      )

def true_ratio(K_1,K_2,err_1,err_2):
    out = K_1/K_2
    error = err_1/K_2+K_1/K_2**2*err_2 
    return out,error

true_ratio_val,ratio_error = true_ratio(K_yfp,K_k12,K_yfp_std,K_k12_std)
print('Ratio from monoculture: {}+-{}'.format(true_ratio_val,ratio_error))




plt.show()














    
