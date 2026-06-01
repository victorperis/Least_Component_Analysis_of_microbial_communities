import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy import odr
from scipy.stats import t
import pickle

figsize_factor = 0.77


'''
Choose whether to plot all datapoints, or only those for 
concentrations (written as the corresponding column index)
of the samples in which the OD is proportional to glucose concentration,
from the OD vs C/C_0 plots.

'''
''' For manuscript main '''
good_concentrations = (7,19); output_filename = 'cross-feeding-experiment-results' #largest good inlcuded, smallest good not included.   
''' For SI '''
##good_concentrations =(0,24); output_filename = 'cross-feeding-experiment-results-SI' #if we want to plot everything


'''
Plot in logarithmic scale ?
'''
logarithmic = False
if logarithmic:
    output_filename  = 'logarithmic-'+output_filename


'''
Remove outliers from the analysis. These outliers were singled out during colony counting due to
unusually low colony counts, and are likely to result from tip malfunction during plating,
and the outliers are reproducible across replicate platings (the same tips were used for all replicate platings).
'''
plating_outliers = [(8,0), #Write the well-position in the spreadsheet for EColi. Acid is the same but row +8
                    (8,1),
                    (10,1),
                    (6,2),
                    (4,3),
                    (15,3),
                    (19,1)]
                    
                    


''' We use the poisson estimator from Michael Martini et al. Microbiology Spectrum 2024 '''
''' Initialize arrays of counts and plating_dilutions. '''
EColi_sum_counts = np.zeros((8,24)); EColi_sum_counts[:,:] = np.nan
EColi_sum_plating_dilutions = np.zeros((8,24)); EColi_sum_plating_dilutions[:,:] = np.nan
Acid_sp_sum_counts = np.zeros((8,24)); Acid_sp_sum_counts[:,:] = np.nan
Acid_sp_sum_plating_dilutions = np.zeros((8,24)); Acid_sp_sum_plating_dilutions[:,:] = np.nan

    

cmap = plt.get_cmap('RdBu')


filename = 'R5.xlsx'
plate = pd.read_excel(filename,
                      sheet_name = None,
                      header=None)

dils = ['-6','-5','-4','-3','-2','-1']
markers = ['x','o']
replicate_names = ['R1','R2','R3','R4','R5'] #R1-R3 counted by stella, R4-R5 counted by dani. All recorded by Victor.
working_dilutions = np.arange(1,25)*1.5

for jj,replicate_name in enumerate(replicate_names):
    plate = pd.read_excel('{}.xlsx'.format(replicate_name),
                      sheet_name = None,
                      header=None)
    for ii,dil in enumerate(dils):
        plating_dilution = 10**float(dil)
        plate_current_plating_dilution = plate[dil]
        plate_current_plating_dilution = plate_current_plating_dilution.to_numpy()
        rows,cols = np.shape(plate_current_plating_dilution)
        for row in range(6):
            for col in range(good_concentrations[0],good_concentrations[1]):
                if (col,row) in plating_outliers:
                        ''' Skip Plating Outliers '''
                        continue
                ''' Counts of each species '''
                Ecoli = plate_current_plating_dilution[row,col]
                Acid_sp = plate_current_plating_dilution[row+8,col]
                try:
                    Ecoli = float(Ecoli)
                    Acid_sp = float(Acid_sp)
                    
                except:
                    ''' If not transformable to float, they are words indicating problem '''
                    print('Words present in ',(col,row))
                    continue

                if np.isnan(Ecoli)==False and np.isnan(Acid_sp)==False:
                    ''' If there are counts (not overcrowded), append to estimator '''
                    if np.isnan(EColi_sum_counts[row,col]):
                        ''' Not initialized, so initialize counts '''
                        EColi_sum_counts[row,col] = Ecoli
                        EColi_sum_plating_dilutions[row,col] = plating_dilution
                        Acid_sp_sum_counts[row,col] = Acid_sp
                        Acid_sp_sum_plating_dilutions[row,col] = plating_dilution
                        
                    else:
                        ''' Add counts to current count'''
                        EColi_sum_counts[row,col] += Ecoli
                        EColi_sum_plating_dilutions[row,col] += plating_dilution
                        Acid_sp_sum_counts[row,col] += Acid_sp
                        Acid_sp_sum_plating_dilutions[row,col] += plating_dilution

V = 0.01 #10 uL  the volume in each plating.
EColi_r_mle = EColi_sum_counts/EColi_sum_plating_dilutions/V
EColi_r_sigma = EColi_r_mle/np.sqrt(EColi_sum_counts)
Acid_sp_r_mle = Acid_sp_sum_counts/Acid_sp_sum_plating_dilutions/V
Acid_sp_r_sigma = Acid_sp_r_mle/np.sqrt(Acid_sp_sum_counts)

markers = ['o','x','^'] #one for each biological replicate

full_dataset = [] #for PCA analysis

'''
Open Figure
'''
figsize = np.array([5.0,3.0])
figsize *= figsize_factor
fig,ax = plt.subplots(figsize=figsize)
for row in range(6): #bc row 7 is blank and row 8 is F7 strain (n.a.)
    for col in range(good_concentrations[0],good_concentrations[1]):
        alpha = 1-(col/24)**0.5
        color = cmap(col/24) #one color per glucose concentration
        if (col,row) in plating_outliers:
            ''' Do not append plating outliers to dataset '''
            continue
        if np.isnan(EColi_r_mle[row,col])==True or np.isnan(Acid_sp_r_mle[row,col])==True:
            print(col,row,'has nans')
        full_dataset.append((EColi_r_mle[row,col],
                      Acid_sp_r_mle[row,col]))
        ax.errorbar(EColi_r_mle[row,col],
                    Acid_sp_r_mle[row,col],
                    xerr=EColi_r_sigma[row,col],
                    yerr=Acid_sp_r_sigma[row,col],
                    linewidth=0,elinewidth=2,
                    color='darkgray',zorder=1,alpha=0.5
                    )
        ax.scatter(EColi_r_mle[row,col],
                    Acid_sp_r_mle[row,col],
                   s=30*figsize_factor,zorder=2,
                   c=color,#marker=marker,
                    alpha=1
                    )

'''
Perform PCA on the full dataset
'''
full_dataset = np.array(full_dataset)
strain_stds = np.std(full_dataset,axis=0,ddof=1)
strain_means = np.mean(full_dataset,axis=0)
z_scored = True
print('Using z_scored abundances in PCA:',z_scored)
if z_scored:
    corr = np.corrcoef(full_dataset.T)
else:
    corr = np.cov(full_dataset.T)
eigvals,eigvects = np.linalg.eigh(corr)
argsort = np.argsort(eigvals)
eigvals = eigvals[argsort]; eigvects = eigvects[:,argsort]

if z_scored:
    eigvects = eigvects/strain_stds[:,None]#rescale according to change from z to n

eigvects = eigvects/np.sqrt(np.sum(eigvects**2,axis=0))[None,:] #rescale vector norm to 1
print('Y_ac/Y_EColi on acid: {}'.format(eigvects[0,0]/eigvects[1,0]))

if z_scored == False:
    eigvects /= 1e8 #else the arrows are huge

arrow_length = np.std(eigvects[:,0]@(full_dataset).T)
print('arrow_length',arrow_length)

if 'SI' not in output_filename:
    ''' Plot arrows for  '''
    ### Least Component
    ax.arrow(
             1.25e8,0.25e9,
##             0.4*1e9*eigvects[0,0],0.4*1e9*eigvects[1,0], #scale arrow length to fit it nicely, but is arbitrary.
             arrow_length*eigvects[0,0],arrow_length*eigvects[1,0], #length is proportional to spread of data.
             facecolor='tab:purple',edgecolor='k',
             width=4e7)

'''
Try plotting the linear relationship from cross-feeding, using
the slope derived from the least component.
'''
slope = -eigvects[0,0]/eigvects[1,0]
x_vect = np.linspace(0,1e9,100)
if logarithmic==False and 'SI' not in output_filename:
    ax.plot(x_vect,x_vect*slope,
            'gray',linestyle=':',
            linewidth=3)


'''
Calculate mean and std for each experimental condition (glucose conc.)
across the 6 replicates, and plot with a square marker.
'''
n0_EColi = 1e100
n0_Acid = 1e100
for col in range(good_concentrations[0],good_concentrations[1]):
    ''' Calculate mean and std. of estimates in a given column.
    Aka.: of replicates in the same conditions (glucose concentration).
    '''
    color = cmap(col/24)
    mean_Ecoli = np.nanmean(EColi_r_mle[:,col])
    mean_Acid = np.nanmean(Acid_sp_r_mle[:,col])
    if mean_Ecoli<n0_EColi: n0_EColi = mean_Ecoli
    if mean_Acid<n0_Acid: n0_Acid = mean_Acid
    std_Ecoli = np.nanstd(EColi_r_mle[:,col])
    std_Acid = np.nanstd(Acid_sp_r_mle[:,col])
    ax.errorbar(mean_Ecoli,
                mean_Acid,
                xerr = std_Ecoli,
                yerr = std_Acid,
                c='k',zorder=5,
                elinewidth=2)
    ax.scatter(mean_Ecoli,mean_Acid,
               color=color,
               s=80*figsize_factor,marker='s',edgecolor='k',zorder=10)

''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
'''
Add null expectation for exponentially growing bacteria under competition for glucose
'''
n0_EColi = 5174488.567990371 #mean across replicates for column 24, the lowest glucose abundance. good indicator of initial CFU/mL since there is virtually no growth
n0_Acid = 5090252.707581226 #idem

### Growth rates and yields, measured in monoculture
with open('../lab_monocultures_find_acid_specialist/Acid_specialist_glucose_yields.dat','rb') as f:
    _ = pickle.load(f)
    f.close()
    Acid_sp_glucose_yield,Acid_sp_glucose_yield_err = _
with open('../lab_monocultures_find_acid_specialist/Acid_specialist_glucose_growth_rate.dat','rb') as f:
    _ = pickle.load(f)
    f.close()
    Acid_sp_growth_rate,Acid_sp_growth_rate_err = _
with open('../competition_data_analysis/monoculture_yields.dat','rb') as f:
    _ = pickle.load(f)
    f.close()
    EColi_glucose_yield,EColi_glucose_yield_err = _[:2] #the two FIRST elements are those of K12
with open('../competition_data_analysis/monoculture_growth_rates.dat','rb') as f:
    _ = pickle.load(f)
    f.close()
    EColi_growth_rate,EColi_growth_rate_err = _[-2:] #the two LAST elements are those of K12
##Acid_sp_growth_rate = 0.037 #in h^-1. From monoculture growth in glucose.
##Acid_sp_growth_rate_err = 0.002
##EColi_growth_rate = 0.151  #in h^-1.
##EColi_growth_rate_err = 0.009


''' We can directly write Acid_sp(EColi) '''
EColi_vect = np.linspace(n0_EColi,1e9,1000)
Acid_sp_from_EColi_in_competition = (EColi_vect/EColi_vect[0])**(Acid_sp_growth_rate/EColi_growth_rate)*n0_Acid
Acid_sp_from_EColi_in_competition_upper = (EColi_vect/EColi_vect[0])**(Acid_sp_growth_rate/EColi_growth_rate/3)*n0_Acid
Acid_sp_from_EColi_in_competition_lower = (EColi_vect/EColi_vect[0])**(Acid_sp_growth_rate/EColi_growth_rate*3)*n0_Acid

if logarithmic==False and 'SI' not in output_filename:
    ax.plot(EColi_vect,
            Acid_sp_from_EColi_in_competition,
            '--',color='gray',linewidth=3)






''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

'''
Add Colorbar
'''
### Add colorbar ###
import matplotlib as mpl
import matplotlib.ticker as mticker

norm = mpl.colors.Normalize(
    vmax=np.log10(3/1.5**0),
    vmin=np.log10(3/1.5**24)
)
conc_max = 3/1.5**good_concentrations[0]
conc_min = 3/1.5**(good_concentrations[1]-1)
print(conc_max,conc_min)

sm = mpl.cm.ScalarMappable(norm=norm, cmap='RdBu_r')
sm.set_array([])  # required for older matplotlib versions

cbar = plt.colorbar(sm, ax=ax)
cbar.set_label(r'Glucose Conc. (% weight)',fontsize=8)
cbar.set_ticks([0,-1,-2,-3])
cbar.set_ticklabels([r'$10^{0}$',r'$10^{-1}$',r'$10^{-2}$',r'$10^{-3}$'],fontsize=8)
ax.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
ax.yaxis.get_major_formatter().set_useMathText(True)
ax.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))
ax.xaxis.get_major_formatter().set_useMathText(True)
if logarithmic:
    ax.set_xscale('log')
    ax.set_yscale('log')
if logarithmic ==False and 'SI' not in output_filename:
    ax.set_xticks(np.array([0,0.25,0.5,0.75])*1e9)
    ax.set_aspect(strain_stds[1]/strain_stds[0]*0.7)
ax.yaxis.get_offset_text().set_fontsize(8)
ax.xaxis.get_offset_text().set_fontsize(8)
ax.set_xlabel(r'$E_K$'+' Strain \n Abundance ('+r'$\times 10^8$ CFU/mL)',fontsize=8)
ax.set_ylabel('A. Guillouiae \n Abundance ('+r'$\times 10^9$ CFU/mL)',fontsize=8)
ax.set_ylim(top=1.25e9)
plt.tight_layout()
plt.savefig('{}.svg'.format(output_filename),format='svg',transparent=True)
plt.savefig('{}.png'.format(output_filename),format='png',dpi=300,transparent=True)


''' Plot Barplot of glucose growth rates of both strains '''
figsize_factor = 0.5
figsize = np.array([1,1.5])
figsize*= figsize_factor
fig,ax = plt.subplots(figsize=figsize)
for i in range(2):
    if i==0: facecolor = 'peru';data = EColi_growth_rate; error=EColi_growth_rate_err
    elif i==1: facecolor='lightcoral';data = Acid_sp_growth_rate; error=Acid_sp_growth_rate_err
    ax.bar(i,
           data,
           facecolor=facecolor,edgecolor='k',
           width=1,linewidth=1,
        )
    ax.errorbar(i,
           data,
           yerr=error,
           ecolor='k',elinewidth=2,
        )

ax.set_xticks([0,1])
ax.set_yticklabels(ax.get_yticks(),fontsize=8)
ax.set_xticklabels(['K12','Ac.Sp.'],fontsize=8)
ax.set_ylabel(r'Growth Rate $(h^{-1})$', fontsize=8)
plt.tight_layout()
plt.savefig('barplot-growth-rates.svg',format='svg',transparent=True)
plt.savefig('barplot-growth-rates.png',format='png',dpi=300,transparent=True)


''' Plot Barplot of glucose growth rates of both strains '''
figsize_factor = 0.5
figsize = np.array([1,1.5])
figsize*= figsize_factor
fig,ax = plt.subplots(figsize=figsize)
gluc_conc = 0.2 #in % weight
for i in range(2):
    if i==0: facecolor = 'peru';data = EColi_glucose_yield/gluc_conc; error=EColi_glucose_yield_err/gluc_conc
    elif i==1: facecolor='lightcoral';data = Acid_sp_glucose_yield/gluc_conc; error=Acid_sp_glucose_yield_err/gluc_conc
    ax.bar(i,
           data,
           facecolor=facecolor,edgecolor='k',
           width=1,linewidth=1,
        )
    ax.errorbar(i,
           data,
           yerr=error,
           ecolor='k',elinewidth=2,
        )

ax.set_xticks([0,1])
ax.set_yticklabels(ax.get_yticks(),fontsize=8)
ax.set_xticklabels(['K12','Ac.Sp.'],fontsize=8)
ax.set_ylabel(r'Yield $(OD/[Gluc.])$',fontsize=8)
plt.tight_layout()
plt.savefig('barplot-yields.svg',format='svg',transparent=True)
plt.savefig('barplot-yields.png',format='png',dpi=300,transparent=True)






plt.show()
