import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import itertools
import pickle


rows = ['A','B','C','D','E','F','G','H']
cols = np.array(np.arange(1,13),dtype=str)

tags= np.array(list(itertools.product(cols,rows)),dtype=str)
for i in range(len(tags)):
    tag = tags[i]
    tags[i] = tag[1]+tag[0]

tags = tags[:,0]

print(tags)



fig,ax = plt.subplots(figsize=(4,3))


mediums = ['glu-0.2','cit-0.2']
F2_data = []; is_F2_plotted = False
for replicate in ['R1','R2','R3']:
##    if replicate !='R1': continue
    
    if replicate=='R1':
        days=['17','18','19']
        month = '06'
        year = '25'
        mediums = ['glu-0.2','cit-0.2']
    elif replicate=='R2':
        days= ['01','02','03']
        month = '07'
        year = '2025'
        mediums = ['Gluc0.2.VictorStella.r2','Citr0.2.VictorStella.r2']
    elif replicate=='R3':
        days = ['02','03','04']
        month = '07'
        year = '25'
        mediums = ['glu0.2-r3-sp-vpy','cit0.2-r3-sp-vpy']

    
    for day in days:
        ### On every day of the experiment
        filenames = []
        for medium in mediums:
            ### Plot OD_citrate vs OD_glucose
            #### We used a different format for the names of the files of each replicate...
            if replicate=='R1':
                filename = 'R1/{}-{}-{}-{}.xlsx'.format(day,month,year,medium)
            elif replicate=='R2':
                filename = 'R2/{}.{}.{}.{}.xlsx'.format(year,month,day,medium)
            elif replicate=='R3':
                filename = 'R3/{}-{}-{}-{}.xlsx'.format(day,month,year,medium)
            filenames.append(filename)

        print(filenames)
        outs = []
        for i,filename in enumerate(filenames):

        ##    filename = '{}-{}-{}-{}.xlsx'.format(day,month,year,medium)

            data = pd.read_excel(filename) #data organized in rows, then cols: A1,A2,...,A12,B1,...


            ''' Convert data to array and reshape with plate shape '''
            data = data.iloc[0,1:]
            data = np.array(data,dtype=float)

            outs.append(data)


        outs = np.array(outs)

        ### Focus on the F2 strain, which is the acid specialist        
        better_citrate = np.array(tags)=='F2'
        chosen_tags = tags[better_citrate]
        x_place = outs[0,:][better_citrate]*1.01
        y_place = outs[1,:][better_citrate]*1.01
        F2_data.append(outs[:,better_citrate])
        
        ax.scatter(outs[0,:],outs[1,:],#label='{}/{}'.format(day,month),
                   s=50,facecolor='none',edgecolor='lightgray')
        
        ax.scatter(outs[0,better_citrate],outs[1,better_citrate],#label='{}/{}'.format(day,month),
                   s=50,facecolor='none',edgecolor='crimson',
                   label='Individual Replicate' if is_F2_plotted==False else None)
        is_F2_plotted = True

F2_mean = np.mean(F2_data,axis=0)
F2_error = np.std(F2_data,axis=0,ddof=1)
ax.errorbar(F2_mean[0],F2_mean[1],
            xerr= F2_error[0],
            yerr= F2_error[1],
            color='gray',)
ax.scatter(F2_mean[0],F2_mean[1],
           s=80,zorder=10,
           marker='s',edgecolor='k',facecolor='tab:red',
           label = 'Mean')



with open('Acid_specialist_glucose_yields.dat','wb') as f:
    pickle.dump((F2_mean[0]-0.08, #subtract control OD
                 F2_error[0]),
                f)
    f.close()

line = np.linspace(0,1.4,100)
ax.plot(line,line,'k--')
ax.set_xlim(0,1.4)
ax.set_ylim(0,1.4)
ax.set_xlabel('OD (Glucose 0.2%)')
ax.set_ylabel('OD (Citrate 0.2%)')
ax.legend(loc='best',frameon=False,
##          title='Day'
          title = 'Acinetobacter guillouiae'
          )
##ax.set_title('{}-{}-{}'.format(day,month,year))
plt.tight_layout()
plt.savefig('F2_growth_on_citrate_and_glucose.png',format='png',dpi=300,transparent=True)
plt.savefig('F2_growth_on_citrate_and_glucose.svg',format='svg',transparent=True)
plt.show()
