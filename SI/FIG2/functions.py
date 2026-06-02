import numpy as np
import pickle
from numba import njit
import multiprocessing as mp
import time
import os

v_sp = 1e-6 #variance of the noise on species dynamics

os.environ["OMP_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"

print("OMP_NUM_THREADS =", os.environ.get("OMP_NUM_THREADS"))

@njit
def species_growth_rates(R_vec,n_vec,
                 growth_choices #Shape (N_R,N_n) matrix of 1s and 0s indicating resources being consumed. At most one resource is consumed per species.
                ):
    ''' Growth of the species as a function of the resource abundances: g_i^\alpha(\vec{R},\vec{n})
    By default, we use nR/(1+R) on all resources.
    '''
    R_part = R_vec/(1+R_vec)
    outer = np.outer(R_part,n_vec) #Shape (N_R,N_n)
    outer = outer*growth_choices
    
    # if np.random.random()<1e-3:
    #     print('growth_rates',np.sum(outer>1e-3))
    #     print(outer[0,:])

    ''' Calculate also the total species growth rates '''
    species_total_growth_rate = np.sum(outer,axis=0) #shape (N_species)
    return outer,species_total_growth_rate

@njit
def resource_depletion_production(R_vec,n_vec,
                   total_species_growth_rate, #Shape (N_species)
                                  inv_yields):

    '''
    Calculates g_i^alpha*1/Y_i^\alpha
    '''
##    N_res,N_sp = np.shape(inv_yields)
##    total_species_growth_rate_array = np.repeat(total_species_growth_rate,N_res).reshape((N_res,N_sp))

    total_species_growth_rate_array = np.zeros((len(R_vec),len(n_vec)),dtype=total_species_growth_rate.dtype)
    for r in range(len(R_vec)):
        total_species_growth_rate_array[r,:] = total_species_growth_rate
        
    resource_changes = total_species_growth_rate_array*inv_yields #is shape (N_res,N_species)

    # if np.random.random()<1e-3:
    #     print('production terms',resource_changes)
    
    return resource_changes

@njit
def d_resources(R_vec,n_vec,
               entry_rates,
                growth_rates,
                inv_yields):

    res_changes = resource_depletion_production(R_vec,n_vec,growth_rates,inv_yields)

    # if np.random.random()<1e-3:
    #     print('res_changes',res_changes)

    consumption_production = -np.sum(res_changes,axis=1) #Shape N_resource

    return entry_rates+consumption_production

@njit
def d_species(n_vec,growth_rates):

    species_growth = np.sum(growth_rates,axis=0)

    return species_growth


def main_simulation(args):
    print('PID',os.getpid())
    

    sim,other_arguments = args

    (how_many_simulations,N_species,N_resources,independent_species_cutoff,nsteps,
     growth_choices,entry_rates,inv_yields,
     dilution,dt,sigma) = other_arguments

    
    ''' Resource entries. Supplied resources are always the same, but their entry rates are random. '''
    entered_res = np.arange(min(independent_species_cutoff,N_resources))[entry_rates[:min(independent_species_cutoff,N_resources)]>0] #the ones that marked with 1.
    #### Constant entry rate
    entry_rates[entered_res] = 1 
    print('entry_rates',
          entry_rates)


    if independent_species_cutoff<N_resources:
        entry_rates[independent_species_cutoff:] = np.abs(np.random.normal(entry_rate_mean,entry_rate_sigma,size=N_resources-independent_species_cutoff))
        
    print('enterd res.',entered_res)
    print('PID',os.getpid(),entry_rates[0])


    species_abus = np.random.random(N_species)
    resource_abus = np.random.random(N_resources)
    

    time0 = time.time()


    ''' Prepare colored noise ''' 
    sigma_outside = sigma #use 1 for the simulations as in the paper CR chemostat.
    sigma_inside = sigma
    thres = 1e-12
    outside_resource_noise_prev = np.random.normal(0,sigma_outside,size=(N_resources))*np.sqrt(dt)*(entry_rates>=thres) #noise for supplied resources is large
    inside_resource_noise_prev = np.random.normal(0,sigma_inside,size=(N_resources))*np.sqrt(dt)*(entry_rates<thres) #noise for cross-fed resources is small
    tau = 0.000000001#autocorrelation time of the resource noise    


    for step in range(nsteps):   
    
        
        ''' Growth rates '''
        growth_rates,total_growth_rate_per_species = species_growth_rates(resource_abus,species_abus,
                     growth_choices #Shape (N_R,N_n) matrix of 1s and 0s indicating resources being consumed. At most one resource is consumed per species.
                    )

                
        d_res = d_resources(resource_abus,species_abus,
                   entry_rates,
                    total_growth_rate_per_species,
                   inv_yields)-resource_abus*dilution
    
        
        d_sp = d_species(species_abus,growth_rates)-species_abus*dilution+1e-12
        
        species_abus += d_sp*dt+np.random.normal(0,v_sp**0.5,size=(N_species))*np.sqrt(dt)

        ''' Get white noise '''
        outside_resource_noise = np.random.normal(0,sigma_outside,size=(N_resources))*np.sqrt(dt)*(entry_rates>=thres)
        inside_resource_noise = np.random.normal(0,sigma_inside,size=(N_resources))*np.sqrt(dt)*(entry_rates<thres)

        outside_resource_noise[resource_abus<thres] = 0
        inside_resource_noise[resource_abus<thres] = 0
        
        resource_abus += d_res*dt+outside_resource_noise+inside_resource_noise
        resource_abus[resource_abus<thres] = 0
        species_abus[species_abus<thres] = 0

        outside_resource_noise_prev = outside_resource_noise
        inside_resource_noise_prev = inside_resource_noise

    print('Sim {} done. Took {:.3f} s'.format(sim,time.time()-time0),flush=True)
    
    return sim,resource_abus,species_abus


'''
Initialize all variables as in the jupyter notebook
'''
''' System constants '''
N_species = 20
N_resources = 1
independent_species_cutoff = N_species #N_species #species beyond this cutoff exist independently of all others and are akin to noise
dt = 0.005
dilution = 0.1
how_many_supplied_resources = 1

set_common_rand_seed = False


dataset_name = 'other_r20_logsigma'


nat_scale = N_species*v_sp



if __name__=='__main__':

    how_many_simulations = 100 #500 for lienar sigmas, 100 for log sigmas
    nsteps=20000 #was using 20 000 for linear sigmas, could not be enough to get rid of transients.
    inv_yields = np.ones(N_species)
    growth_choices = np.ones(N_species)


    ''' Resource entries. Selects only which resources are supplied from the outside. But actual entry rates are random in each process. '''
    entry_rates = np.zeros(N_resources)
    entered_res = np.array([0]); entry_rates[0] = 1
    
    entry_rates[entered_res] = 1
        
    ''' Keep the name "traces" for compatibility with option A, even if these is cross-sample data '''
    species_traces_full = np.zeros((how_many_simulations,N_species))
    resource_traces_full = np.zeros((how_many_simulations,N_resources))

    sims = np.arange(how_many_simulations)

    
    sigmas = np.sqrt(np.logspace(-8,-1,100,base=10))
    for ii,sigma in enumerate(sigmas):
        print('Doing simulations for sigma ',ii,sigma)
        

        simulation_arguments = (how_many_simulations,N_species,N_resources,independent_species_cutoff,nsteps,
                                 growth_choices,entry_rates,inv_yields,
                                 dilution,dt,sigma)
        
        args = [(sim,simulation_arguments) for sim in sims]    

        print('Initializing pool..')
        pool = mp.Pool(10)

        print('Going into multiprocecssing...')
        time0 = time.time()
        with pool as p:
            results = p.map(main_simulation,args)

            
        print('Multiprocessing finished. Time taken: {:.3f}s'.format(time.time()-time0))
        time.sleep(5)

        for result in results:
             sim,resource_abus,species_abus = result

             species_traces_full[sim,:] = species_abus
             resource_traces_full[sim,:] = resource_abus
        
        species_traces_full = np.array(species_traces_full)
        resource_traces_full = np.array(resource_traces_full)

        print('Final arrays defined. Saving them to output...')

        with open('species_traces_full.{}.dat'.format(ii),'wb') as f:
            pickle.dump(species_traces_full,f)
            f.close()
        with open('resource_traces_full.{}.dat'.format(ii),'wb') as f:
            pickle.dump(resource_traces_full,f)
            f.close()
    
