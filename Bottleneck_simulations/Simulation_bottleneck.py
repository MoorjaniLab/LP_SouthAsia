#!/usr/bin/env python
# coding: utf-8

# **Bottleneck simulations to mimick Steppe admixture in South Asians (to test for LCT frequency in South Asians)**

import tskit
import msprime
import numpy as np
import math
import matplotlib.pyplot as plt
import pandas as pd
import demes
import demesdraw
import time
import random

#####
## **Functions**
#####

##My function, input: msprime simulation with samples from introgressed population + 1 Outgroup (here 'Steppe')
##ouput: for each migrating segment, the length and the number of Individuals who carried it (number of leaves below the migration node)
def get_migrating_tracts_count(ts):
    outgroup_id = [p.id for p in ts.populations() if p.metadata['name']=='Steppe'][0] ##find the id of the outgroup
    migrating_tracts = []
    n_leaves=[]
    # Get all tracts that migrated into the outgroup population (backwards)
    for migration in ts.migrations():
        if migration.dest == outgroup_id:
            mig_node=migration.node
            #find the number of leaves below the node
            for tree in ts.trees():
                if tree.interval.left>migration.right:
                    break
                elif tree.interval.left>=migration.left and tree.interval.right<=migration.right:
                    migrating_tracts.append((tree.interval.left, tree.interval.right))
                    n_leaves.append(len(list(tree.leaves(mig_node))))
    return np.array(migrating_tracts),np.array(n_leaves)


#Get ancestry tract for one individual when multiples are samples
def get_migrating_tracts_oneind(ts,ind1):
    outgroup_id = [p.id for p in ts.populations() if p.metadata['name']=='Steppe'][0] ##find the id of the outgroup
    migrating_tracts = []
    n_leaves=[]
    # Get all tracts that migrated into the neanderthal population (backwards)
    for migration in ts.migrations():
        if migration.dest == outgroup_id:
            mig_node=migration.node
            ##find if ind is leave of node
            for tree in ts.trees():
                if tree.interval.left>=migration.left and tree.interval.right<=migration.right and ind1 in list(tree.leaves(mig_node)):
                    migrating_tracts.append((tree.interval.left, tree.interval.right))
    return np.array(migrating_tracts)


###Count % genome above thresold, from output of get_migrating_tracts_count():
def percentage_above_limit(freqs,length_tracts,lim):
    den=0
    num=0
    for i in range(len(freqs)):
        if (np.isnan(length_tracts[i])==False):
            den+=length_tracts[i]
            if (freqs[i]>lim):
                num+=length_tracts[i]
            
    return num/den


##compute value of percentile from count_per_bin (in frequence) and bin_width:
def percentile_value(count_per_bin, bin_width, percentile):
    cumm=0
    i=0
    while(cumm+count_per_bin[i]<percentile):
        cumm+=count_per_bin[i]
        i+=1     
    
    return max((i-1)*bin_width+((percentile-cumm)/count_per_bin[i])*bin_width,0)

##Compare migrating tracts, and merge them (add counts if similar)
def merge_similar_migrating_tracts(array_migrating_tracts,array_count):
    df = pd.DataFrame(data=migrating_tracts, columns=['start','end'])
    df['counts']=counts
    aggregation_functions = {'counts': 'sum', 'start': 'first', 'end': 'first'}
    df_new = df.groupby(df['start']).aggregate(aggregation_functions)
    df_new['lengths']= df_new['end']-df_new['start']
    return df_new



# -----------------------------------------------------------------------------------------------------
# Parameters for demography
# -----------------------------------------------------------------------------------------------------

# Times are provided in years, so we convert into generations.
gen_time = 28.0 #generation time
rec_rate = 1.2e-8 #recombination rate
mutation_rate = 1.2e-8 #mutation rate

N_SAsia = 10_000
N_Steppe=10_000
N_ancWE=10_000
T_europa_asia = 30_000 / gen_time

# Admixture times and proportions from archaic humans
admixtureproportion = 0.078 #For Toda-like for Gujjar: 0.21
T_Steppe_into_SAS = 3_500 / gen_time

demography_simple = msprime.Demography() #initializing an msprime "demography" object
# This is the "trunk" population that we merge other populations into
demography_simple.add_population(
    name="SAHG",
    description="Ancestral South Asians",
    initial_size=N_SAHG,
    initially_active=True,
);

demography_simple.add_population(
    name="SouthAsia",
    description="SouthAsian",
    initial_size=N_SAsia,
    initially_active=True,
)    
demography_simple.add_population(
    name="Steppe",
    description="Steppe",
    initial_size=N_Steppe,
    initially_active=True,
)

# Split Steppe and South Asians
demography_simple.add_population_split(
    time=T_europa_asia, derived=["SouthAsia"], ancestral="SAHG"
)

demography_simple.add_population_split(
    time=T_europa_asia, derived=["Steppe"], ancestral="SAHG"
)

# Add Steppe geneflow for one generation (pulse model)
demography_simple.add_migration_rate_change(time = T_Steppe_into_SAS, rate = admixtureproportion, dest = "Steppe", source = "SouthAsia")
demography_simple.add_migration_rate_change(time = T_Steppe_into_SAS+1, rate = 0, dest = "Steppe", source = "SouthAsia")


demography_simple.sort_events()



#plot the demography
graph_simple = demography_simple.to_demes()
fig, ax = plt.subplots(figsize=(10,6))  # use plt.rcParams["figure.figsize"]
demesdraw.tubes(graph_simple, ax=ax, seed=1)
plt.show()


####
## Tests scenario and percentage of admixture
####

##Simulation to test scenario and functions
n_sample_A=5
samples_A = msprime.SampleSet(n_sample_A,population="SouthAsia")
samples_S = msprime.SampleSet(1,population="Steppe")
sample_list_simple=[samples_A,samples_S]

##Recombination map 
#to model chromosome 14
full_map_chr14=msprime.RateMap.read_hapmap("formsprime_genetic_map_chr14.txt" )


# Simulate chr14
length_chr14=106880170 ##as in chr14
sequence_length=length_chr14
elapsed = time.time()
ts = msprime.sim_ancestry(
    ploidy=1, ## test to change
    samples=sample_list_simple, scrip_for_python.sh
    demography=demography_simple,
    record_migrations=True,
    record_provenance=False,
    sequence_length=sequence_length, 
    #recombination_rate=rec_rate,
    recombination_rate=full_map_chr14,
    random_seed=None)

print(time.time() - elapsed)


##compute percentage of outgroup ancestry per indiv in different simulations (test)
elapsed=time.time()
pct_n=[]
for i in range(10):
    ts = msprime.sim_ancestry(
        ploidy=1,
        samples=sample_list_simple, 
        demography=demography_simple,
        record_migrations=True, 
        sequence_length=sequence_length,
        recombination_rate=full_map_chr14,
        record_provenance=False,
        random_seed=None)
    mig_tracts_ind0=get_migrating_tracts_oneind(ts,0)
    if (len(mig_tracts_ind0)):
        sum_mig0=sum(mig_tracts_ind0[:,1]-mig_tracts_ind0[:,0])
    else:
        sum_mig0=0
    
    pct_n.append(sum_mig0/sequence_length)
    
print("Nb ite:",len(pct_n),"mean outgroup ancestry %",np.mean(pct_n))
print(time.time() - elapsed)

#Distribution of % of outgroup ancestry per individual
plt.figure(figsize=(10,6))
plt.hist(pct_n,facecolor='g')
plt.xlabel('Fraction of Bases with Outgroup ancestry');plt.ylabel('Number of Simulated Chromosomes')
plt.show()

##Multiple samples in a simulation

n_sample_A=500
samples_A = msprime.SampleSet(n_sample_A,population="SouthAsia",ploidy=1)
samples_S = msprime.SampleSet(1,population="Steppe",ploidy=1)
sample_list_simple=[samples_A,samples_S]

# Simulate
elapsed = time.time()
#sequence_length=100_000_000
ts = msprime.sim_ancestry(
    ploidy=1, ## test to change
    samples=sample_list_simple, 
    demography=demography_simple,
    record_migrations=True,
    sequence_length=sequence_length, 
    recombination_rate=full_map_chr14,
    #recombination_rate=rec_rate,
    record_provenance=False,
    random_seed=None)

print(time.time() - elapsed)


# In[430]:


##Frequency of outgroup tracts
##Get segments and nb of individuals carring them
elapsed = time.time()
migrating_tracts, counts = get_migrating_tracts_count(ts)
print(time.time() - elapsed)

## Some migrating tracts overlapping, for now only add counts to each other
df_migrating_tracts=merge_similar_migrating_tracts(migrating_tracts,counts)

freqs=df_migrating_tracts['counts']/n_sample_A ##convert nb in frequency 

##Length of segments
length_tracts=df_migrating_tracts['lengths']
sum_mig=sum(length_tracts) ##total length

##For genome-wide distribution need to add segment without ancestry
no_outgroup=sequence_length-sum_mig 
length_tracts=np.append(length_tracts,no_outgroup)
freqs=np.append(freqs,0)

print(len(length_tracts),len(freqs))



###Generate table for histogram of the percentage of genome having %of outgroup frequency 
bin_width=0.01
count_per_bin=np.zeros(int(1/bin_width)+1)
for i in range(len(freqs)):
    nb_bin=math.floor(freqs[i]/bin_width)
    count_per_bin[nb_bin]+=length_tracts[i]/sequence_length
    
count_per_bin=count_per_bin/sum(count_per_bin)

bin_scale=np.arange(0,1.0001,bin_width)

##plot
plt.figure(figsize=(10,6))
plt.plot(bin_scale,count_per_bin)
plt.xlabel('Genome-Wide Outgroup Frequency');plt.ylabel('Distribution')
plt.show()


###
## Bootstrap constant scenario
###
n_simu=50
n_boot=50
bin_width=0.01

lim_togive=0.65 ##frequency in the pop we simulate
admixtureproportion=0.078 ## admixture porportion in the pop we want to simulate


##simulate and keep info
for i_simu in range(n_simu):
    sum_mig=0
    if (i_simu%5==0):
        print(i_simu)
    ts = msprime.sim_ancestry(
        ploidy=1,
        samples=sample_list_simple, 
        demography=demography_simple,
        record_migrations=True, 
        sequence_length=sequence_length,
        #recombination_rate=rec_rate,
        recombination_rate=full_map_chr14,
        record_provenance=False,
        random_seed=None)

    migrating_tracts, counts = get_migrating_tracts_count(ts)
    if (len(migrating_tracts)>0):
        ## Some migrating tracts overlapping, for now only add counts to each other
        df_migrating_tracts=merge_similar_migrating_tracts(migrating_tracts,counts)
        freqs=df_migrating_tracts['counts']/n_sample_A ##convert nb in frequency 
        ##Length of segments
        length_tracts=df_migrating_tracts['lengths']
        sum_mig=sum(length_tracts) ##total length
        
    ##for global on the genome (remove last variable to only consider archaic segments)
    no_archaic=sequence_length-sum_mig
    length_tracts=np.append(length_tracts,no_archaic)
    freqs=np.append(freqs,0)
    if (i_simu==0):
        df_for_boot=pd.DataFrame(length_tracts,columns =["length_tracts"+str(i_simu)])
        df_for_boot.insert(1,"freqs"+str(i_simu),freqs)
    else:
        df_for_boot_add=pd.DataFrame(length_tracts,columns =["length_tracts"+str(i_simu)])
        df_for_boot_add.insert(1,"freqs"+str(i_simu),freqs)
        df_for_boot = pd.concat([df_for_boot,df_for_boot_add],axis=1)
            

###Bootstrap
#sample i_simu simulation, n_boot times
##collect gw_95 and above20 param
gw_95_vec=np.zeros(n_boot)
gw_above20_vec=np.zeros(n_boot)
gw_above_limit_vec=np.zeros(n_boot)

for ite_n_boot in range(n_boot):
    count_per_bin=np.zeros(int(1/bin_width)+1)

    gw_above_limit=0
    gw_above20_limit=0

    for ite_i_simu in range(i_simu):
        bt_i=random.randint(0,n_simu-1)
        bin_for_distri=[math.floor(float(x)) for x in df_for_boot["freqs"+str(bt_i)]/bin_width if str(x) != 'nan']
        for i_bin in range(len(bin_for_distri)):
            count_per_bin[bin_for_distri[i_bin]]+=df_for_boot["length_tracts"+str(bt_i)][i_bin]

        gw_above_limit+=percentage_above_limit(df_for_boot["freqs"+str(bt_i)],df_for_boot["length_tracts"+str(bt_i)],lim_togive)/i_simu
        gw_above20_limit+=percentage_above_limit(df_for_boot["freqs"+str(bt_i)],df_for_boot["length_tracts"+str(bt_i)],0.2)/i_simu

    count_per_bin=count_per_bin/sum(count_per_bin)

    gw_above_limit_vec[ite_n_boot]=gw_above_limit
    gw_95_vec[ite_n_boot]=percentile_value(count_per_bin, bin_width, 0.95)
    gw_above20_vec[ite_n_boot]=gw_above20_limit

all_results = pd.DataFrame([["Constant",N_SAsia,0,np.mean(gw_95_vec),np.var(gw_95_vec),np.mean(gw_above20_vec),np.var(gw_above20_vec),np.mean(gw_above_limit_vec),np.var(gw_above_limit_vec)]],columns=['Scenario','Ne', 't', 'mean_gw95','var_gw_95','mean_gw_above20','var_gw_above20','mean_gw_abovelim','var_gw_abovelim'])
gw_95= pd.DataFrame(gw_95_vec, columns=["constant"])
gw_above20= pd.DataFrame(gw_above20_vec, columns=["constant"])
gw_above_lim= pd.DataFrame(gw_above_limit_vec, columns=["constant"])


####
## Bottleneck simple
####
size_bt_vec=[5000,1000,100,90,70,50,30,10]

##Bottleneck
start_bt=20
end_bt_vec=[30,40,50,70,120] ## can change to fit pop too


##
# Bootstrap/simu parameters
##
n_simu=50
n_boot=50
bin_width=0.01

for size_bt in size_bt_vec:
    print("size_bt=",size_bt)
    for start_bt in end_bt_vec:
        print("start_bt=",start_bt, "end_bt=", end_bt)
        
        ##scenario
        demography_bt = msprime.Demography() #initializing an msprime "demography" object
        # This is the "trunk" population that we merge other populations into
        demography_bt.add_population(
            name="SAHG",
            description="Ancestral South Asians",
            initial_size=N_SAHG,
            initially_active=True,
        );
        
        demography_bt.add_population(
            name="SouthAsia",
            description="SouthAsian",
            initial_size=N_SAsia,
            initially_active=True,
        )    
        demography_bt.add_population(
            name="Steppe",
            description="Steppe",
            initial_size=N_Steppe,
            initially_active=True,
        )
        
        # Split Steppe and South Asians
        demography_bt.add_population_split(
            time=T_europa_asia, derived=["SouthAsia"], ancestral="SAHG"
        )
        demography_bt.add_population_split(
            time=T_europa_asia, derived=["Steppe"], ancestral="SAHG"
        )
        # Add Steppe geneflow for one generation (pulse model)
        demography_bt.add_migration_rate_change(time = T_Steppe_into_SAS, rate = admixtureproportion, dest = "Steppe", source = "SouthAsia")
        demography_bt.add_migration_rate_change(time = T_Steppe_into_SAS+1, rate = 0, dest = "Steppe", source = "SouthAsia")
        
        # Add Bottleneck in Asian population
        demography_bt.add_population_parameters_change(
            time=start_bt, population="SouthAsia", initial_size=N_SAsia
        )
        
        demography_bt.add_population_parameters_change(
            time=end_bt, population="SouthAsia", initial_size=size_bt
        )
        
        demography_bt.sort_events()
        
        ##simulations
        ##simulate and keep info
        for i_simu in range(n_simu):
            sum_mig=0
            if (i_simu%5==0):
                print(i_simu)
            ts = msprime.sim_ancestry(
                ploidy=1,
                samples=sample_list_simple, 
                demography=demography_bt,
                record_migrations=True, 
                sequence_length=sequence_length,
                recombination_rate=full_map_chr14,
                record_provenance=False,
                random_seed=None)
        
            migrating_tracts, counts = get_migrating_tracts_count(ts)
            if (len(migrating_tracts)>0):
                ## Some migrating tracts overlapping, for now only add counts to each other
                df_migrating_tracts=merge_similar_migrating_tracts(migrating_tracts,counts)
                freqs=df_migrating_tracts['counts']/n_sample_A ##convert nb in frequency 
                ##Length of segments
                length_tracts=df_migrating_tracts['lengths']
                sum_mig=sum(length_tracts) ##total length
                
            ##for global on the genome 
            no_outgroup=sequence_length-sum_mig
            length_tracts=np.append(length_tracts,no_outgroup)
            freqs=np.append(freqs,0)
            if (i_simu==0):
                df_for_boot=pd.DataFrame(length_tracts,columns =["length_tracts"+str(i_simu)])
                df_for_boot.insert(1,"freqs"+str(i_simu),freqs)
            else:
                df_for_boot_add=pd.DataFrame(length_tracts,columns =["length_tracts"+str(i_simu)])
                df_for_boot_add.insert(1,"freqs"+str(i_simu),freqs)
                df_for_boot = pd.concat([df_for_boot,df_for_boot_add],axis=1)
                    
        ##Bootstrap
        ###Bootstrap
        #sample i_simu simulation, n_boot times
        ##collect gw_95 and above20 param
        gw_95_vec=np.zeros(n_boot)
        gw_above20_vec=np.zeros(n_boot)
        gw_above_limit_vec=np.zeros(n_boot)
        
        for ite_n_boot in range(n_boot):
            print("n_boot=",ite_n_boot)
            count_per_bin=np.zeros(int(1/bin_width)+1)
        
            gw_above_limit=0
            gw_above20_limit=0
            gw_above_limit=0
        
            for ite_i_simu in range(i_simu):
                bt_i=random.randint(0,n_simu-1)
                bin_for_distri=[math.floor(float(x)) for x in df_for_boot["freqs"+str(bt_i)]/bin_width if str(x) != 'nan']
                for i_bin in range(len(bin_for_distri)):
                    count_per_bin[bin_for_distri[i_bin]]+=df_for_boot["length_tracts"+str(bt_i)][i_bin]
        
                gw_above_limit+=percentage_above_limit(df_for_boot["freqs"+str(bt_i)],df_for_boot["length_tracts"+str(bt_i)],lim_togive)/i_simu
                gw_above20_limit+=percentage_above_limit(df_for_boot["freqs"+str(bt_i)],df_for_boot["length_tracts"+str(bt_i)],0.2)/i_simu
        
            count_per_bin=count_per_bin/sum(count_per_bin)

            gw_above_limit_vec[ite_n_boot]=gw_above_limit
            gw_95_vec[ite_n_boot]=percentile_value(count_per_bin, bin_width, 0.95)
            gw_above20_vec[ite_n_boot]=gw_above20_limit
            
        i=len(all_results)
        all_results.loc[i+1]=[str(size_bt)+"bt"+str(end_bt-start_bt),size_bt,end_bt-start_bt,np.mean(gw_95_vec),np.var(gw_95_vec),np.mean(gw_above20_vec),np.var(gw_above20_vec),np.mean(gw_above_limit_vec),np.var(gw_above_limit_vec)]  
        gw_95.insert(1,str(size_bt)+"bt"+str(end_bt-start_bt),gw_95_vec)
        gw_above20.insert(1,str(size_bt)+"bt"+str(end_bt-start_bt),gw_above20_vec)
        gw_above_lim.insert(1,str(size_bt)+"bt"+str(end_bt-start_bt),gw_above_limit_vec)


####
##write FINAL
####
sum_results="results_metrics_per_scenario_start20bt"
outfile_gw95="gw_95_start20bt"
outfile_gw20="gw_above20_start20bt"
outfile_gw_limit="gw_above_limit_start20bt"

all_results.to_csv(sum_results+str(n_simu)+"_"+str(n_boot)+".txt", index=None, sep=' ', mode='w')
gw_95.to_csv(outfile_gw95+str(n_simu)+"_"+str(n_boot)+".txt", index=None, sep=' ', mode='w')
gw_above20.to_csv(outfile_gw20+str(n_simu)+"_"+str(n_boot)+".txt", index=None, sep=' ', mode='w')
gw_above_lim.to_csv(outfile_gw_limit+str(n_simu)+"_"+str(n_boot)+".txt", index=None, sep=' ', mode='w')




