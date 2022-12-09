import sys
import os
import msprime
import numpy as np
import scipy
from IPython.display import SVG
import matplotlib.pyplot as plt
import matplotlib.pyplot as ml
import demesdraw


print('Simulated parameters:')
print('number of windows:', sys.argv[1])
print('T1:', sys.argv[2])
print('Ne after T1:', sys.argv[3])
print('T2:', sys.argv[4])
print('Ne after T2:', sys.argv[5])
print('T2:', sys.argv[6])
print('Ne after T2:', sys.argv[7])
print('output file:', sys.argv[8])


def demogr_model_0_1(t1,t2,t3,NeT1,NeT2,NeT3):
    demography = msprime.Demography()

    Ne = 1e4

    #initializing populations
    demography.add_population(name="gra", initial_size=Ne*1.2)
    demography.add_population(name="alb", initial_size=Ne*40)
    demography.add_population(name="sam", initial_size=Ne)
    demography.add_population(name="agu", initial_size=Ne)
    demography.add_population(name="T1", initial_size=NeT1)
    demography.add_population(name="T2", initial_size=NeT2)
    demography.add_population(name="T3", initial_size=NeT3)

    #adding split times
    demography.add_population_split(time=t1, derived=["gra","sam"], ancestral="T1")
    demography.add_population_split(time=t2, derived=["T1","alb"], ancestral="T2")
    demography.add_population_split(time=t3, derived=["T2","agu"], ancestral="T3")

    #setting up gene flow
    #demography.set_migration_rate("alb", "agu", mig_rate)

    #setting up admixture
    #demography.add_admixture(time=7, derived="alb", ancestral=["alb_anc", "agu"], proportions=[0.25, 0.75])

    #Simultanious test
    ml.rcParams['figure.figsize'] = (8.0, 5.0)
    graph = msprime.Demography.to_demes(demography)
    fig, ax = plt.subplots()  # use plt.rcParams["figure.figsize"]
    demesdraw.tubes(graph, ax=ax, seed=1)
    plt.savefig('demography'+sys.argv[2]+'_'+sys.argv[3]+'_'+sys.argv[4]+'_'+sys.argv[5]+'_'+sys.argv[6]+'_'+sys.argv[7]+'.pdf')


    return demography


def simulate_windows(num_replicates,demography):
    ancestry_reps = msprime.sim_ancestry(samples={"gra": 12, "alb": 14, "sam": 4, "agu": 12},
                                         demography=demography, sequence_length=10000,
                                         num_replicates=num_replicates, ploidy=1)
    for ts in ancestry_reps:
        mutated_ts = msprime.sim_mutations(ts, rate=0.9e-8)
        yield mutated_ts


def ts_newick(path,genealogies):
    #file = open(path, "w")

    #if model == 'one_way':
    #    genealogies=tern_simulate_P3P2_s()

    #if model == 'two_way':
    #    genealogies=tern_simulate_twoway_s()

    for replicate_index, ts in enumerate(genealogies):
        for t in ts.trees():
            newick=t.newick(precision=1)
            replace_strings_tens=ts_newick_rename_dict_tens()
            replace_strings_hun=ts_newick_rename_dict_hun()
            for word in replace_strings_tens.items():
                newick = newick.replace(str(word[0]), str(word[1]))
                for word in replace_strings_hun.items():
                    newick = newick.replace(str(word[0]), str(word[1]))
            #print(newick+"\n",flush=True)
            #file.write(newick+"\n")

    #file.close()

def ts_newick_rename_dict_hun():
    sample_names = dict()

    #{"gra": 12, "alb": 14, "sam": 4, "agu": 12}

    nP1=12
    nP2=14
    nP3=4
    n0=12

    #Total number of empirical samples
    ntotal=nP1+nP2+nP3+n0

    #Creating list of newick names
    ts_name=list()
    for i in range(1,ntotal+1):
        ts_name.append(","+str(i)+":")

    #Creating list of population samples
    nP1_names=list()
    for i in range(1,nP1+1):
        nP1_names.append(",gra_"+str(i)+":")

    nP2_names=list()
    for i in range(1,nP2+1):
        nP2_names.append(",alb_"+str(i)+":")

    nP3_names=list()
    for i in range(1,nP3+1):
        nP3_names.append(",sam_"+str(i)+":")

    n0_names=list()
    for i in range(1,n0+1):
        n0_names.append(",agu_"+str(i)+":")

    pop_names=nP1_names+nP2_names+nP3_names+n0_names

    #Using dictionary comprehension to convert lists to dictionary
    sample_names = {ts_name[i]: pop_names[i] for i in range(len(ts_name))}

    return sample_names

def ts_newick_rename_dict_tens():
    sample_names = dict()

    nP1=12
    nP2=14
    nP3=4
    n0=12

    #Total number of empirical samples
    ntotal=nP1+nP2+nP3+n0

    #Creating list of newick names
    ts_name=list()
    for i in range(1,ntotal+1):
        ts_name.append("("+str(i)+":")

    #Creating list of population samples
    nP1_names=list()
    for i in range(1,nP1+1):
        nP1_names.append("("+"gra_"+str(i)+":")

    nP2_names=list()
    for i in range(1,nP2+1):
        nP2_names.append("("+"alb_"+str(i)+":")

    nP3_names=list()
    for i in range(1,nP3+1):
        nP3_names.append("("+"sam_"+str(i)+":")

    n0_names=list()
    for i in range(1,n0+1):
        n0_names.append("("+"agu_"+str(i)+":")

    pop_names=nP1_names+nP2_names+nP3_names+n0_names

    #Using dictionary comprehension to convert lists to dictionary
    sample_names = {ts_name[i]: pop_names[i] for i in range(len(ts_name))}

    return sample_names


#ts_newick(sys.argv[8],int(sys.argv[1]),demogr_model_0_1(float(sys.argv[2]),float(sys.argv[4]),float(sys.argv[6]),float(sys.argv[3]),float(sys.argv[5]),float(sys.argv[7])))
print("Creating model")
demography=demogr_model_0_1(float(sys.argv[2]),float(sys.argv[4]),float(sys.argv[6]),float(sys.argv[3]),float(sys.argv[5]),float(sys.argv[7]))

print("Simulating windows")
genealogies=simulate_windows(int(sys.argv[1]),demography)

print("Parsing windows")
#ts_newick(sys.argv[8],int(sys.argv[1]),demogr_model_0_1(float(sys.argv[2]),float(sys.argv[4]),float(sys.argv[6]),float(sys.argv[3]),float(sys.argv[5]),float(sys.argv[7])))

file = open(sys.argv[8], "w")
for replicate_index, ts in enumerate(genealogies):
    for t in ts.trees():
        newick=t.newick(precision=1)
        replace_strings_tens= ts_newick_rename_dict_tens()
        replace_strings_hun=ts_newick_rename_dict_hun()
        for word in replace_strings_tens.items():
            newick = newick.replace(str(word[0]), str(word[1]))
            for word in replace_strings_hun.items():
                newick = newick.replace(str(word[0]), str(word[1]))
        #print(newick+"\n",flush=True)
        file.write(newick+"\n")

    #file.close()

#print("Simulated windows")
#ts_newick(sys.argv[8],genealogies)
print("Converted to newick")

#conda activate dasha
#cd GitHub/BW_wagtail_conundrum/
