import os
from class_individu import individu
import configparser
import pandas as pd
import numpy as np
import collections as col
from pylab import *
import errno
import csv
import matplotlib.pyplot as plt

RNAPs_genSC = 0.1
###########################################################
#                       Functions                         #
###########################################################

# Read the config files
def read_config_file(path):
    config = configparser.ConfigParser()
    # to preserve capital letters
    config.optionxform = str
    config.read(path)
    return config

###################### Reading files ######################

# you can combine those two functions
def load_gff(filename):
    gff_df_raw = pd.read_table(filename, sep='\t', comment='#', header=0)
    return gff_df_raw

def load_tab_file(filename):
    data = pd.read_table(filename, sep='\t', header=0)
    return data

def str2num(s):
    s[s == '+'] = 1 #True
    s[s == '-'] = -1 #False
    return s

def get_tr_nbr_csv(csv_file):
    csv_tr_nbr = pd.read_csv(csv_file, sep=';', header=None)
    tr_nbr = csv_tr_nbr.values
    return tr_nbr.flatten()

######################## Others ###########################

# Get the genome size from the header of gff file (befor renaming it)
def get_genome_size(gff_df):
    genome_size = int(gff_df.columns[4]) - int(gff_df.columns[3])
    return genome_size

# Rename the header (columns names)
def rename_gff_cols(gff_df):
    names=["seqid", "source", "type","start","end","score","strand","phase","attributes"]
    gff_df.columns = names
    return gff_df

# Whether the gene is on the + strand or - strand
def in_forward(tr_id):
    if strands[tr_id] == 1. :
        return True
    else:
        return False

# Get the transciption unit with the list of transcripts
def get_TU(TUindex_list):
    TU_dict = col.defaultdict(list)
    for index, TUindex in enumerate(TUindex_list):
        TU_dict[TUindex].append(index)
    return TU_dict

# calculate the initiation rate
def f_init_rate(tr_prob, sig, sigma_t, epsilon, m):
    tr_prob_sig = tr_prob * np.exp((1/(1+np.exp((sig-sigma_t)/epsilon)))*m)
    return tr_prob_sig

# Get the list of all the possible transcripts
def get_tr_info(tss, tts, TU_tts, Kon, Poff):
    this_TU_tts = []
    tr_id = []
    tr_start = []
    tr_end = []
    tr_strand = []
    tr_size = []
    tr_rate = []
    sum_Kon = np.sum(Kon)

    j = 0 # trancript id indice
    for i in tss.index.values: # All TSSs
        # get the TU of this tss
        TU_id = tss['TUindex'][i]
        # the list of TTS that are in the same TU of this tss_id (i)
        # TU_tts ex : defaultdict(list, {0: [1150, 2350], 1: [6250]})
        # On prend tt les tts qui existent dans la meme UT du tss choisi
        this_TU_tts = TU_tts[TU_id] # pour 0 => [1150, 2350]
        # + or -
        if tss['TUorient'][i] == '+' :
            # go right
            # tr_rate ======>  [ 0.1875   0.03125  0.00625  0.025    0.0625   0.25     0.125    0.3125 ]
            k = TU_id # TTS id index : k start from the first position of each TU
            proba_rest = 1
            while proba_rest > 0 :
                if tss['TSS_pos'][i] < tts['TTS_pos'][k]:
                    tr_id.append(j)
                    tr_strand.append(1)
                    tr_start.append(tss['TSS_pos'][i])
                    # after getting the TSSs, we shall (in every loop) generate a new tr_end
                    tr_end.append(tts['TTS_pos'][k])
                    # the probability to choose a specific transcript
                    tr_rate.append(Kon[i] * (Poff[k] * proba_rest))
                    proba_rest = (1 - Poff[k]) * proba_rest
                    j += 1
                k += 1
        else:
            # go leftB
            k = TU_id #0 #len(this_TU_tts)# TTS id index ### [0 1 2 3 4 5]
            proba_rest = 1
            while proba_rest > 0 : # and k < len(this_TU_tts) : # >= 0 : #
                if tts['TTS_pos'][k] < tss['TSS_pos'][i] : # and tts['TUindex'][k] == TU_id :
                    tr_id.append(j)
                    tr_strand.append(-1)
                    tr_start.append(tss['TSS_pos'][i])
                    # after getting them, we shall (in every loop) generate a new tr_end
                    # tr_end.append(this_TU_tts[k])
                    tr_end.append(tts['TTS_pos'][k])
                    # the probability to choose a specific transcript
                    tr_rate.append(Kon[i] * (Poff[k] * proba_rest))
                    proba_rest = (1 - Poff[k]) * proba_rest
                    j += 1
                k += 1
    tr_size = np.abs(np.array(tr_start) - np.array(tr_end))
    ts_beg_all_trs = np.zeros(len(tr_id), dtype=int64)
    ts_remain_all = np.around(tr_size)
    return (tr_id, tr_strand, tr_start, tr_end, tr_rate, tr_size, ts_beg_all_trs, ts_remain_all)

def f_prob_init_rate(init_rate, sum_init_rate, DELTA_T):
    return (1-np.exp(-sum_init_rate*DELTA_T)) * (init_rate/sum_init_rate)

def f_prob_unhooked_rate(sum_Kon, DELTA_T, RNAPs_unhooked_nbr):
    return np.exp(-sum_Kon*DELTA_T)/RNAPs_unhooked_nbr

# Get the transciption unit with the list of tts belonging to TU.
def get_TU_tts(tss, tts):
    TU_tts = col.defaultdict(list)
    for index, TUindex in enumerate(tss['TUindex'].values):
        TU_tts[TUindex].append(tts['TTS_pos'][index])
    return TU_tts

def calc_sigma(Barr_sigma, GYRASE_CONC, k_GYRASE, x0_GYRASE, GYRASE_CTE, TOPO_CONC, k_TOPO, x0_TOPO, TOPO_CTE, DELTA_T): #RNAPs_genSC

    d_sigma = (-GYRASE_CONC*1/(1+np.exp(-k_GYRASE*(Barr_sigma-x0_GYRASE)))*GYRASE_CTE + TOPO_CONC*1/(1+np.exp(k_TOPO*(Barr_sigma-x0_TOPO)))*TOPO_CTE) * DELTA_T
    Barr_sigma += d_sigma

    return Barr_sigma

def plotGenome(ind,title):
    fig = plt.figure()
    ax=plt.axes()

    pB= ax.scatter(ind.Barr_fix, np.ones(10), s=20, c='r',label="Barrieres")
    pTSS = ax.scatter(ind.TSS_pos, np.ones(10), s=20, c='g',label="TSS")
    pTTS = ax.scatter(ind.TTS_pos, np.ones(10), s=20, c='b',label="TTS")
    pline = ax.plot(np.array(range(ind.genome_size)), np.ones(ind.genome_size), alpha=0.3, c='black')
    for i in range(0,len(ind.TSS_pos)) :
        if ind.strands[i]>0 :
            ax.arrow(ind.TSS_pos[i]+2.5, 1, ind.TTS_pos[i]-ind.TSS_pos[i]-5.5, 0, head_width=0.01, head_length=1, fc='k', ec='k')
        else :
            ax.arrow(ind.TTS_pos[i]-2.5, 1, -(ind.TTS_pos[i]-ind.TSS_pos[i]-5.5), 0, head_width=0.01, head_length=1, fc='k', ec='k')
    plt.legend(handles=[pB,pTSS,pTTS])
    plt.ylim([0.5,1.5])
    for i, txt in enumerate(ind.noms_genes) :
        plt.annotate(txt, (min(ind.TSS_pos[i],ind.TTS_pos[i])+5,0.95))
    plt.title('Genome '+title)
    plt.axis('off')

    plt.show()

###########################################################
#              Main process (Simulation)                  #
###########################################################

def simulation():
    # Get params file
    INI_file=input("Nom du fichier de paramètres : ")
    config = read_config_file(INI_file)

    # Get inputs infos from the config file
    GFF_file = config.get('INPUTS', 'GFF')
    TSS_file = config.get('INPUTS', 'TSS')
    TTS_file = config.get('INPUTS', 'TTS')
    Prot_file = config.get('INPUTS', 'BARR_FIX')

    # Get values from the config file
    m = config.getfloat('GLOBAL', 'm')
    sigma_t = config.getfloat('GLOBAL', 'sigma_t')
    epsilon = config.getfloat('GLOBAL', 'epsilon')

    SIGMA_0 = config.getfloat('SIMULATION', 'SIGMA_0')
    DELTA_X = config.getfloat('SIMULATION', 'DELTA_X')
    DELTA_T = config.getfloat('SIMULATION', 'DELTA_T')
    RNAPS_NB = config.getint('SIMULATION', 'RNAPS_NB')
    ITERATIONS_NB = config.getfloat('SIMULATION', 'ITERATIONS_NB')
    OUTPUT_STEP = config.getfloat('SIMULATION', 'OUTPUT_STEP')

    GYRASE_CONC = config.getfloat('SIMULATION', 'GYRASE_CONC')
    TOPO_CONC = config.getfloat('SIMULATION', 'TOPO_CONC')
    TOPO_CTE = config.getfloat('SIMULATION', 'TOPO_CTE')
    GYRASE_CTE = config.getfloat('SIMULATION', 'GYRASE_CTE')
    TOPO_EFFICIENCY = config.getfloat('SIMULATION', 'TOPO_EFFICIENCY')
    k_GYRASE = config.getfloat('SIMULATION', 'k_GYRASE')
    x0_GYRASE = config.getfloat('SIMULATION', 'x0_GYRASE')
    k_TOPO = config.getfloat('SIMULATION', 'k_TOPO')
    x0_TOPO = config.getfloat('SIMULATION', 'x0_TOPO')

    p_inv = config.getfloat('SIMULATION', 'p_inv')
    p_indel = config.getfloat('SIMULATION', 'p_indel')
    p_keep  = config.getfloat('SIMULATION', 'p_keep')
    nb_iter = int(config.getfloat('SIMULATION', 'nb_iter'))

    gff_df_raw = load_gff(GFF_file)
    tss = load_tab_file(TSS_file)
    tts = load_tab_file(TTS_file)
    prot = load_tab_file(Prot_file)
    
    # Genome size
    genome_size = get_genome_size(gff_df_raw)
    gff_df = rename_gff_cols(gff_df_raw)
    params = {'sigma_t':sigma_t,'epsilon':epsilon, 'SIGMA_0':SIGMA_0, 'DELTA_X':DELTA_X, \
    'DELTA_T':DELTA_T, 'RNAPS_NB':RNAPS_NB, 'ITERATIONS_NB':ITERATIONS_NB, \
    'OUTPUT_STEP':OUTPUT_STEP, 'GYRASE_CONC':GYRASE_CONC, 'TOPO_CONC':TOPO_CONC, \
    'TOPO_CTE':TOPO_CTE, 'GYRASE_CTE':GYRASE_CTE, 'TOPO_EFFICIENCY':TOPO_EFFICIENCY, \
    'k_GYRASE':k_GYRASE, 'x0_GYRASE':x0_GYRASE, 'k_TOPO':k_TOPO, 'x0_TOPO':x0_TOPO, 'm':m}

    # Read Envir
    f_env = open("environment.dat", "r")
    env = f_env.readlines()
    genes = np.zeros(10, dtype=float)
    compteur = 0
    for line in env:
        genes[compteur] = line.split()[1]
        compteur += 1


    # Create individu
    ind = individu(gff_df, tss, tts, prot, genome_size, DELTA_X, genes, p_keep)
    
    # Initialize fitness
    genes_level = start_transcribing(params,ind)
    ind.fitness = ind.calcul_fitness(genes_level)

    # Plot genome initial
    plotGenome(ind,'Initial')

    # Simulations
    events = [0] # to save all events index
    genes_expr = [genes_level] # to save all genes levels
    fitnesses = [ind.fitness] # to save all fitnesses
    
    for i in range(1, nb_iter+1):
        
        print(i)
        
        inversion=False
        insertion=False
        deletion=False
        
        # Inversion
        p=np.random.rand()
        if p<p_inv:
            ind.inversion()
            inversion=True

        # Insertion Deletion
        p=np.random.rand()
        if p<p_indel:
            ev = ind.indel()
            if ev==1:
                insertion=True
            else:
                deletion=True
                
        # Get index of event
        if not inversion and not insertion and not deletion: events.append(0)
        elif inversion and not insertion and not deletion: events.append(1)
        elif inversion and insertion: events.append(2)
        elif inversion and deletion: events.append(3)
        elif not inversion and insertion: events.append(4)
        elif not inversion and deletion: events.append(5)

        # Set data frame attributs of ind (with values calculated in inversion and/or insertion/deletion
        ind.upgrade_new_pd_dataframe()

        # Simulate transcription and get genes level expression
        genes_lev = []
        for j in range(0,5) :
            genes_lev.append(start_transcribing(params,ind))
        genes_level = np.mean(genes_lev,axis=0)

        if i==nb_iter:
            genes_expr.append(genes_level)

        # Set new fitness
        ind.new_fitness=ind.calcul_fitness(genes_level)
        
        # Choose individual to keep according to fitnesses
        ind.choice_indiv()
        fitnesses.append(ind.fitness)


    # Plot final genome
    plotGenome(ind,'Final')

    # Plot Fitnesses over time with color for each event
    fig = plt.figure()
    ax = plt.axes()
    colormap = np.array(['grey', 'red', 'yellow','green', 'blue', 'm'])
    labels = np.array(['No event', 'Inversion', 'Inversion + insertion', 'Inversion + deletion', 'Insertion', 'Deletion'])
    X = np.array(range(len(fitnesses)))
    Y = np.array(fitnesses)
    ax.scatter(X, Y, s=20, c=colormap[events])
    ax.plot(X, Y, alpha=0.3, c='black')
    plt.title("Fitness en fonction du temps")
    plt.xlabel("Iterations")
    plt.ylabel("Fitness")
    plt.show()

    return(fitnesses,events,genes_expr)

###########################################################
#         Transcription Process (Simulation)              #
###########################################################
def start_transcribing(p, ind):
    # TSS_pos
    # TSS_pos = (tss['TSS_pos'].values/DELTA_X).astype(int)

    # Kon
    Kon = ind.newtss['TSS_strength'].values

    # Poff
    Poff = ind.newtts['TTS_proba_off'].values

    # Dict of transcription units with the list of tts belonging to TU.
    TU_tts = get_TU_tts(ind.newtss, ind.newtts)

    # The RNAPs id
    RNAPs_id = np.full(p['RNAPS_NB'], range(0, p['RNAPS_NB']), dtype=int)

    # The position of RNAPs
    RNAPs_pos = np.full(p['RNAPS_NB'], NaN) #np.zeros(RNAPS_NB, dtype=int)

    # RNAPs_last_pos
    RNAPs_last_pos = np.full(p['RNAPS_NB'], NaN) #np.zeros(RNAPS_NB, dtype=int)

    # Strands orientation
    strands=ind.newstrands

    # list of all possible transcripts
    tr_id, tr_strand, tr_start, tr_end, tr_rate, tr_size, ts_beg_all_trs, ts_remain_all = get_tr_info(ind.newtss, ind.newtts, TU_tts, Kon, Poff)

    # convert all variables to numpy array
    tr_id = np.array(tr_id)
    tr_strand = np.array(tr_strand)

    tr_start = np.array(tr_start)/p['DELTA_X']
    tr_start = tr_start.astype(int64)

    tr_end = np.array(tr_end)/p['DELTA_X']
    tr_end = tr_end.astype(int64)

    tr_rate = np.array(tr_rate)

    tr_size = np.array(tr_size)/p['DELTA_X']
    tr_size = tr_size.astype(int64)

    ts_beg_all_trs = np.array(ts_beg_all_trs)

    ts_remain_all = np.array(ts_remain_all)/p['DELTA_X']
    ts_remain_all = ts_remain_all.astype(int64)

    # The number of times transcripts has been transcribed
    tr_nbr = np.zeros(len(tr_id), dtype=int)

    #genome=ind.newgenome

    #Barr_fix = np.copy(ind.newBarr_fix)
    Barr_fix = (ind.newBarr_fix/p['DELTA_X']).astype(int)
    TSS_pos = (ind.newTSS_pos/p['DELTA_X']).astype(int)
    TTS_pos = (ind.newTTS_pos/p['DELTA_X']).astype(int)
    genome =  int(ind.newgenome/p['DELTA_X'])


    # just for the echo we can assign it directely
    Barr_pos = np.copy(Barr_fix)
    Dom_size = np.ediff1d(Barr_pos)
    Dom_size = np.append(Dom_size, genome-Barr_fix[-1]+Barr_fix[0]) # !! change Barr_fix to Barr_pos case : O | |

    Barr_type = np.full(len(Barr_fix), 0, dtype=int)
    Barr_sigma = np.full(len(Barr_fix), p['SIGMA_0'])

    # here we need to make an Barr_ts_remain
    # to track the position of each RNAPol
    # each position in Barr_ts_remain is associated with the same position in Barr_pos
    Barr_ts_remain = np.full(len(Barr_fix), NaN) # The Barr_ts_remain of fixed barr is NaN


    ######### Variables used to get the coverage ##########

    id_shift_fwd = list(range(1, genome))
    id_shift_fwd.append(0)
    id_shift_fwd = np.array(id_shift_fwd)
    id_shift_bwd = list(range(0, genome-1))
    id_shift_bwd.insert(0, genome-1)
    id_shift_bwd = np.array(id_shift_bwd)

    cov_bp = np.arange(0, ind.genome_size, p['DELTA_X'])
    cov_bp = np.resize(cov_bp, genome)

    #mmmm sigma = np.full(genome, SIGMA_0)

    ###########################################################
    #                 initiation of values                    #
    ###########################################################

    # save the time when RNApoly is starting trasncribing a specific transcript
    tr_times = col.defaultdict(list)

    # numpy array where will save all RNAPs info
    save_RNAPs_info = np.full([p['RNAPS_NB'], 2, int(p['ITERATIONS_NB']/p['DELTA_T'])], np.nan) # nbr d'ele (cols)

    # the same for transcripts info
    save_tr_info = np.full([len(tr_id), 2, int(p['ITERATIONS_NB']/p['DELTA_T'])], np.nan)

    # in those variables, we will save/append info in each time step to save them as --> all_res ;-)
    save_Barr_sigma = list()
    save_Dom_size = list()
    save_mean_sig_wholeGenome = list()

    ########### Go !

    RNAPs_unhooked_id = np.copy(RNAPs_id)

    RNAPs_strand = np.full(p['RNAPS_NB'], NaN)
    ts_beg = np.full(p['RNAPS_NB'], NaN)
    ts_remain = np.full(p['RNAPS_NB'], NaN)
    # RNAPs_tr will contain the id of the picked transcript
    RNAPs_tr = np.full(p['RNAPS_NB'], -1, dtype=(int64))
    # get the TSSs ids
    tss_id = ind.newtss.index.values

    # in the case of RNAP_NBR = 0
    RNAPs_hooked_id = []

    for t in range(0,int(p['ITERATIONS_NB']/p['DELTA_T'])):
        # we need to know each TSS belong to which Domaine
        TSS_pos_idx = np.searchsorted(Barr_pos, TSS_pos)

        # after knowing the domaine of each TSS we can get sigma
        sigma_tr_start = Barr_sigma[TSS_pos_idx-1]

        # get the initiation rates
        if (len(tr_rate)!= len(sigma_tr_start)) :
            print("Longueur tr_rate : "+str(len(tr_rate))+". Longueur sigma_tr_start : "+str(len(sigma_tr_start)))
        init_rate = f_init_rate(tr_rate, sigma_tr_start, p['sigma_t'], p['epsilon'], p['m'])

        sum_init_rate = np.sum(init_rate)

        prob_init_rate = f_prob_init_rate(init_rate, sum_init_rate, p['DELTA_T'])

        if np.size(RNAPs_unhooked_id)!=0:
            # get the unhooked rates
            prob_unhooked_rate = f_prob_unhooked_rate(sum_init_rate, p['DELTA_T'], len(RNAPs_unhooked_id))
            # craete the numpy array
            prob_unhooked_rate = np.full(len(RNAPs_unhooked_id), prob_unhooked_rate)
            all_prob = np.concatenate([prob_init_rate, prob_unhooked_rate])

            # create the numpy array that will contains [ nTSS , Unhooked RNAPS ]
            tss_and_unhooked_RNAPs = np.concatenate([tss_id, np.full(len(RNAPs_unhooked_id), -1, dtype=int)])

            picked_tr = np.random.choice(tss_and_unhooked_RNAPs, len(RNAPs_unhooked_id), replace=False, p=all_prob) #RNAPs_unhooked_id

            # This is the KEY !
            picked_tr_hooked_id = picked_tr[np.where(picked_tr!=-1)[0]]
            picked_tr_unhooked_id = picked_tr[np.where(picked_tr==-1)[0]]

            new_RNAPs_hooked_id = RNAPs_unhooked_id[np.where(picked_tr!=-1)[0]] ### Change

            RNAPs_tr[new_RNAPs_hooked_id] = picked_tr[picked_tr!=-1]
            RNAPs_strand[new_RNAPs_hooked_id] = tr_strand[picked_tr[np.where(picked_tr!=-1)]]

            # The new position of each polymerase
            # if there is no RNAP already at this position
            RNAPs_pos[new_RNAPs_hooked_id] = tr_start[picked_tr[np.where(picked_tr!=-1)]].astype(int)

            # take the position and use them to get the index in which u will insert them in Barr_pos array
            Barr_pos_RNAPs_idx = np.searchsorted(Barr_pos, RNAPs_pos[new_RNAPs_hooked_id])

            #after getting the idx, we start inserting
            Barr_pos = np.insert(Barr_pos, Barr_pos_RNAPs_idx, RNAPs_pos[new_RNAPs_hooked_id])
            Dom_size = np.ediff1d(Barr_pos)
            Dom_size = np.append(Dom_size, genome-Barr_pos[-1]+Barr_pos[0])
            Barr_type = np.insert(Barr_type, Barr_pos_RNAPs_idx, RNAPs_strand[new_RNAPs_hooked_id])

            # Now Sigma
            Barr_sigma = np.insert(Barr_sigma, Barr_pos_RNAPs_idx, Barr_sigma[Barr_pos_RNAPs_idx-1])

            # RNAPs_last_pos
            RNAPs_last_pos[new_RNAPs_hooked_id] = tr_end[picked_tr_hooked_id]
            ts_beg[new_RNAPs_hooked_id] = 0
            ts_remain[new_RNAPs_hooked_id] = ts_remain_all[picked_tr_hooked_id]
            Barr_ts_remain = np.insert(Barr_ts_remain, Barr_pos_RNAPs_idx, ts_remain[new_RNAPs_hooked_id])
            RNAPs_hooked_id = np.where(RNAPs_tr!=-1)[0]

        ts_beg[RNAPs_hooked_id] += 1
        ts_remain[RNAPs_hooked_id] -= 1

        # save the time when RNApoly FINISHS trasncribing a specific transcript
        for x in RNAPs_tr[np.where(ts_remain==0)] :
            tr_times[x].append(t*p['DELTA_T']) # + 0.5

        tr_nbr[RNAPs_tr[np.where(ts_remain==0)]]+=1

        Barr_ts_remain[np.where(Barr_type == -1)]-=1
        Barr_ts_remain[np.where(Barr_type == 1)]-=1

        # Get the index of RNAPs to remove
        rm_RNAPs_idx = np.where(Barr_ts_remain == 0)[0]

        # recover sigma value of the removed position
        removed_sigma = Barr_sigma[rm_RNAPs_idx]
        removed_dom_size = Dom_size[rm_RNAPs_idx]

        # recover the old_dom_size : the size of the previous domaine before combination/merging
        old_dom_size = Dom_size[rm_RNAPs_idx-1]
        old_sigma = Barr_sigma[rm_RNAPs_idx-1]


        # update Dom_size
        #Dom_size[rm_RNAPs_idx-1] += removed_dom_size
        # or
        Dom_size = np.ediff1d(Barr_pos)
        Dom_size = np.append(Dom_size, genome-Barr_fix[-1]+Barr_fix[0])

        Barr_sigma[rm_RNAPs_idx-1] = (old_dom_size*old_sigma+removed_dom_size*removed_sigma)/(old_dom_size+removed_dom_size)

        # and remove them
        Barr_pos = np.delete(Barr_pos, rm_RNAPs_idx)
        Barr_type = np.delete(Barr_type, rm_RNAPs_idx)
        Barr_ts_remain = np.delete(Barr_ts_remain, rm_RNAPs_idx)
        Barr_sigma = np.delete(Barr_sigma, rm_RNAPs_idx)
        Dom_size = np.delete(Dom_size, rm_RNAPs_idx)

        # update the RNAPs_tr array
        RNAPs_tr[np.where(ts_remain==0)] = -1
        # update the RNAPs_unhooked_id based on RNAPs_tr
        RNAPs_unhooked_id = np.where(RNAPs_tr==-1)[0]

        # reset the arrays
        RNAPs_strand[RNAPs_unhooked_id] = NaN
        RNAPs_pos[RNAPs_unhooked_id] = NaN
        RNAPs_last_pos[RNAPs_unhooked_id] = NaN
        ts_beg[RNAPs_unhooked_id] = NaN
        ts_remain[RNAPs_unhooked_id] = NaN

        Barr_pos[np.where(Barr_type == -1)]-=1
        Barr_pos[np.where(Barr_type == 1)]+=1

        # Update the position of polymerases still transcribing
        RNAPs_pos[np.where(RNAPs_strand == 1)]+=1
        RNAPs_pos[np.where(RNAPs_strand == -1)]-=1


        # Update the Dom_size (+1 or -1)
        Dom_size = np.ediff1d(Barr_pos)
        Dom_size = np.append(Dom_size, genome-Barr_pos[-1]+Barr_pos[0])

        # UPDATE SIGMA
        # R_plus_pos : the ids of RNA pol in the + strand
        R_plus_pos = np.where(Barr_type == 1)[0].astype(int)
        # R_minus_pos : the ids of RNA pol in the - strand
        R_minus_pos = np.where(Barr_type == -1)[0].astype(int)

        #### Extract all types of domaines (Those are ids of domaines)
        # Barr_type_ahead to make the extraction circular ;)
        Barr_type_ahead = np.roll(Barr_type, -1)
        # O +
        Barr_Dom_RPlus = np.where((Barr_type==0) & (Barr_type_ahead==1))
        # O -
        Barr_Dom_RMinus = np.where((Barr_type==0) & (Barr_type_ahead==-1))
        # O O
        Barr_Dom_Barr = np.where((Barr_type==0) & (Barr_type_ahead==0))
        # + +
        RPlus_Dom_RPlus = np.where((Barr_type==1) & (Barr_type_ahead==1))
        # - -
        RMinus_Dom_RMinus = np.where((Barr_type==-1) & (Barr_type_ahead==-1))
        # + -
        RPlus_Dom_RMinus = np.where((Barr_type==1) & (Barr_type_ahead==-1))
        # - +
        RMinus_Dom_RPlus = np.where((Barr_type==-1) & (Barr_type_ahead==+1))
        # - O
        RMinus_Dom_Barr = np.where((Barr_type==-1) & (Barr_type_ahead==0))
        # + O
        RPlus_Dom_Barr = np.where((Barr_type_ahead==0) & (Barr_type==+1))

        #### And then correct the value of Sigma in each case (before/after)
        corr_sig_Barr_Dom_RPlus = (Dom_size[Barr_Dom_RPlus]-1)/(Dom_size[Barr_Dom_RPlus]) # Sigma decrease x1
        corr_sig_Barr_Dom_RMinus = (Dom_size[Barr_Dom_RMinus]+1)/(Dom_size[Barr_Dom_RMinus]) # Sigma increase x1
        corr_sig_Barr_Dom_Barr = (Dom_size[Barr_Dom_Barr])/(Dom_size[Barr_Dom_Barr]) # Sigma FIX
        corr_sig_RPlus_Dom_RPlus = (Dom_size[RPlus_Dom_RPlus])/(Dom_size[RPlus_Dom_RPlus]) # Sigma FIX
        corr_sig_RMinus_Dom_RMinus = (Dom_size[RMinus_Dom_RMinus])/(Dom_size[RMinus_Dom_RMinus]) # Sigma FIX
        corr_sig_RPlus_Dom_RMinus = (Dom_size[RPlus_Dom_RMinus]+2)/(Dom_size[RPlus_Dom_RMinus]) # Sigma increase x2
        corr_sig_RMinus_Dom_RPlus = (Dom_size[RMinus_Dom_RPlus]-2)/(Dom_size[RMinus_Dom_RPlus]) # Sigma decrease x2
        corr_sig_RMinus_Dom_Barr = (Dom_size[RMinus_Dom_Barr]-1)/(Dom_size[RMinus_Dom_Barr]) # Sigma decrease x1
        corr_sig_RPlus_Dom_Barr = (Dom_size[RPlus_Dom_Barr]+1)/(Dom_size[RPlus_Dom_Barr]) # Sigma increase x1

        ### Multiply Sigma *= Corr (Each sigma value correspond to an specific domaine)
        Barr_sigma[Barr_Dom_RPlus] *= corr_sig_Barr_Dom_RPlus
        Barr_sigma[Barr_Dom_RMinus] *= corr_sig_Barr_Dom_RMinus
        Barr_sigma[Barr_Dom_Barr] *= corr_sig_Barr_Dom_Barr
        Barr_sigma[RPlus_Dom_RPlus] *= corr_sig_RPlus_Dom_RPlus
        Barr_sigma[RMinus_Dom_RMinus] *= corr_sig_RMinus_Dom_RMinus
        Barr_sigma[RPlus_Dom_RMinus] *= corr_sig_RPlus_Dom_RMinus
        Barr_sigma[RMinus_Dom_RPlus] *= corr_sig_RMinus_Dom_RPlus
        Barr_sigma[RMinus_Dom_Barr] *= corr_sig_RMinus_Dom_Barr
        Barr_sigma[RPlus_Dom_Barr] *= corr_sig_RPlus_Dom_Barr

        ### Now calculate the SC generated in each domaine
        # RNAPs_genSC_all : contains an array of RNAPs_genSC that should be added or substracted from each domaine
        RNAPs_genSC_all = RNAPs_genSC/Dom_size

        # Now update the value of sigma
        Barr_sigma[Barr_Dom_RPlus] -= RNAPs_genSC_all[Barr_Dom_RPlus]
        Barr_sigma[Barr_Dom_RMinus] += RNAPs_genSC_all[Barr_Dom_RMinus]
        Barr_sigma[RPlus_Dom_RMinus] += 2*RNAPs_genSC_all[RPlus_Dom_RMinus]
        Barr_sigma[RMinus_Dom_RPlus] -= 2*RNAPs_genSC_all[RMinus_Dom_RPlus]
        Barr_sigma[RMinus_Dom_Barr] -= RNAPs_genSC_all[RMinus_Dom_Barr]
        Barr_sigma[RPlus_Dom_Barr] += RNAPs_genSC_all[RPlus_Dom_Barr]


        # Now calc_sigma
        Barr_sigma = calc_sigma(Barr_sigma, p['GYRASE_CONC'], p['k_GYRASE'], p['x0_GYRASE'], \
        p['GYRASE_CTE'], p['TOPO_CONC'], p['k_TOPO'], p['x0_TOPO'], p['TOPO_CTE'], p['DELTA_T'])

        mean_sig_wholeGenome = np.sum(Barr_sigma*Dom_size)/genome

        # Update the initiation rate
        init_rate = f_init_rate(tr_rate, sigma_tr_start, p['sigma_t'], p['epsilon'], p['m'])


    return(tr_nbr/sum(tr_nbr))


###########################################################
#                          GO !                           #
###########################################################
F,e,g=simulation()
print("\nList of events :")
print(e)
print("\nList of fitnesses :")
print(F)
print("\nList of genes_level :")
print(g)

