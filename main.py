import os
import class_individu
import configparser
import pandas as pd
import numpy as np
import collections as col
from pylab import *
import errno
import csv


RNAPs_genSC = 0.1
###########################################################
#                       Functions                         #
###########################################################

def create_config_file(config, config_file_path, TSS_file, SIGMA_0, RNAPS_NB): # DELTA_X, D, J_0, SIGMA_0, RNAPS_NB,
    # Create the directory
    # output_dir = "D_%d/delta_x_%d" %(D, DELTA_X)
    #os.makedirs(output_dir, exist_ok=True)
    # Create the config file
    config.set('INPUTS','TSS', str(TSS_file))
    #config.set('GLOBAL','D', str(D))
    #config.set('GLOBAL', 'J_0', str(J_0))
    #config.set('SIMULATION','DELTA_X', str(DELTA_X))
    config.set('SIMULATION','SIGMA_0', str(SIGMA_0))
    config.set('SIMULATION','RNAPS_NB', str(RNAPS_NB))
    with open(config_file_path, 'w') as configfile:
        config.write(configfile)

# Read the config files
def read_config_file(path):
    config = configparser.ConfigParser()
    # to preserve capital letters
    config.optionxform = str
    config.read(path)
    return config

# Read the config files and return the values of each variable
# this function will be useful when we are in another script
def read_config_file_v2(path):
    config = configparser.ConfigParser()
    # to preserve capital letters
    config.optionxform = str
    config.read(path)

    # get inputs infos from the config file
    GFF_file = config.get('INPUTS', 'GFF')
    TSS_file = config.get('INPUTS', 'TSS')
    TTS_file = config.get('INPUTS', 'TTS')
    Prot_file = config.get('INPUTS', 'BARR_FIX')

    # get values from the config file
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
    #SIGMA_0 = 0 #((-np.log(((GYRASE_CONC*GYRASE_CTE)/TOPO_CONC*TOPO_CTE)-1))/k)+x_0
    #$print("SIGMA_0 --> ", SIGMA_0)

    return GFF_file, TSS_file, TTS_file, Prot_file, m, sigma_t, epsilon, SIGMA_0, DELTA_X, DELTA_T, RNAPS_NB, ITERATIONS_NB, OUTPUT_STEP, GYRASE_CONC, TOPO_CONC, TOPO_CTE, GYRASE_CTE, TOPO_EFFICIENCY, k_GYRASE, x0_GYRASE, k_TOPO, x0_TOPO


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
            k = 0 #TU_id #0 #len(this_TU_tts)# TTS id index ### [0 1 2 3 4 5]
            proba_rest = 1
            while proba_rest > 0 and k < len(this_TU_tts) : # >= 0 : #
                if tts['TTS_pos'][k] < tss['TSS_pos'][i] : # and tts['TUindex'][k] == TU_id :
                    tr_id.append(j)
                    tr_strand.append(-1)
                    tr_start.append(tss['TSS_pos'][i])
                    # after getting them, we shall (in every loop) generate a new tr_end
                    tr_end.append(this_TU_tts[k])
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

###################### Saving files #######################

def save_files(output_dir,
                Barr_pos, Barr_type, Dom_size, Barr_ts_remain, Barr_sigma,
                tr_nbr, tr_times, save_RNAPs_info, save_tr_info,
                save_Barr_sigma, save_Dom_size, save_mean_sig_wholeGenome,
                DELTA_X, RNAPs_genSC,
                RNAPs_tr, RNAPs_pos, RNAPs_unhooked_id,
                init_rate, Kon, RNAPS_NB, SIGMA_0, GYRASE_CONC, TOPO_CONC):

    if RNAPs_genSC != 0:
        output_dir += "/withSC_Kon_%.06f/RNAPn_%s/Sig0_%s/Gyrase_%s_TopoI_%s/" %(Kon[0], RNAPS_NB, SIGMA_0, GYRASE_CONC, TOPO_CONC)
    else:
        output_dir += "/withoutSC_Kon_%.06f/RNAPn_%s/Sig0_%s/Gyrase_%s_TopoI_%s/" %(Kon[0], RNAPS_NB, SIGMA_0, GYRASE_CONC, TOPO_CONC)

    # make sure that the output direcory exists, and create one if it's not
    os.makedirs("%s/resume_sim" %output_dir, exist_ok=True)
    os.makedirs("%s/all_res" %output_dir, exist_ok=True)

    # save tr_nbr
    tr_nbr = pd.Series(tr_nbr)
    tr_nbr.to_csv("%s/save_tr_nbr.csv" %output_dir, sep=';', index=False) #, header = ["Diffusion = " + str(D)])

    # convert tr_times dict to pandas serie
    tr_times = pd.DataFrame.from_dict(tr_times, orient='index') #pd.DataFrame([tr_times]) #pd.Series(tr_times)
    # save the tr_times to csv file
    tr_times.to_csv("%s/save_tr_times.csv" %output_dir, sep=';', index=True, header = False)

    # Save last info
    np.savez("%s/resume_sim/resume_sim_RNAPs.npz" %output_dir, RNAPs_tr = RNAPs_tr,
                                                               RNAPs_pos = RNAPs_pos,
                                                               RNAPs_unhooked_id = RNAPs_unhooked_id)

    np.savez("%s/resume_sim/resume_sim_tr.npz" %output_dir, tr_nbr = tr_nbr,
                                                            init_rate = init_rate)

    np.savez("%s/resume_sim/resume_sim_Barr.npz" %output_dir, Barr_pos = Barr_pos,
                                                              Barr_type = Barr_type,
                                                              Dom_size = Dom_size,
                                                              Barr_ts_remain = Barr_ts_remain,
                                                              Barr_sigma = Barr_sigma)

    # Save all info
    np.savez("%s/all_res/save_RNAPs_info" %output_dir, RNAPs_info = save_RNAPs_info)
    np.savez("%s/all_res/save_tr_info" %output_dir, tr_info = save_tr_info)
    np.savez("%s/all_res/save_sigma_info" %output_dir, Barr_sigma_info = save_Barr_sigma, Dom_size_info = save_Dom_size, mean_sig_wholeGenome = save_mean_sig_wholeGenome)

###########################################################
#         Transcription Process (Simulation)              #
###########################################################

def simulation():
    config = read_config_file(INI_file)

    # get inputs infos from the config file
    GFF_file = config.get('INPUTS', 'GFF')
    TSS_file = config.get('INPUTS', 'TSS')
    TTS_file = config.get('INPUTS', 'TTS')
    Prot_file = config.get('INPUTS', 'BARR_FIX')

    # get values from the config file
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
    #SIGMA_0 = 0 #((-np.log(((GYRASE_CONC*GYRASE_CTE)/TOPO_CONC*TOPO_CTE)-1))/k)+x_0
    #$print("SIGMA_0 --> ", SIGMA_0)

    # define the output directory
    os.makedirs(output_dir, exist_ok=True)

    # path to the input files (remove the "params.ini" from the path)
    pth = INI_file[:-10]
    gff_df_raw = load_gff(pth+GFF_file)
    tss = load_tab_file(pth+TSS_file)
    tts = load_tab_file(pth+TTS_file)
    prot = load_tab_file(pth+Prot_file)
    genome_size = get_genome_size(gff_df_raw)
    gff_df = rename_gff_cols(gff_df_raw)
    TU_tts =  get_TU_tts(tss, tts)
    params = [sigma_t, epsilon, SIGMA_0, DELTA_X, DELTA_T, RNAPS_NB, ITERATIONS_NB, OUTPUT_STEP, GYRASE_CONC, TOPO_CONC, TOPO_CTE, GYRASE_CTE, TOPO_EFFICIENCY, k_GYRASE, x0_GYRASE, k_TOPO, x0_TOPO, m]
    i = individu(gff_df_raw, tss, tts, prot, TU_tts, params)
    #faire calculs
    for i in range(0, nb_iter):
        #faire des trucs
