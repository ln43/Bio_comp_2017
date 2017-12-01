import sys
import simulation as sim
import numpy as np

# INI_file=sys.argv[1]
# output_dir=sys.argv[2]

INI_file="params.ini"
output_dir="OUTPUTDIR"

sim.start_transcribing(INI_file, output_dir)



TSS_pos=np.array([ 16,  66, 116, 166, 216, 266, 316, 366, 416, 466])
strands=np.array([1, 1, 1, 1, 1, 1, 1, 1, 1, 1], dtype=object)
TTS_pos=np.array([ 33,  83, 133, 183, 233, 283, 333, 383, 433, 483])
genome=499

def inversion(TSS_pos,TTS_pos,strands,genome):
    inv1=np.random.randint(0,genome)
    while set(np.ndarray.tolist(np.where(inv1>TSS_pos)[0])).intersection(np.ndarray.tolist(np.where(inv1<TTS_pos)[0])):
        inv1=np.random.randint(0,genome)
        
    inv2=np.random.randint(0,genome)
    while set(np.ndarray.tolist(np.where(inv2>TSS_pos)[0])).intersection(np.ndarray.tolist(np.where(inv2<TTS_pos)[0])):
        inv2=np.random.randint(0,genome)   
        
    if inv1>inv2:
        S=set(np.ndarray.tolist(np.where(inv1<TSS_pos)[0])).union(np.ndarray.tolist(np.where(inv2>TTS_pos)[0]))
    else:
        S=set(np.ndarray.tolist(np.where(inv1<TSS_pos)[0])).intersection(np.ndarray.tolist(np.where(inv2>TTS_pos)[0]))
    
    strands[list(S)]=-strands[list(S)]
    return strands
    
strands=inversion(TSS_pos,TTS_pos,strands,genome)
print(strands)


################## STRUCTURE

# Indiv 1 (tousgenesidentiques)
# Simulation indiv1
# 
# for i in range(1,niter):
#     prob_inv
#     if prob_inv :
#         -> inversion indiv 1
#     prob_insertdel
#     if prob_insertdel :
#         -> insertion/deletion indiv 1

# Si pas de inversion ni insert/del, passer direct à l'itération suivante !

#     -> Ecritures fichier individus 2
#     -> Simulation indiv 2
#     -> calcul fitness
#     if meilleur fitness
#         indiv 1 = indiv 2
        
        
