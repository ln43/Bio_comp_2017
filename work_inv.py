import numpy as np

TSS_pos=np.array([ 16,  66, 116, 166, 216, 266, 316, 366, 416, 466])
strands=np.array([1, 1, 1, 1, 1, 1, 1, 1, 1, 1], dtype=object)

TTS_pos=np.array([ 33,  83, 133, 183, 233, 283, 333, 383, 433, 483])

# TTS_pos  = (tts['TTS_pos'].values/DELTA_X).astype(int)


Barr_fix=np.array([  0,  50, 100, 150, 200, 250, 300, 350, 400, 450])
genome=499

# Get point to invert not in coding part
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