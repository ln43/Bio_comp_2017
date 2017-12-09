import numpy as np

TSS_pos=np.array([ 16,  66, 116, 166, 216, 266, 316, 366, 416, 466])
strands=np.array([1, 1, 1, 1, 1, 1, 1, 1, 1, 1], dtype=object)

TTS_pos=np.array([ 33,  83, 133, 183, 233, 283, 333, 383, 433, 483])

genes=np.array(["g1","g2","g3","g4","g5","g6","g7","g8","g9","g10"])

# TTS_pos  = (tts['TTS_pos'].values/DELTA_X).astype(int)

genome=499

# Get point to invert not in coding part
inv1=np.random.randint(0,genome)
while set(np.ndarray.tolist(np.where(inv1>TSS_pos)[0])).intersection(np.ndarray.tolist(np.where(inv1<TTS_pos)[0])):
    inv1=np.random.randint(0,genome)
inv2=np.random.randint(0,genome)
while set(np.ndarray.tolist(np.where(inv2>TSS_pos)[0])).intersection(np.ndarray.tolist(np.where(inv2<TTS_pos)[0])):
    inv2=np.random.randint(0,genome)   

    
if inv1>inv2:
    S1=np.ndarray.tolist(np.where(inv2>TTS_pos)[0])
    S2=np.ndarray.tolist(np.where(inv1<TSS_pos)[0])
    S=S2+S1
    
    newS1=S[0:len(S1)]
    newS2=S[len(S1):]
    newS1.reverse()
    newS2.reverse()
    
    newgenes=np.copy(genes)
    newgenes[S1]=genes[newS1]
    newgenes[S2]=genes[newS2]
    genes=newgenes
    
    newstrands=np.copy(strands)
    newstrands[S1]=-strands[newS1]
    newstrands[S2]=-strands[newS2]
    strands=newstrands
    
else:
    S=list(set(np.ndarray.tolist(np.where(inv1<TSS_pos)[0])).intersection(np.ndarray.tolist(np.where(inv2>TTS_pos)[0])))
    S2=np.copy(S)
    S.reverse()
    genes[S2]=genes[S]
    strands[S2]=-strands[S]

