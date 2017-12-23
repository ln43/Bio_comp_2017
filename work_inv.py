import numpy as np
import matplotlib.pyplot as plt

TSS_pos=np.array([ 16,  66, 116, 166, 216, 266, 316, 366, 416, 466])
newTSS_pos = np.copy(TSS_pos)

strands=np.array([1, 1, 1, 1, 1, 1, 1, 1, 1, 1], dtype=object)
newstrands=np.copy(strands)


TTS_pos=np.array([ 33,  83, 133, 183, 233, 283, 333, 383, 433, 483])
newTTS_pos = np.copy(TTS_pos)

noms_genes=np.array(["g1","g2","g3","g4","g5","g6","g7","g8","g9","g10"])
newnoms_genes=np.copy(noms_genes)


Barr_fix=np.array([  0,  50, 100, 150, 200, 250, 300, 350, 400, 450])
newBarr_fix=np.copy(Barr_fix)

genome=499

# Get point to invert not in coding part
inv1=np.random.randint(0,genome)
while set(np.ndarray.tolist(np.where(inv1>=TSS_pos)[0])).intersection(np.ndarray.tolist(np.where(inv1<=TTS_pos)[0])):
    inv1=np.random.randint(0,genome)
inv2=np.random.randint(0,genome)
while set(np.ndarray.tolist(np.where(inv2>=TSS_pos)[0])).intersection(np.ndarray.tolist(np.where(inv2<=TTS_pos)[0])):
    inv2=np.random.randint(0,genome)   

    
if inv1>inv2:
    print("ext")
    S1=np.ndarray.tolist(np.where(inv2>TTS_pos)[0])
    S2=np.ndarray.tolist(np.where(inv1<TSS_pos)[0])
    S=S2+S1
    if len(S1)>0 or len(S2)>0 :
        newS1=S[0:len(S1)]
        newS2=S[len(S1):]
        newS1.reverse()
        newS2.reverse()
        
        newnoms_genes=np.copy(noms_genes)
        newnoms_genes[S1]=noms_genes[newS1]
        newnoms_genes[S2]=noms_genes[newS2]
    
        newstrands=np.copy(strands)
        newstrands[S1]=-strands[newS1]
        newstrands[S2]=-strands[newS2]
        
        if len(S1)>0 and len(S2)>0 :
            ecart2=inv1+inv2-TTS_pos[S1[-1]]-TSS_pos[S2[0]]
            newTTS_pos[S2]=TTS_pos[S2] + ecart2
            newTSS_pos[S2]=TSS_pos[S2] + ecart2
            
            ecart1=TTS_pos[S1[-1]]-(inv1+inv2-TSS_pos[S2[0]])
            newTSS_pos[S1]=TSS_pos[S1]-ecart1
            newTTS_pos[S1]=TTS_pos[S1]-ecart1
        else :
            if len(S1)>0:
                ecart1=TTS_pos[S1[len(S1)-1]]-(inv2-(genome-inv1+TSS_pos[0]))
                newTSS_pos[S1]=TSS_pos[S1]-ecart1
                newTTS_pos[S1]=TTS_pos[S1]-ecart1
            else :
                if len(S2)>0:
                    ecart2=inv1+inv2+genome-TTS_pos[S2[-1]]-TSS_pos[S2[0]]
                    newTTS_pos[S2]=TTS_pos[S2] + ecart2
                    newTSS_pos[S2]=TSS_pos[S2] + ecart2
                    
        translation=int((genome+newTSS_pos[0]-newTTS_pos[-1])/2) - newTSS_pos[0]
        newTSS_pos=newTSS_pos+translation
        newTTS_pos=newTTS_pos+translation
        newBarr_fix[0]=0
        newBarr_fix[1:]=newTTS_pos[0:-1]+(newTSS_pos[1:]-newTTS_pos[0:-1])/2+1
    
else:
    S=list(set(np.ndarray.tolist(np.where(inv1<TSS_pos)[0])).intersection(np.ndarray.tolist(np.where(inv2>TTS_pos)[0])))
    if len(S)>0 :
        S2=np.copy(S)
        S.reverse()
        newnoms_genes[S2]=noms_genes[S]
        newstrands[S2]=-strands[S]
    
        ecart=inv1+(inv2-TTS_pos[S2[-1]])-TSS_pos[S2[0]]
        newTSS_pos[S2] = TSS_pos[S2]+ecart
        newTTS_pos[S2] = TTS_pos[S2]+ecart
        newBarr_fix[0]=0
        newBarr_fix[1:]=newTTS_pos[0:-1]+(newTSS_pos[1:]-newTTS_pos[0:-1])/2+1
    
    
# PLOT GENOME
fig = plt.figure()   
ax1 = fig.add_subplot(211)
ax1.scatter(Barr_fix, np.ones(10), s=40, c='r')
ax1.scatter(TSS_pos, np.ones(10), s=40, c='g')
ax1.scatter(TTS_pos, np.ones(10), s=40, c='b')
plt.legend(['Barrieres','TSS','TTS'])
plt.ylim([-1,2])
for i, txt in enumerate(noms_genes) :
    ax1.annotate(txt, (TSS_pos[i]+5,0.5))
plt.title('Before inversion')

ax2 = fig.add_subplot(212)
ax2.scatter(newBarr_fix, np.ones(10), s=40, c='r')
ax2.scatter(newTSS_pos, np.ones(10), s=40, c='g')
ax2.scatter(newTTS_pos, np.ones(10), s=40, c='b')
plt.legend(['Barrieres','TSS','TTS'])
plt.ylim([-1,2])
for i, txt in enumerate(newnoms_genes) :
    ax2.annotate(txt, (newTSS_pos[i]+5,0.5))
plt.title('After inversion')

plt.show()
