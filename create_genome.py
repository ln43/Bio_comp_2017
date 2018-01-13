import numpy as np
import os

num_indiv=input("Entrer le numero du genome a creer : ")

Barr_fix=np.arange(1,33000,3000)
TSS_pos=np.zeros(10)

TSS_pos=[np.random.randint(Barr_fix[i]+500,Barr_fix[i]+1500) for i in range(10)]
TTS_pos=[np.random.randint(TSS_pos[i]+200,Barr_fix[i+1]-500) for i in range(10)]    
strands=np.random.choice([0,1],10,p=[0.2,0.8])
signe=["-","+"]

os.makedirs("INDIV_"+str(num_indiv), exist_ok=True)
os.chdir("INDIV_"+str(num_indiv))

# WRITE GFF
f = open("indiv"+str(num_indiv)+".gff", "a")
f.write("##gff-version 3\n")
f.write("#!gff-spec-version 1.20\n")
f.write("#!processor NCBI annotwriter\n")
f.write("##sequence-region tousgenesidentiques 1 30000\n")
f.write("indiv"+str(num_indiv)+"	RefSeq	region	1	30000	.	+	.	ID=id1;Name=newindiv\n")

for i in range(10) :
    f.write("indiv"+str(num_indiv)+"	RefSeq	gene	"+str(TSS_pos[i])+"	"+str(TTS_pos[i])+"	.	"+signe[strands[i]]+"	.	IS=g"+str(i+1)+";Name=g"+str(i+1)+"\n")
    
f.close()

# WRITE prot.dat
f = open("prot"+str(num_indiv)+".dat", "a")
f.write("prot_name	prot_pos\n")
for i in range(10) :
    f.write("hns	"+str(Barr_fix[i])+"\n")

f.close()

# WRITE TSS.dat
f = open("TSS"+str(num_indiv)+".dat", "a")
f.write("TUindex	TUorient	TSS_pos	TSS_strength\n")
for i in range(10) :
    f.write(str(i)+"	"+signe[strands[i]]+"	"+str(TSS_pos[i])+"	"+".02"+"\n")
f.close()

# WRITE TTS.dat
f = open("TTS"+str(num_indiv)+".dat", "a")
f.write("TUindex	TUorient	TTS_pos	TTS_proba_off\n")
for i in range(10) :
    f.write(str(i)+"	"+signe[strands[i]]+"	"+str(TTS_pos[i])+"	"+"1."+"\n")
f.close()


os.chdir("..")