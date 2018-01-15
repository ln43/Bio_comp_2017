import numpy as np
import pandas as pd

class individu():
	def __init__(self, gff_df, tss, tts, prot, genome_size, DELTA_X, genes,p_keep):
		# Data frame, contains all informations from tss.dat, tts.dat
		self.tts = tts
		self.newtts = pd.DataFrame.copy(tts)

		self.tss = tss
		self.newtss = pd.DataFrame.copy(tss)

		self.DELTA_X=int(DELTA_X) #Spatial discretization in bp


		# Genome Size
		self.genome_size = genome_size
		self.newgenome = self.genome_size

		# Position Vectors at time t and time t+1 (new...)
		self.TSS_pos = tss['TSS_pos'].values
		self.newTSS_pos = np.copy(self.TSS_pos)

		self.TTS_pos = tts['TTS_pos'].values
		self.newTTS_pos = np.copy(self.TTS_pos)

		self.Barr_fix = prot['prot_pos'].values
		self.newBarr_fix = np.copy(self.Barr_fix)

		self.strands = self.str2num(gff_df['strand'].values)
		self.newstrands = np.copy(self.strands)
		
		self.noms_genes = np.arange(1, len(genes) + 1)
		self.newnoms_genes = np.copy(self.noms_genes)
		
		
		self.genes_level_envir = genes # Ideal gene expression profil
		
		self.fitness = 0
		self.new_fitness = 0

		self.pkeep=p_keep # tolerated difference between fitness(t+1) and fitness(t) if fitness(t+1)>fitness(t)

	def str2num(self, s):
		s[s == '+'] = 1 #True
		s[s == '-'] = -1 #False
		return s
		
	# Make inversion of one part of the genome
	def inversion(self):
		# Invert part of chromosome
		
		# Discretize positions
		TTS_pos_disc=(np.copy(self.TTS_pos)/self.DELTA_X).astype(int)
		TSS_pos_disc=(np.copy(self.TSS_pos)/self.DELTA_X).astype(int)
		Barr_fix_disc=(np.copy(self.Barr_fix)/self.DELTA_X).astype(int)
		genome_disc=int(self.genome_size/self.DELTA_X)

		# Copy position in new vector to modify at time t+1		
		newTTS_pos_disc=np.copy(TTS_pos_disc)
		newTSS_pos_disc=np.copy(TSS_pos_disc)
		newBarr_fix_disc=np.copy(Barr_fix_disc)

		# Extract position of coding part
		coding=[]
		for i in range(len(TTS_pos_disc)):
			coding.extend(range(min(TTS_pos_disc[i],TSS_pos_disc[i]),max(TTS_pos_disc[i],TSS_pos_disc[i])))
		coding=np.unique([coding[i]+2 for i in range(len(coding))]+coding+[coding[i]-2 for i in range(len(coding))]) # avoid position too closed to coding part (genes must be separated by at least 2*DELTA_X bases

		# Get inversion positions
		inv1=np.random.randint(0,genome_disc)
		while inv1 in coding: # test if inv1 belongs to non-coding part
			inv1=np.random.randint(0,genome_disc)
	
		inv2=np.random.randint(0,genome_disc)
		while inv2 in coding:  # test if inv2 belongs to non-coding part
			inv2=np.random.randint(0,genome_disc)
	
		# Get position min and max between TSS and TTS for each gene
		T_temp=np.array([TTS_pos_disc,TSS_pos_disc])
		Tmax=T_temp.max(axis=0) # = TTS if strand='+' and TSS if strands = '-'
		Tmin=T_temp.min(axis=0) # = TTS if strand='-' and TSS if strands = '+'
	
	
		if inv1>inv2: # Invert extern part (extremities)
	
			# Get index of genes to invert
			S1=np.ndarray.tolist(np.where(inv2>Tmax)[0])
			S2=np.ndarray.tolist(np.where(inv1<Tmin)[0])
			S=S2+S1
	
			if len(S1)>0 or len(S2)>0 : # if gene exists in selected part
				
				# Get new genes index when inverted
				newS1=S[0:len(S1)]
				newS2=S[len(S1):]
				newS1.reverse()
				newS2.reverse()
	
				# Invert name
				self.newnoms_genes=np.copy(self.noms_genes)
				self.newnoms_genes[S1]=self.noms_genes[newS1]
				self.newnoms_genes[S2]=self.noms_genes[newS2]
	
				# Invert strands
				self.newstrands=np.copy(self.strands)
				self.newstrands[S1]=-self.strands[newS1]
				self.newstrands[S2]=-self.strands[newS2]
	
				# Invert TSS_pos and TTS_pos
				if len(S1)>0 and len(S2)>0 : # gene on both extremities
					ecart2=inv1+inv2-Tmax[S1[-1]]-Tmin[S2[0]]
					newTTS_pos_disc[S2]=TTS_pos_disc[S2] + ecart2
					newTSS_pos_disc[S2]=TSS_pos_disc[S2] + ecart2
	
					newTSS_pos_disc[S1]=TSS_pos_disc[S1] + ecart2
					newTTS_pos_disc[S1]=TTS_pos_disc[S1] + ecart2
				else :
					if len(S1)>0: # gene at begginning only
						ecart1=Tmax[S1[-1]]-(inv2-(genome_disc-inv1+Tmin[0]))
						newTSS_pos_disc[S1]=TSS_pos_disc[S1]-ecart1
						newTTS_pos_disc[S1]=TTS_pos_disc[S1]-ecart1
					elif len(S2)>0: # gene at the end only
						ecart2=inv1+inv2+genome_disc-Tmax[S2[-1]]-Tmin[S2[0]]
						newTTS_pos_disc[S2]=TTS_pos_disc[S2] + ecart2
						newTSS_pos_disc[S2]=TSS_pos_disc[S2] + ecart2
				
				# Translate all genes not to have genes cut
				translation=int((genome_disc+min(newTSS_pos_disc[0],newTTS_pos_disc[0])-max(newTTS_pos_disc[-1],newTSS_pos_disc[-1]))/2) - min(newTSS_pos_disc[0],newTTS_pos_disc[0])
				newTSS_pos_disc=newTSS_pos_disc+translation+1
				newTTS_pos_disc=newTTS_pos_disc+translation+1
	
	
		else: # Invert intern part
		
			# Get index of genes to invert
			S=list(set(np.ndarray.tolist(np.where(inv1<Tmin)[0])).intersection(np.ndarray.tolist(np.where(inv2>Tmax)[0])))
			S.sort()
			
			if len(S)!=0: # if gene exists in selected part
			
				# Get new genes index when inverted
				S2=np.copy(S)
				S.reverse()
	
				# Invert name
				self.newnoms_genes[S2]=self.noms_genes[S]
				
				# Invert strands
				self.newstrands[S2]=-self.strands[S]
	
				# Invert TSS_pos and TTS_pos
				ecart=inv1+inv2-Tmax[S2[-1]]-Tmin[S2[0]]
				newTSS_pos_disc[S2] = TSS_pos_disc[S2]+ecart
				newTTS_pos_disc[S2] = TTS_pos_disc[S2]+ecart

		# Set TSS at right position and TTS at left position if strands = '+' and conversely
		newT_temp=np.array([newTTS_pos_disc,newTSS_pos_disc])
		newTmax=newT_temp.max(axis=0)
		newTmin=newT_temp.min(axis=0)
		newTTS_pos_disc[self.newstrands<0]=newTmin[self.newstrands<0]
		newTTS_pos_disc[self.newstrands>0]=newTmax[self.newstrands>0]
		newTSS_pos_disc[self.newstrands<0]=newTmax[self.newstrands<0]
		newTSS_pos_disc[self.newstrands>0]=newTmin[self.newstrands>0]

		# Create new barriers between two genes (at equal distance)
		newBarr_fix_disc=[0]
		for i in range(0,len(self.Barr_fix)-1):
				newBarr_fix_disc.append(int((min(newTSS_pos_disc[i+1],newTTS_pos_disc[i+1]) +max(newTSS_pos_disc[i],newTTS_pos_disc[i]))/2))
		newBarr_fix_disc= np.asarray(newBarr_fix_disc)
		
		# Convert position in true unit 
		self.newBarr_fix=np.copy(newBarr_fix_disc)*self.DELTA_X
		self.newTSS_pos=np.copy(newTSS_pos_disc)*self.DELTA_X
		self.newTTS_pos=np.copy(newTTS_pos_disc)*self.DELTA_X

	# Make Insertion or Deletion of one part of the genome
	def indel(self) :
		
		# Extract position of coding base
		coding=[]
		for i in range(len(self.TTS_pos)):
			coding.extend(range(min(self.TTS_pos[i],self.TSS_pos[i]),max(self.TTS_pos[i],self.TSS_pos[i])))

		prob=np.random.rand(1) # Proba to choose between insertion and deletion
		
		if prob<0.5 : # Insertion
		
			# Position on Insertion
			ind=np.random.randint(0,self.genome_size)
			while ind in coding: # test if ind belongs to non-coding part
				ind=np.random.randint(0,self.genome_size)
			
			# Set new positions for TSS, TTS and Barr_fix (translation on left part)
			self.newTSS_pos[np.where(ind<self.newTSS_pos)[0]]+=self.DELTA_X
			self.newTTS_pos[np.where(ind<self.newTTS_pos)[0]]+=self.DELTA_X
			self.newBarr_fix[np.where(ind<self.newBarr_fix)[0]]+=self.DELTA_X
			
			self.newgenome+=self.DELTA_X # Set new genome size
			
			return(1) # Index of 'Insertion' for the final plot
			
		else : # Deletion
			coding=np.unique(np.concatenate((coding,np.arange(0+2*self.DELTA_X))))
			coding=np.unique([coding[i]+2*self.DELTA_X for i in range(len(coding))]+coding) # Avoid to have gene which are separated with less than 2*DELTA_X
			
			# Position on Deletion
			ind=np.random.randint(0,self.genome_size)
			while ind in coding: # test if ind belongs to non-coding part
				ind=np.random.randint(0,self.genome_size)

			# Set new positions for TSS, TTS and Barr_fix (translation on right part)
			self.newTSS_pos[np.where(ind<self.newTSS_pos)[0]]-=self.DELTA_X
			self.newTTS_pos[np.where(ind<self.newTTS_pos)[0]]-=self.DELTA_X
			#self.newBarr_fix[np.where(ind<self.newBarr_fix)[0]]-=self.DELTA_X
			
			newBarr_fix=[0]
			for i in range(0,len(self.Barr_fix)-1):
				newBarr_fix.append(int((min(self.newTSS_pos[i+1],self.newTTS_pos[i+1]) +max(self.newTSS_pos[i],self.newTTS_pos[i]))/2))
			self.newBarr_fix= np.asarray(newBarr_fix)
			
			self.newgenome-=self.DELTA_X # Set new genome size
			
			return(0) # Index of 'Deletion' for the final plot

	# Modify data frame argument after inversion or insertion/deletion
	def upgrade_new_pd_dataframe(self) :
		
		self.newtss.TSS_pos=self.newTSS_pos
		self.newtts.TTS_pos=self.newTTS_pos
		signes={1:"+",-1:"-"}
		strands_sign=[]
		for i in self.newstrands :
			strands_sign.append(signes[i])
		self.newtss.TUorient=strands_sign
		self.newtts.TUorient=strands_sign


	# Calculate new fitness comparing genes expression profil to the ideal one
	def calcul_fitness(self,genes_level) :
		return (sum(abs((self.genes_level_envir-genes_level[np.argsort(self.newnoms_genes)]))/genes_level[np.argsort(self.newnoms_genes)]))

	def update_fitness(self,genes_level) :
		self.new_fitness=self.calcul_fitness(genes_level)
		
	# Choose to keep or not the new individual created
	# We want to minimize fitness with a probability to keep the wrong genome anyway
	def choice_indiv(self) :
		
		if self.fitness<self.new_fitness :
			p = int(self.new_fitness-self.fitness>=self.pkeep) # if new fitness is worst, we will keep new genome if the difference between fitnesses at time t and t+1 is small (<pkeep)
		else :
			p=0

		if p==1 : # keep the old genome
			self.newTTS_pos = np.copy(self.TTS_pos)
			self.newgenome = self.genome_size
			self.newTSS_pos = np.copy(self.TSS_pos)
			self.newBarr_fix = np.copy(self.Barr_fix)
			self.newstrands = np.copy(self.strands)
			self.newnoms_genes = np.copy(self.noms_genes)
			self.new_fitness = self.fitness
		else: # keep the new genome, even if it doesn't improve the fitness
			self.TTS_pos = np.copy(self.newTTS_pos)
			self.genome_size = self.newgenome
			self.TSS_pos = np.copy(self.newTSS_pos)
			self.Barr_fix = np.copy(self.newBarr_fix)
			self.strands = np.copy(self.newstrands)
			self.noms_genes = np.copy(self.newnoms_genes)
			self.fitness = self.new_fitness
