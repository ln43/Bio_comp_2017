class individu():
	def __init__(self, gff_df, tss, tts, prot, genome_size, DELTA_X, genes):
		self.tts = tts
		self.tss = tss
		self.gff_df_raw=gff_df
		self.prot = prot

		self.TTS_pos = (tts['TTS_pos'].values/DELTA_X).astype(int)
		self.newTTS_pos = np.copy(self.TTS_pos)

		self.genome_size = genome_size
		self.genome = int(genome_size/DELTA_X)
		self.newgenome = genome

		self.TSS_pos = (tss['TSS_pos'].values/DELTA_X).astype(int)
		self.newTSS_pos = np.copy(self.TTS_pos)

		self.Barr_fix = (prot['prot_pos'].values/DELTA_X).astype(int)
		self.newBarr_fix = np.copy(Barr_fix)

		self.strands = self.str2num(gff_df['strand'].values)
		self.newstrands = np.copy(strands)
		self.noms_genes = range(1, len(genes) + 1)
		self.genes = genes
		self.newgenes = np.copy(self.genes)
		self.fitness = self.calcul_fitness(genes)
		self.new_fitness = 0

	def str2num(self, s):
		s[s == '+'] = 1 #True
		s[s == '-'] = -1 #False
		return s

	def inversion(self):
		inv1=np.random.randint(0,self.genome)
		while set(np.ndarray.tolist(np.where(inv1>self.TSS_pos)[0])).intersection(np.ndarray.tolist(np.where(inv1<self.TTS_pos)[0])): # test if inv1 belongs to non-coding part
			inv1=np.random.randint(0,self.genome)

		inv2=np.random.randint(0,self.genome)
		while set(np.ndarray.tolist(np.where(inv2>self.TSS_pos)[0])).intersection(np.ndarray.tolist(np.where(inv2<self.TTS_pos)[0])):  # test if inv2 belongs to non-coding part
			inv2=np.random.randint(0,self.genome)

		if inv1>inv2: # Invert extern part
			S1=np.ndarray.tolist(np.where(inv2>self.TTS_pos)[0])
			S2=np.ndarray.tolist(np.where(inv1<self.TSS_pos)[0])
			S=S2+S1
			
			newS1=S[0:len(S1)]
			newS2=S[len(S1):]
			newS1.reverse()
			newS2.reverse()
			
			self.newgenes=np.copy(self.genes)
			self.newgenes[S1]=self.genes[newS1]
			self.newgenes[S2]=self.genes[newS2]
			
			self.newstrands=np.copy(self.strands)
			self.newstrands[S1]=-self.strands[newS1]
			self.newstrands[S2]=-self.strands[newS2]
			
		else: # Invert intern part
			S=list(set(np.ndarray.tolist(np.where(inv1<self.TSS_pos)[0])).intersection(np.ndarray.tolist(np.where(inv2>self.TTS_pos)[0])))
			S2=np.copy(S)
			S.reverse()
			self.newgenes[S2]=self.genes[S]
			self.newstrands[S2]=-self.strands[S]
		
	def indel(self) :
		ind=np.random.randint(0,self.genome)
		while set(np.ndarray.tolist(np.where(ind>self.TSS_pos)[0])).intersection(np.ndarray.tolist(np.where(ind<self.TTS_pos)[0])): # test if ind belongs to non-coding part
			ind=np.random.randint(0,self.genome)
		prob=np.random.rand(1)
		if prob<0.5 : # Insertion
			print("insert")
			self.newTSS_pos[np.where(ind<self.TSS_pos)[0]]+=1
			self.newTTS_pos[np.where(ind<self.TTS_pos)[0]]+=1
			self.newBarr_fix[np.where(ind<self.Barr_fix)[0]]+=1
			self.newgenome+=1
		else : # Deletion
			print("deletion")
			self.newTSS_pos[np.where(ind<self.TSS_pos)[0]]-=1
			self.newTTS_pos[np.where(ind<self.TTS_pos)[0]]-=1
			self.newBarr_fix[np.where(ind<self.Barr_fix)[0]]-=1
			self.newgenome-=1

