class individu():
	def __init__(self, gff_df, tss, tts, prot, genome_size, DELTA_X):
#params : liste trucs csts
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
	
		if inv1>inv2:
			S=set(np.ndarray.tolist(np.where(inv1<self.TSS_pos)[0])).union(np.ndarray.tolist(np.where(inv2>self.TTS_pos)[0]))
		else:
			S=set(np.ndarray.tolist(np.where(inv1<self.TSS_pos)[0])).intersection(np.ndarray.tolist(np.where(inv2>self.TTS_pos)[0]))
	
		self.strands[list(S)]=-self.strands[list(S)]

