class individu():
	def __init__(self, gff_df, tss, tts, prot, TU_tts, genome_size, DELTA_X):
#params : liste trucs csts
		self.tts = tts
		self.tss = tss
		self.TTS_pos = (tts['TTS_pos'].values/DELTA_X).astype(int)
		self.genome_size = genome_size
		self.genome = int(genome_size/DELTA_X)
		self.TSS_pos = (tss['TSS_pos'].values/DELTA_X).astype(int)
		self.barr_fix = (prot['prot_pos'].values/DELTA_X).astype(int)
		self.strands = self.str2num(gff_df['strand'].values)

	def str2num(self, s):
    		s[s == '+'] = 1 #True
    		s[s == '-'] = -1 #False
    		return s
