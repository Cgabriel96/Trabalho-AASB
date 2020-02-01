
from re import finditer

def iubToRE(iub):
	'''
	converta a representação de uma enzima de restrição na expressão regular correspondente
	'''
	dic = {
			'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T',
			'R': '[GA]', 'Y': '[CT]', 'M': '[AC]', 'K': '[GT]', 'S': '[GC]', 'W': '[AT]',
			'B': '[CGT]', 'D': '[AGT]', 'H': '[ACT]', 'V': '[ACG]', 'N': '[ACGT]'
	}
	site = iub.replace("^", "")
	regexp = ""
	for c in site:
		regexp += dic[c]
	#print(regexp)
	return regexp


def cutPositions (enzyme, sequence):
	'''
	determina as posições onde uma enzima de restrição “corta” a sequência; retorna uma lista de índices
	'''
	cutpos = enzyme.find('^')    # retorna primeira ocorrência, (onde há corte dentro do padrão)
	regexp = iubToRE(enzyme)
	matches = finditer(regexp, sequence)
	locs = [] 			# lista de índices onde a enzima de restrição corta a sequência
	for m in matches:
		locs.append(m.start() + cutpos)
	return locs


def cutSubsequences(locs, sequence):    # sub-sequências resultantes do corte da sequência
	'''
	determina as subsequências resultantes do corte da sequência seq, usando os índices na lista locs
	'''
	res = []
	positions = locs
	positions.insert(0, 0) 		# insere 0 na posição 0 da lista
	positions.append(len(sequence)) 		# insere len da seq no fim da lista
	for i in range(len(positions)-1):
		res.append(sequence[positions[i]:positions[i+1]])
	# só dá uma subsequência?
	return res


if __name__ == '__main__':
	iubToRE("G^AMTV")
	pos = cutPositions('G^ATTC', 'GTAGAAGATTCTGAGATCGATTC')
	cut = cutSubsequences(pos, 'GTAGAAGATTCTGAGATCGATTC')
	print(cut)