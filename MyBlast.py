from Bio.Blast import NCBIWWW
from Bio import SeqIO

class MyBlast:
    '''
    simplified blast...
        |_ não usa matriz de scoring: considera apenas matches perfeitos entre a query e as sequências da BD
        |_ Critérios simples usados para a extensão dos hits
        |_ Score será a contagem do nº de matches

    procura sequências mais similares em bases de dados
    BLAST - Basic Local Alignment Search Tool
    Procura palavras do mesmo tamanho nas sequências da BD, que comparadas com as palavras (e.g. 3 AAs ou 5-15 bases de
    DNA) da sequência query tenham alta similaridade
    '''

    def __init__(self, filename=None, w=3):
        '''
        Construtor a partir do nome do ficheiro da base de dados e do valor de w
        '''
        if filename is not None:
            self.readDatabase(filename)
        else:
            self.db = []
        self.w = w
        self.map = None

    def readDatabase(self, filename):
        '''
        assumindo que uma base de dados é apenas um conjunto de seqs num ficheiro de texto, uma por linha
        '''
        self.db = []
        f = open(filename)
        for line in f:
            self.db.append(line.rstrip())
        f.close()
        
    def addSequenceDB(self, seq):
        if seq not in self.db:              # para eliminar a redundancia da bd???
            self.db.append(seq)
        
    def buildMap (self, query):
        '''
        criar um mapa do query
        dada a query(seq) e o w(tamanho de palavras), vai-se mapear cada ocorrencia de subseqs de tamanho w
            |_ o resultado sera um dicionario (key - subseqs; value - lista de indices(pos onde ocorre a subseq)
        '''
        self.map = {}  
        for i in range(len(query)-self.w+1):        # todas as palavras de tamanho w da query
            subseq = query[i:i+self.w]          # palavra de tamanho w
            if subseq in self.map:
                self.map[subseq].append(i)
            else:
                self.map[subseq] = [i]

    def getHits (self, seq, query):
        '''
        descobrir todas as co-ocorrências de palavras de tamanho W entre a query e a sequência
        (so sao considerados matches perfeitos de tamanho w)
        o resultado é uma lista de hits (cada hit sendo um tuplo)
            |_ 1o elemento é a pos da palavra na query; 2o elemento é a pos da palavra na seq
        '''
        if self.map is None: 
            self.buildMap(query)
        res = []
        for i in range(len(seq)-self.w+1):
            subseq = seq[i:i+self.w]
            if subseq in self.map:
                l = self.map[subseq]
                for ind in l:
                    res.append((ind, i))
        return res
        
    def extendsHit (self, seq, hit, query):
        '''
        cada hit sera extendido em ambas as direçoes
        um hit sera extendido enquanto o no de matches nessa direçao for pelo menos metade do tamanho da extensao
        seq - sequencia da bd
        o resultado sera um tuplo com 4 valores:
            |_ inicio do hit extendido na query
            |_ inicio do hit extendido na seq
            |_ tamanho do hit extendido
            |_ score do hit extendido (= no de matches)
        '''
        stq, sts = hit[0], hit[1]           # stq - inicio do hit na query; sts - inicio do hit na seq
        matfw = 0       # no de matches a frente
        k = 0           # pos andadas em frente
        bestk = 0
        while 2*matfw >= k and stq+self.w+k < len(query) and sts+self.w+k < len(seq):
            if query[stq+self.w+k] == seq[sts+self.w+k]: 
                matfw += 1
                bestk = k+1
            k += 1
        size = self.w + bestk       # tamanho do hit extendido
    
        k = 0           # pos andadas para tras
        matbw = 0       # no de matches para tras
        bestk = 0
        while 2*matbw >= k and stq > k and sts > k:
            if query[stq-k-1] == seq[sts-k-1]: 
                matbw += 1
                bestk = k+1
            k += 1
        size += bestk
        # o stq-bestk e sts-bestk representam o indice onde começa a seq de hit extendida
        return stq-bestk, sts-bestk, size, self.w+matfw+matbw
        
    def hitBestScore(self, seq, query):
        '''
        recupera todos os hits e extende cada um deles
        score e calculado como o no de matches
        seleciona o maior score; caso empate - seleciona com menor tamanho
        '''
        hits = self.getHits(seq, query)
        bestScore = -1.0
        best = ()
        for h in hits:
            ext = self.extendsHit(seq, h, query)
            score = ext[3]              # o score do hit extendido (numero de matches)
            if score > bestScore or (score == bestScore and ext[2] < best[2]):
                bestScore = score
                best = ext
        return best

    def bestAlignment (self, query):
        '''
        calcula o melhor alinhamento da query com uma das seqs da BD
        '''
        self.buildMap(query)
        bestScore = -1.0
        res = (0, 0, 0, 0, 0)
        for k in range(0, len(self.db)):            # o k é o indice da seq na lista de bd
            bestSeq = self.hitBestScore(self.db[k], query)
            if bestSeq != ():
                score = bestSeq[3]  
                if score > bestScore or (score == bestScore and bestSeq[2] < res[2]):
                    bestScore = score
                    res = bestSeq[0], bestSeq[1], bestSeq[2], bestSeq[3], k
        if bestScore < 0:
            return ()
        else:
            return res


class Blastp:

    def __init__(self, input_file, output_file, database="nr"):
        self.__input_file = input_file
        self.__output_file = output_file
        self.__database = database

    def blastp(self):
        record = SeqIO.read(self.__input_file, "fasta")
        print("Blast a decorrer. Aguarde.")
        ficheiro = open(self.__output_file, "w")
        blast_result = NCBIWWW.qblast("blastp", self.__database, record.format("fasta"))
        ficheiro.write(blast_result.read() + "\n")
        print("O resultado do Blast exportado para o ficheiro %s." % (self.__output_file))
        ficheiro.close()
        blast_result.close()


# --------------------------\\-----------------------------

# TESTS

def test1():
    mb = MyBlast('seqBlast.txt', 11)
    query = 'gacgcctcgcgctcgcgcgctgaggcaaaaaaaaaaaaaaaaaaaatcggatagctagctgagcgctcgatagcgcgttcgctgcatcgcgtatagcgctgaagctcccggcgagctgtctgtaaatcggatctcatctcgctctatcct '
    r = mb.bestAlignment(query)
    print(r)

# (1, 38, 149, 108, 3)


def test2():
    mb = MyBlast('seqBlast.txt', 11)
    query2 = 'cgacgacgacgacgaatgatg'
    r = mb.bestAlignment(query2)
    print(r)       

# (0, 0, 21, 21, 4)


if __name__ == '__main__':
    test1()
    print()
    test2()
