
from MySeq import MySeq
from MatrixNum import MatrixNum

class MyMotifs:
    '''
    impolementa a estrutura de dados e metodos precisos para criar a PWM (position weighted matrix) para determinar
    diferentes representaçoes deterministicas dos motifs e para determinar a probabilidade de ocorrencia do motif ao
    longo da seq
    '''

    def __init__(self, seqs):
        self.seqs = seqs        # objetos classe MySeq
        self.tam = len(seqs[0])
        self.alfabeto = seqs[0].alfabeto()
        self.calculaContagens()
        self.criaPWM()

    def calculaContagens(self):
        '''
        cria matriz de contagens a partir das instancias (seqs)
        '''
        self.counts = createMatZeros(len(self.alfabeto), self.tam)      # sera que e o MatrixNum? probably | createMatZeros()
        for s in self.seqs:
            for i in range(self.tam):
                lin = self.alfabeto.index(s[i])
                self.counts[lin][i] += 1

    def criaPWM(self):
        '''
        cria uma matriz pwm a partir das contagens
        '''
        if self.counts is None:
            self.calculaContagens()
        self.pwm = createMatZeros(len(self.alfabeto), self.tam)      # | createMatZeros()
        for i in range(len(self.alfabeto)):
            for j in range(self.tam):
                self.pwm[i][j] = float(self.counts[i][j]) / len(self.seqs)

    def consenso(self):
        '''
        define o consenso que e dado pelos caracteres mais comuns em cada posiçao (coluna) do motif
        '''
        res = ''
        for j in range(self.tam):
            maxcol = self.counts[0][j]
            maxcoli = 0
            for i in range(1, len(self.alfabeto)):
                if self.counts[i][j] > maxcol:
                    maxcol = self.counts[i][j]
                    maxcoli = i
            res += self.alfabeto[maxcoli]
        return res

    def probabSeq(self, seq):
        '''
        da a probabilidade de uma dada seq inst do mesmo tamanho das instancias ser gerada pelo motif
        '''
        res = 1.0
        for i in range(self.tam):
            lin = self.alfabeto.index(seq[i])
            res *= self.pwm[lin][i]
        return res

    def mostProbableSeq(self, seq):
        '''
        da a posicao mais provavel de ocorrencia do motif numa seq de tamanho maior do que o motif
        '''
        maximo = -1.0
        maxind = -1
        for k in range(len(seq)-self.tam):
            p = self.probabSeq(seq[k:k+self.tam])
            if p > maximo:
                maximo = p
                maxind = k
        return maxind

    def score(self):
        score = 0
        for i in range(len(self.pwm[0])):
            maxi = -1
            for j in range(len(self.pwm)):
                if self.pwm[j][i] > maxi:
                    maxi = self.pwm[j][i]
            score += maxi
        return score


def createMatZeros(rows, cols):         # nao devia de ser da class MatrixNum?
    mat = []                            # da erro nao sei pq
    for i in range(rows):
        mat.append([])
        for j in range(cols):
            mat[i].append(0.0)
    return mat


def printmat(mat):                      # metodo da MatrixNum
    for r in mat:
        print(r)
    print()

# --------------------------\\-----------------------------

# TESTS


def test():
    seq1 = MySeq('AAAGTT')
    seq2 = MySeq('CACGTG')
    seq3 = MySeq('TTGGGT')
    seq4 = MySeq('GACCGT')
    seq5 = MySeq('AACCAT')
    seq6 = MySeq('AACCCT')
    seq7 = MySeq('AAACCT')
    seq8 = MySeq('GAACCT')
    lseqs = [seq1, seq2, seq3, seq4, seq5, seq6, seq7, seq8]
    motifs = MyMotifs(lseqs)
    printmat(motifs.counts)
    printmat(motifs.pwm)
    print(motifs.consenso())
    print(motifs.alfabeto)
    print()
    #print(motifs.probabSeq('AAACCT'))
    #print(motifs.probabSeq('ATACAG'))
    print()
    #print(motifs.mostProbableSeq('CTATAAACCTTACATC'))
    print(motifs.score())


if __name__ == '__main__':
    test()