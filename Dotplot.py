import sys


class Dotplot:
    '''
    representa matriz de pontos com atributos das seq a colocar nas linhas e colunas (seq1, seq2)
    e a matriz de pontos (lista de listas ou objeto MatrixNum) - mat
    '''

    def __init__(self, seq1, seq2):
        self.seq1 = seq1
        self.seq2 = seq2
        self.mat = []                         # é uma matriz onde vao as seq1 e seq2
        for i in range(len(seq1)):
            self.mat.append([])               # cria uma lista vazia
            for j in range(len(seq2)):
                self.mat[i].append(0)         # faz a matriz mas poe tudo a zeros
        # ajuda para ver visualmente
        # print(self.mat)

    def buildMatSimple(self):
        '''
        constroi a matriz de pontos
        '''
        for i in range(len(self.seq1)):
            for j in range(len(self.seq2)):
                if self.seq1[i] == self.seq2[j]:
                    self.mat[i][j] = 1
                else:
                    self.mat[i][j] = 0

    def printmat(self):
        '''
        imprime a matriz
        '''
        for letter in self.seq2:
            sys.stdout.write('  '+letter)
        sys.stdout.write('\n')
        for i in range(len(self.mat)):
            sys.stdout.write(self.seq1[i])
            for j in range(len(self.mat[i])):
                if self.mat[i][j] >= 1:
                    sys.stdout.write(' * ')
                else:
                    sys.stdout.write('   ')
            sys.stdout.write('\n')

    def createMat(self, seq1, seq2):            # cria uma matriz de zeros
        # tem que ser quadrada a matriz??
        self.seq1 = seq1
        self.seq2 = seq2
        self.mat = []
        for i in range(len(seq1)):
            self.mat.append([])
            for j in range(len(seq2)):
                self.mat[i].append(0)
        return self.mat

    def extendedDotPlot(self, seq1, seq2, window, stringency):
        mat = self.createMat(seq1, seq2)
        start = int(window/2)
        for i in range(start, len(seq1)-start):
            for j in range(start, len(seq2)-start):
                matches = 0
                l = j - start
                for k in range(i-start, i+start+1):
                    if seq1[k] == seq2[l]:
                        matches += 1
                    l += 1
                if matches >= stringency:
                    mat[i][j] = 1
                else:
                    mat[i][j] = 0
        return mat


# --------------------------\\\-----------------------------

# TESTS

def test():
    x = Dotplot('ACTTGA', 'AGGCTA')
    x.buildMatSimple()
    x.printmat()
    print('\n')
    x.extendedDotPlot('ACTCCCGAATGTGA', 'ACTCCCGAATGTGA', 3, 2)  # seq vs self
    x.printmat()


if __name__ == '__main__':
    test()



