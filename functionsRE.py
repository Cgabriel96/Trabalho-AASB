from re import search, finditer, sub, match
from MySeq import MySeq
from MyBlast import MyBlast

'''
Funções para procurar ERs em sequências
'''


def procuraPadraoER(seq, pad):
    mo = search(pad, seq)
    if mo is not None:
        # caso encontre o padrão na sequência, devolve a posição de início da primeira ocorrência
        return mo.span()[0]
    else:
        return -1


def procuraTodasOcsER(seq, pad):
    mos = finditer(pad,
                   seq)  # Return an iterator yielding MatchObject instances over all non-overlapping matches for the RE pattern in string.
    res = []
    for x in mos:
        res.append(x.span()[0])  # o primeiro indice representa qual a posiçao onde começa o padrao
    return res


def procuraTodasOcsERcomOverlap(seq, pad):
    # '?=pad' - Specifies position using a pattern. Doesn't have a range.
    return procuraTodasOcsER(seq, '(?=' + pad + ')')


def searchER():
    seq = input('Sequencia: ')
    pad = input('Padrao: ')

    res = procuraPadraoER(seq, pad)
    if res >= 0:
        print('Padrao descoberto na posicao: ', res)
    else:
        print('Padrao nao encontrado')

    todos_res = procuraTodasOcsER(seq, pad)
    if len(todos_res) > 0:
        print('Padrao descoberto nas posicoes: ', todos_res)
    else:
        print('Padrao nao encontrado')

    todos_ov = procuraTodasOcsERcomOverlap(seq, pad)
    if len(todos_ov) > 0:
        print('Padrao descoberto nas posicoes: ', todos_ov)
    else:
        print('Padrao nao encontrado')


def readFASTA(filename):
    fh = open(filename)
    name = fh.readline()[1:-1]
    sequence = ''
    for line in fh:
        sequence += sub('\n', '', line)     # sub - substititui white spaces (\s) por '' em cada linha...penso eu
    #print('Nome: ' + name)
    #print(name)
    #print('Sequência: ' + sequence)
    #print(sequence)
    fh.close()
    return sequence


def readTEXT(filename):
    file = open(filename)
    linhas = file.readlines()
    seq = ''
    if linhas[0][0] != '>':
        for li in linhas:
            seq += sub('\n', '', li)
    else:
        for li in range(1, len(linhas)):
            seq += sub('\n', '', linhas[li])
    file.close()
    #print(seq)
    return seq


def readfasta(filename):        # precisa de ser otimizado
    seq = filename
    fh = open(seq, "r")
    line = fh.readline()
    id_meta = []
    meta = ""
    sequence = ""
    while line:
        line = line.rstrip("\n")
        if ">" in line:
            meta = line
        else:
            sequence += sub('\s', '', line)
        line = fh.readline()
    print(meta)
    print(sequence)
    seqo = search(r'(^\w*\S*)(\s\D(.*))', meta)
    nome = seqo.group(2)
    id_meta.append(nome)
    id_ = seqo.group(1)
    id_meta.append(id_)
    #accnumber = seqo.group(0)
    #id_meta.append(accnumber)
    numerodenuc = len(sequence)
    id_meta.append(numerodenuc)
    print(nome)     # da mal
    print(id_)
    #print(accnumber)
    print(numerodenuc)
    fh.close()
    return id_meta


if __name__ == '__main__':
    #searchER()
    seq = readFASTA('sequence.fasta')
    myseq = MySeq(seq)
    print()
    prot = myseq.orfs()
    for prota in prot:
        print(prota)
    print()
    protaaa = myseq.maiorProteinaORFs()
    myseq.maiorProteinaORFs().printseq()

    mb = MyBlast('blast_test.fasta', 11)
    r = mb.bestAlignment(protaaa)
    print(r)
