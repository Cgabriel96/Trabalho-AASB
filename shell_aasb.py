
from cmd import *
import pickle
import sys
from MySeq import MySeq
from MatrixNum import MatrixNum
from SubstMatrix import SubstMatrix
from Dotplot import Dotplot
from AlignSeq import AlignSeq
from MyAlign import MyAlign
from MultipleAlign import MultipleAlign
from MyBlast import MyBlast
from MyBlast import Blastp
from BinaryTrees import BinaryTree
from ClustHier import ClustHier
from upgma import UPGMA
from MyMotifs import MyMotifs
import functionsRE
import iubToRE
import Prosite
import re

try:
    pickle_in = open('dict.pickle', 'rb')
    dicionario = pickle.load(pickle_in)
except:
    dicionario = {}
    pickle_out = open('dict.pickle', 'wb')
    pickle.dump(dicionario, pickle_out)


class ProtaShell(Cmd):
    intro = '''
                                   _______    _______    ________  _________  _________
                                  /       \  /       \  /       / /__   ___/ /        /
                                 /   _____/ /     ___/ /       /    /  /    /        /
                                /  /       /  /\  \   /       /    /  /    /  ___   /
                               /__/       /__/ /__/  /_______/    /__/    /__/  /__/
                                  
++====================================================================================================================++
||                                                                                                                    ||            
||    novaseq -> introduzir uma sequência de proteína manualmente                                                     ||
||    remover_seq -> remove a sequência especificada na base de dados local                                           ||
||    importar_seq -> importar uma sequência de uma ficheiro (.fasta ou .txt)                                         ||
||    blast_db -> blast da sequência proteica contra a base de dados local                                            ||
||    blast_ex -> blast da sequência proteica contra a base de dados externa (ex: swissprot, nr)                      ||
||    freq_simbolo <padrão> -> efetua a contagem de certos padrões na seq de aa                                       ||
||    procura <file_input file_output database> -> efetua a procura de padrões na seq de aa usando ER                 ||
||    arvore <file_input> -> cria uma árvore filogenética a partir do AM contra bd local, usando o algoritmo UPGMA    ||
||    alinhamento_multiplo <file_input> -> efetua o alinhamento múltiplo de sequências homólogas                      ||
||    sair -> fecha o programa                                                                                        ||
||                                                                                                                    ||
||                          Pressione < ? comando > para mais informações acerca do comando                           ||
||                 (qualquer ficheiro importado/exportado DEVE estar na mesma working directory!)                     ||
++====================================================================================================================++ 
    '''
    prompt = '@Prota > '

    #del dicionario['seq2']

    def do_novaseq(self, arg):
        '- Introduz uma nova sequência na base de dados local: novaseq <id_seq> <seq de aa>'
        try:
            lista_arg = arg.split()
            num_args = len(lista_arg)
            if num_args == 2:
                idseq = lista_arg[0]
                print(idseq)
                seq = lista_arg[1]
                print(seq)
                alin = MySeq(seq, tipo='protein')
                if alin.validaER():
                    print(f'Sequência com id:{idseq} foi adicionada à base de dados.')
                    dicionario[idseq] = {'seq': seq}
                    dicti = {'tamanho': alin.__len__(), 'maior proteína': alin.maiorProteina().seq, 'polaridade e carga': alin.validaprot()}
                    dicionario[idseq].update(dicti)
                else:
                    print("\nA sequência não é de aminoácidos!")
            else:
                print("\nNúmero de argumentos inválido!")
            print(dicionario)
        except:
            print("\nErro ao adicionar a sequência!")

    def do_novafeature(self, arg):
        '- Introduz uma nova feature da sequência na base de dados local: novafeature <id_seq> <nome_feature> <feature>'
        try:
            lista_arg = arg.split()
            num_args = len(lista_arg)
            if num_args == 3:
                idseq = lista_arg[0]
                featname = lista_arg[1]
                feature = lista_arg[2]
                dicti = {featname: feature}
                dicionario[idseq].update(dicti)
                print(dicionario)
            else:
                print("\nNúmero de argumentos inválido!")
        except:
            print("\nErro ao adicionar a feature!")

    def do_editarfeature(self, arg):
        '- Edita o conteúdo existente da feature selecionada: editarfeature <id_seq> <nome_feature> <feature>'
        try:
            lista_arg = arg.split()
            num_arg = len(lista_arg)
            if num_arg == 3:
                idseq = lista_arg[0]
                feature = lista_arg[1]
                new_value = lista_arg[2]
                if dicionario[idseq][feature] is not None:
                    dicionario[idseq][feature] = new_value
                    print("\nFeature editada!!\n", dicionario[idseq])
                else:
                    print("ID ou Feature não existente!")
            else:
                print("\nNúmero de argumentos inválido!! Ex. seq1 funcao transportador transmembranar ")
        except:
            print("\nErro ao editar a feature!")

    def do_procuraid(self, arg):
        '- Procura uma sequência na base de dados local através do id: procuraid <id_seq>'
        try:
            lista_arg = arg.split()
            num_arg = len(lista_arg)
            if num_arg == 1:
                idseq = lista_arg[0]
                if idseq in dicionario.keys():
                    print(dicionario[idseq])
                else:
                    print("\nID da sequência não encontrado")
            else:
                print("\nNúmero de argumentos inválido!! Ex. seq1 ")
        except:
            print("\nErro a procurar o ID")

    def do_freqsimbol(self, arg):
        '- Procura padrões normais e/ou expressões regulares nas sequências, dá return do número de ocorrências: freqsimbol <pattern_ou_ER>'
        try:
            lista_arg = arg.split()
            num_arg = len(lista_arg)
            if num_arg == 2:
                padrao = lista_arg[0]
                dados = lista_arg[1]
                if dados in dicionario.keys():
                    x = re.findall(padrao.upper(), dicionario[dados]["seq"])
                    print(len(x))
                else:
                    if dados.upper() == "DB":
                        cont = 0
                        for key in dicionario.keys():
                            y = re.findall(padrao.upper(), dicionario[key]["seq"])
                            cont += len(y)
                        print("A frequência do padrão na base de dados local é de: ", cont)
                    else:
                        print("Base de dados não encontrada")
            else:
                print("\nNúmero de argumentos inválido!!")
        except:
            print("\nErro ao verificar a frequência")

    def do_export(self, arg):
        '- Exporta uma sequência, dado o id, ou exporta as sequências todas da base de dados para um ficheiro .txt: export <file_name> <seq_ou_bd>'
        try:
            lista_arg = arg.split()
            num_arg = len(lista_arg)
            if num_arg == 2:
                file = lista_arg[0] + '.txt'
                opc = lista_arg[1]
                f = open(file, 'w+')
                if opc in dicionario.keys():
                    seq = dicionario[opc]['seq']
                    print('>' + opc + '\n' + seq, file=f)
                else:
                    if opc.upper() == "DB":
                        for key in dicionario.keys():
                            seq = dicionario[key]['seq']
                            print('>' + key + '\n' + seq, file=f)
                    else:
                        print("Base de dados não encontrada")
                f.close()
            else:
                print("\nNúmero de argumentos inválido!!")
        except:
            print("\nErro ao exportar a(s) frequência(s)!")

    def do_importseq(self, arg):
        '- Importa uma sequência a partir de um ficheiro dado, usando o nome sem extensão como seq_id: importseq <file_id>'
        try:
            lista_arg = arg.split()
            num_arg = len(lista_arg)
            if num_arg == 2:
                filename = lista_arg[0]
                form = lista_arg[1]
                if form.upper() == 'TXT':
                    filename2 = filename + '.txt'
                    s = functionsRE.readTEXT(filename2)
                    seq = MySeq(s, tipo='protein')
                    if seq.validaER():
                        dicionario[filename] = {'seq': seq.seq}
                        dicti = {'tamanho': seq.__len__(), 'maior proteína': seq.maiorProteina().seq, 'polaridade e carga': seq.validaprot()}
                        dicionario[filename].update(dicti)
                elif form.upper() == 'FASTA':
                    filename2 = filename + '.fasta'
                    s = functionsRE.readFASTA(filename2)
                    seq = MySeq(s, tipo='protein')
                    if seq.validaER():
                        dicionario[filename] = {'seq': seq.seq}
                        dicti = {'tamanho': seq.__len__(), 'maior proteína': seq.maiorProteina().seq,  'polaridade e carga': seq.validaprot()}
                        dicionario[filename].update(dicti)
                else:
                    print('Formato inválido, selecione FASTA ou TXT')
        except:
            print('Erro ao importar a sequência')

    def do_blastlocal(self, arg):
        '- Realiza um Bast de uma sequência importada, através de um ficheiro ou do input do utilizador, contra todas as seqs na base de dados, fazendo return do id da seq mais semelhante: blastlocal'
        try:
            option = int(input('''
    Selecione o número correspondente ao gene a analisar: 
    1 -> sequência 
    2 -> de um ficheiro
            
    Opção: '''))
            blast = MyBlast()
            for ke in dicionario.keys():
                seq = dicionario[ke]['seq']
                blast.addSequenceDB(seq)
            if option == 1:
                query = input('\nSequência:\n')
                r = blast.bestAlignment(query)
                seqdb = blast.db[r[4]]
                for kei in dicionario.keys():
                    if dicionario[kei]['seq'] == seqdb:
                        print('\nSequência mais parecida na base de dados da sequência introduzida: '+kei+'\n')
            if option == 2:
                file = input('\nFicheiro: ')
                query = functionsRE.readTEXT(file)
                r = blast.bestAlignment(query)
                seqdb = blast.db[r[4]]
                for kei in dicionario.keys():
                    if dicionario[kei]['seq'] == seqdb:
                        print('\nSequência mais parecida na base de dados da sequência introduzida: '+kei+'\n')
        except:
            print('Erro ao realizar o blast local')

    def do_arvore(self, arg):
        '''- Cria uma árvore de todas as seqs na base de dados, usando o algoritmo UPGMA: arvore

        -- é necessário pelo menos trêss seqs para formar uma àrvore, e estas não podem ter underscore --'''
        lista_seqs = []
        for id in dicionario.keys():
            #print(dicionario[id]['seq'])
            id_seq = MySeq(dicionario[id]['seq'], tipo='protein')
            lista_seqs.append(id_seq)
            #print(id)
            #print(lista_seqs)
        sm = SubstMatrix()
        sm.createFromMatchPars(0, -1, lista_seqs[0].alfabeto())
        alseq = AlignSeq(sm, -1)
        up = UPGMA(lista_seqs, alseq)
        #up.matdist.printmat()
        arv = up.run()
        arv.printtree()
        print()

    def do_alinmulti(self, arg):
        '- Realiza um alinhamento múltiplo com todas as seqs da base de dados: alimulti'
        try:
            lista_seqs = []
            for id in dicionario.keys():
                #print(dicionario[id]['seq'])
                id_seq = MySeq(dicionario[id]['seq'], tipo='protein')
                lista_seqs.append(id_seq)
            sm = SubstMatrix()
            sm.loadFromFile('blosum62.mat', '\t')
            alseq = AlignSeq(sm, -8)
            ma = MultipleAlign(lista_seqs, alseq)
            alinm = ma.alignConsensus()
            print(alinm)
            #print(ma.scoreSP(alinm))
        except:
            print('Erro ao realizar alinhamento múltiplo\n')

    def do_prosite(self, arg):
        pass

    def do_tamanhodic(self, arg):
        '- Devolve o número de seqs na base de dados: tamanhodic'
        print(str(len(dicionario))+' sequência(s)')

    def do_dicionario(self, arg):
        '- Devolve o id, a sequência e as features de todas as seqs na base de dados: dicionario'
        for ki in dicionario:
            print(ki)
            for j in dicionario[ki]:
                if j == 'features':
                    for i in dicionario[ki][j]:
                        print(j)
                        print(i + ': ' + dicionario[ki][j][i])
                else:
                    print(j + ': ' + str(dicionario[ki][j]))
            print()

    def do_gravar(self, arg):
        '- Grava o progresso da base de dados: gravar'
        pickle_out = open('dict.pickle', 'wb')
        pickle.dump(dicionario, pickle_out)
        pickle_out.close()
        print(dicionario)

    def do_sair(self, arg):
        '- Sair do programa, gravando tudo na base de dados: sair'
        pickle_out = open('dict.pickle', 'wb')
        pickle.dump(dicionario, pickle_out)
        pickle_out.close()
        #print(dicionario)
        print('Desligar shell...')
        sys.exit()



if __name__ == '__main__':
    sh = ProtaShell()
    sh.cmdloop()