
                                   _______    _______    ________  _________  _________
                                  /       \  /       \  /       / /__   ___/ /        /
                                 /   _____/ /     ___/ /       /    /  /    /        /
                                /  /       /  /\  \   /       /    /  /    /  ___   /
                               /__/       /__/ /__/  /_______/    /__/    /__/  /__/


Este ficheiro contém informações acerca das funções utilizadas na Shell do programa: PROTA

O objetivo deste programa é a manipulação/ armazenamento de sequências de aminoácidos em uma
base de dados.

---
NOTA!!
QUALQUER FICHEIRO A SER UTILIZADO POR ESTE PROGRAMA DEVE ENCONTRAR-SE NA DIRECTORIA DO MESMO.
---

O programa inicia criando um dicionário em pickle caso este não existe, caso este exista
ele carrega-o para a base de dados com toda a informação depositado previamente.
Posteriormente é apresentado um menu contendo todos os comandos possíveis do programa,
bem como uma breve explicação do comando. Optou-se por este método em vez do comando help
do cmd, por este não apresentar imediatamente a informação de cada comando, facilitando
visualmente ao utilizador o uso do programa.

--- NOVASEQ ---
A primeira função do Script é designada por novaseq e permite introduzir uma sequência
de aminoácidos manualmente contendo 2 argumentos (1- o ID pretendido para a sequência,
2- a sequência). Após a inserção desta sequência na base de dados são calculadas algumas
features como a maior proteína (caso seja uma sequência com várias proteínas), o seu tamanho
e a frequência de polaridade/ carga.

--- NOVAFEATURE ---
Dado um argumento com ID da sequência, nome da Feature a adicionar e o seu respectivo valor.
É acrescentado à base de dados.

--- EDITARFEATURE ---
Dado um argumento com ID da sequência, nome da Feature a editar, e o novo valor para a Feature.
Actualiza a Feature na base de dados com o novo valor.

--- PROCURAID ---
Dado um argumento com o ID da sequência, é apresentada toda a informação respectiva a essa
sequência.

--- FREQSIMBOL ---
Dado um argumento com o padrao a procurar, e o ID da sequência ou "BD" caso pretenda todas
as sequências da base de dados local.
Retorna a frequência desse padrão na determinada sequência/ ou BD.

--- EXPORT ---
Dado um argumento com o nome do ficheiro em que ficará guardada a informação, e o ID da
sequência ou BD (para a base de dados local).
Cria um ficheiro externo em formato .txt com o ID da sequência e na linha seguinte a sequência.
Caso seja a base de dados, cria espaços entre cada sequência.

--- IMPORTSEQ ---
Dado um argumento com um nome de ficheiro externo já existente e o seu formato (.fasta ou .txt)
Importa para a base de dados local a sequência presente nesse ficheiro e atribui como ID o nome
do ficheiro (sem a extensão).

--- BLASTLOCAL ---
Função sem argumentos, inicia por perguntar ao utilizador se pretende introduzir uma sequência
manualmente ou importar um ficheiro contendo a sequência para a qual pretende realizar o blast
local.
Retorna o ID da sequência da base de dados local, mais semelhante à sequência dada.

--- ARVORE ---
Função sem argumentos que itera todas as sequências da base de dados local para uma lista.
Os parametros de uma Matriz de substituição com os valores: Match = 3, Mismatch = -1
Para o alinhamento o gap = -1
O algoritmo UPGMA cria uma matriz de distâncias criada.
Retorna uma árvore binária apartir da matriz de distâncias.

--- ALINMULTI ---
Função sem argumentos que itera todas as sequências da base de dados local para uma lista.
Utiliza a matriz de substituição BLOSUM62 e define a penalização para gaps = -8.
Retorna um alinhamento múltiplo entre todas as sequências da base de dados local.

--- TAMANHODIC ---
Função sem argumentos, retorna o número de sequências presente na base de dados local.

--- DICIONARIO ---
Função sem argumentos, retorna o dicionário contendo toda a informação da base de dados.

--- GRAVAR ---
Função sem argumentos, grava o estado da base de dados no presente momento.
O utilizador pode encerrar a plataforma de python sem necessitar sair da shell para gravar.

--- SAIR ---
Função sem argumentos, sai da shell do programa e grava automaticamente o estado da base de
dados até ao momento.