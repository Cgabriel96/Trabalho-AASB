
from re import search

'''
x: representa qualquer aminoácido
[ ]: representa lista de aminoácidos alternativos
x(2): representa dois aa’s quaisquer
os “-” são usados para separar as várias posições
'''


def findZyncFinger(seq):
    '''
    Função para procurar este domínio numa proteína com sequência seq
    '''
    regexp = "C.H.[LIVMFY]C.{2}C[LIVMYA]"       # Procurar este domínio numa proteína
    mo = search(regexp, seq)
    if mo is not None:
        return mo.span()[0]
    else:
        return -1


def findProsite(seq, profile):
    '''
    função que trata os diversos símbolos para converter domínios da Prosite em ERs
    '''
    regexp = profile.replace('-', '')
    regexp = regexp.replace('x', '.')
    regexp = regexp.replace('(', '{')
    regexp = regexp.replace(')', '}')
    mo = search(regexp, seq)
    if mo is not None:
        return mo.span()[0]
    else:
        return -1

    
def test():
    seq = "HKMMLASCKHLLCLKCIVKLG"
    print(findZyncFinger(seq))
    print()
    print(findProsite(seq, "C-x-H-x-[LIVMFY]-C-x(2)-C-[LIVMYA]"))


if __name__ == '__main__':
    test()
