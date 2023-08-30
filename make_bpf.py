import sys
import argparse

### Make BPF ###

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Ingrese posiciones ej(-p 1,4,5)')
    #parser.add_argument('pdb', action='store')
    parser.add_argument('-p', '--pos', action='store',type=str)
    args = parser.parse_args()

    #In = open(args.pdb ,'r').readlines()
    In = open('interaction_sites.pdb' ,'r').readlines()
    Out = open('ligand_bias.bpf' , 'w')

    try:
        Entrada_Sistema = (args.pos)
        Corte = Entrada_Sistema.split(',')
        Interes = []
        for j in (Corte):
            Interes.append(int(j))
    except:
        print('Ingrese posiciones ej(-p 1,4,5)')
        exit()


    Vset = -2.00
    r = 1.20

    Out.write('x    y    z     Vset    r     type\n')

    for lines in In:
        if int(lines[23:27]) in Interes: 
            Out.write('{}     {}     {}     {}\n'.format(lines[30:55].strip(),Vset,r,(lines[17:20]).lower()))
    
