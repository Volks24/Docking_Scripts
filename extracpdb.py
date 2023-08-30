import sys
import  os


if __name__ == '__main__':
    
    
    if (len(sys.argv)) <= 1:
        print('falta archivos dlg')
        sys. exit()
    else:
        file_name = sys.argv[1]
        if 'dlg' not in file_name:
            print('no es un archivo *.dlg')
            sys. exit()


    #### Obtengo PDB y RMSD ####

    input_file = open('{}'.format(file_name),'r').readlines()
    rank_value = 1
    ranks = {}

    for lines in input_file:
        if lines[0:6] == 'DOCKED':
            mark = (lines.split(' '))[1]
            if mark == 'USER':
                if 'Run' in lines:
                    rank = int(lines[22:25])
                    pdbqt_out = open('Run{}.pdbqt'.format(rank),'w')
            if mark in ['REMARK' , 'HETATM' , 'ATOM' ,'ROOT', 'ENDROOT' , 'BRANCH' ,'ENDBRANCH' ,'TORSDOF' , 'TER']:
                pdbqt_out.write(lines[8:])
            if mark in ['ENDMDL']:
                pdbqt_out.close()
        if 'RANKING' in lines:
            ranks[int(lines[13:20])] = rank_value
            rank_value = rank_value + 1
    pdbqt_out.close()

    # Renombrar el archivo
    for keys in ranks:  
        os.rename('Run{}.pdbqt'.format(keys), 'Rank{}.pdbqt'.format(ranks[keys]))

    for j in range(1,rank_value):
        print(j)
        os.system('obabel -ipdbqt Rank{}.pdbqt -opdb -O rank{}.pdb'.format(j,j))

    file_out = open('{}.pdb'.format(file_name.split('.')[0]) ,'w')    

    for j in range(1,rank_value):
        file_in = open('rank{}.pdb'.format(j) , 'r').readlines()
        for lines in file_in:
            file_out.write(lines)
        
    file_out.close()