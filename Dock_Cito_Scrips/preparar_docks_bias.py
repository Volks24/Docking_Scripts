import os
import pandas as pd
import shutil
from os import path


def Generar_Mapas_Bias(ruta_carpeta):
    
    AC = 0

    archivos_en_directorio = os.listdir(ruta_carpeta)
    file_name = 'receptor.maps.fld'
    archivos_biased_map = [archivo for archivo in archivos_en_directorio if archivo.endswith('.biased.map')]

    if 'receptor.AC.biased.map' in archivos_biased_map:
        AC = 1

    archivos_sin_extension = [archivo.replace('.biased.map', '.map') for archivo in archivos_en_directorio if archivo.endswith('.biased.map')]

    File = open(ruta_carpeta+'/'+file_name , 'r').readlines()
    num_label = 0
    for i, linea in enumerate(File):
        if "label=" in linea:
            num_label = num_label + 1
            indice_ultima_label = i
        elif linea[16:].split(' ')[0] in archivos_sin_extension:
            File[i] = File[i].replace('.map', '.biased.map')


    if AC == 1:
        Line_AC = 'label=AC-affinity	# component label for variable x'
        File.insert(indice_ultima_label + 1, Line_AC)    
        File.append('variable {} file=receptor.AC.biased.map filetype=ascii skip=6'.format(num_label+1))


    with open(ruta_carpeta+'/receptor.biased_maps.fld', 'w') as archivo_out:
        archivo_out.writelines(File)


if __name__ == '__main__':

    ### Armo DF con los ligandos ###
    os.chdir('Ligandos')
    ligando_df = pd.DataFrame(columns=['Ligand', 'Atom number'])
    lista_ligandos = os.listdir()
    for lig in lista_ligandos:
        number_of_atom = 0
        lig = lig.split('.')[0]
        ligando = open('{}.pdb'.format(lig,'r')).readlines()
        for lines in ligando:
            if ('HETATM' in lines) or ('ATOM'in lines):
                number_of_atom = number_of_atom + 1
        ligando_df.loc[len(ligando_df.index)] = [lig,number_of_atom]
        print('# atomos de {} : {}'.format(lig,len(ligando)))

    print(ligando_df)
    os.chdir('..')


    #### Preparo Archivos a Usar ###
    path_complete = os.getcwd()
    Lista_Ligandos = os.listdir(path_complete+'/Ligandos/')
    Lista_Proteinas = os.listdir(path_complete+'/Receptores/')
    grillas = pd.read_csv('grillas.csv' , sep=';' , names=['Receptor' , 'Coor' ,'Size'])

    grilla_chica = '50,50,50' 
    grilla_grande = '60,60,60'

    for Prot in Lista_Proteinas:
        Prot = Prot.split('.')[0]
        try:
            os.mkdir('{}'.format(Prot))
        except FileExistsError:
            pass
        for lig in Lista_Ligandos:
            lig = lig.split('.')[0]
            try:
                os.mkdir('{}/{}_{}'.format(Prot,Prot,lig))
            except FileExistsError:
                pass

            Size = ligando_df.query('Ligand == @lig')
            sub_set = grillas.query('Receptor == @Prot')
            Corte = (list(Size['Atom number']))
            print(Corte[0])
            npts = sub_set['Size'].tolist()
            Grid = sub_set['Coor'].tolist() 

            ### Caso ###
            print('Receptor {} Ligando {} \n'.format(Prot,lig))
        
            ### Muevo los archivos base ###
            shutil.copy('Receptores/{}.pdb'.format(Prot) , '{}/{}_{}'.format(Prot,Prot,lig))
            shutil.copy('Ligandos/{}.pdb'.format(lig) , '{}/{}_{}'.format(Prot,Prot,lig))
            shutil.copy('Bias/{}.bpf'.format(Prot) , '{}/{}_{}/ligand_bias.bpf'.format(Prot,Prot,lig))

            ### Paso a la Carpeta de Interes ###
            os.chdir('{}/{}_{}'.format(Prot,Prot,lig))
        
            ### Pasos Docking ###
            
            ### 1° Genero Ligando pdbqt Protonado ###
            os.system('pythonsh $MGLUTIL/prepare_ligand4.py -l {}.pdb -A hydrogens -o ligand.pdbqt'.format(lig))
            
            ### 2° Genero receptor pdbqt Protonado ###
            os.system('pythonsh $MGLUTIL/prepare_pdb_split_alt_confs.py -r {}.pdb'.format(Prot))
            if path.exists(Prot+'_A.pdb'):
                os.system('pythonsh $MGLUTIL/prepare_receptor4.py -A hydrogens -r {}_A.pdb -o receptor.pdbqt -U deleteAltB'.format(Prot))
            else:
                os.system('pythonsh $MGLUTIL/prepare_receptor4.py -A hydrogens -r {}.pdb -o receptor.pdbqt -U deleteAltB'.format(Prot))
            
            ### 3° GRILLA ###
            if int(Corte[0]) < 35:
                os.system('pythonsh $MGLUTIL/prepare_gpf4.py  -l  ligand.pdbqt -p npts="{}" -p gridcenter="{}" -r receptor.pdbqt -o receptor.gpf'.format(grilla_chica,Grid[0]))
            else:
                os.system('pythonsh $MGLUTIL/prepare_gpf4.py  -l  ligand.pdbqt -p npts="{}" -p gridcenter="{}" -r receptor.pdbqt -o receptor.gpf'.format(grilla_grande,Grid[0]))
            
            ### 3° MAPS ###
            os.system('autogrid4 -p receptor.gpf -l receptor.glg')
            
            ### 4° Generar DPF (Configura las opciones del docking) ###
            os.system('pythonsh $MGLUTIL/prepare_dpf42.py  -p  ga_run=100  -l  ligand.pdbqt  -r receptor.pdbqt -o receptor.dpf') 

            ### 5° Bias ###
        
            os.system('pythonsh  $MGLUTIL/contrib/adbias/bias2pdb.py ligand_bias.bpf')

            os.system('pythonsh  $MGLUTIL/contrib/adbias/prepare_bias.py -b ligand_bias.bpf -g receptor.gpf -d receptor.dpf')

            ## Mod maps ##

            path_actual = os.getcwd()

            Generar_Mapas_Bias(path_actual)

            ### Vuelvo al inicio ###
            os.chdir(path_complete)


##### Creo Slurm #####

salida = open('Citocromos_Docking_GPU.slurm' , 'w')
salida.write('#!/bin/bash\n')
salida.write('#SBATCH  --job-name=Dock_GPU         # Job name\n')
salida.write('#SBATCH --nodes=1                   # Run all processes on a single node\n')
salida.write('#SBATCH --ntasks=1                  # Run a single task\n')
salida.write('#SBATCH --cpus-per-task=1           # Number of CPU cores per task\n')
salida.write('#SBATCH -p gpu\n')
salida.write('#SBATCH --gres=gpu:1\n')
salida.write('#SBATCH --exclude=nodo13,nodo12\n')
salida.write('##SBATCH  --nodelist=nodoX\n')
salida.write('##SBATCH -e slurm-%j.err         ##Log del error\n')
salida.write('##SBATCH -o slurm-%j.out         ##Log salida normal\n')
salida.write('\n')

for Prot in Lista_Proteinas:
    Prot = Prot.split('.')[0]
    for lig in Lista_Ligandos:
        lig = lig.split('.')[0]
        salida.write('cd {}/{}_{}\n'.format(Prot,Prot,lig))
        salida.write('time srun /grupos/programs/autodock/AutoDock-GPU/bin/autodock_gpu_128wi  -L ligand.pdbqt -M receptor.biased_maps.fld --output-cluster-poses 1 -R ligand.pdbqt --nrun 100 --rmstol 3 -N docking_gpu\n')
        salida.write('cd ..\n')
        salida.write('cd ..\n')

salida.close

