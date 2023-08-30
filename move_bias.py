## Python
## Mueve los archivos para el bias

import os
import shutil
            
try:
    os.mkdir('Bias')
except:
    pass

Archivos = (os.listdir())

bias_files = []

for arc in Archivos:
    if 'map' in arc:
        bias_files.append(arc)
    if 'gpf' in arc:
        bias_files.append(arc)
    if 'pdbqt' in arc:
        bias_files.append(arc)
    if 'dpf' in arc:
        bias_files.append(arc)
            
## Armo los pdb para el bias ##
for files in bias_files:
	shutil.copy(files, "Bias/"+files)	
os.chdir('Bias')
os.system('obabel -ipdbqt ligand.pdbqt -opdb -Oligand_H.pdb')
os.system('obabel -ipdbqt receptor.pdbqt -opdb -Oreceptor_H.pdb')