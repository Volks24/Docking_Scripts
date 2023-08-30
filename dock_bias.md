###### Bias #####

##1 Crear carpeta con los Archivos necesarios del docking normal para aplicar el move_bias
python move_bias.py

##2 Marcar los aminoacido de interes en el receptor
### Sitio de Interaccion ###
## -r pos1,pos2,etc

pythonsh $MGLUTIL/contrib/adbias/ideal_interaction_sites.py -i receptor_H.pdb -c A -r 244,248,252,101,297

## Se crea ideal_interactions.pdb ##

## 3 Crear BPF ##
## Buscar los puntos de interes (3/4) y crear *.bpf

python make_bpf -p 3,5,6

### Crear Archivo bpf y conertir ###

pythonsh  $MGLUTIL/contrib/adbias/bias2pdb.py ligand_bias.bpf

### Modificar archivos dock para bias ###

pythonsh  $MGLUTIL/contrib/adbias/prepare_bias.py -b ligand_bias.bpf -g receptor.gpf -d receptor.dpf

### Autodock ###

## normal
autodock4 -p ligand_dock.biased.dpf -l ligand_dock.dlg 
## GPU
# modificar archvivo maps
python bias_gpu_maps.py

# Correr GPU

adt_gpu -L ligand_random.pdbqt -M receptor.maps.fld --output-cluster-poses 1 -R ligand.pdbqt --nrun 100
