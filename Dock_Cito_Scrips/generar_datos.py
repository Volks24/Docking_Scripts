import os
from os import path
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from statistics import mean



def buscar_archivos_por_nombre(carpeta, nombre_deseado):
    archivos_encontrados = []

    # Obtener la lista de archivos en la carpeta
    archivos_en_carpeta = os.listdir(carpeta)

    # Filtrar los archivos por nombre
    for archivo in archivos_en_carpeta:
        if archivo.startswith(nombre_deseado):
            archivos_encontrados.append(archivo)

    return (archivos_encontrados)

def Armado_DataFrame(input_file,case,Number_Atoms):
    enegy_df = pd.DataFrame(columns=['proteina','ligand','rank', 'lowest energy', 'run', 'mean_energy', 'cluster_pob' , 'Atom', 'normalizer_energy'])
    for j in range(0,len(input_file)):
        if 'CLUSTERING HISTOGRAM' in input_file[j]:       
            flag = True
            pos = j + 10
            while flag == True:
                if input_file[pos][0:5] != '_____':
                    data = input_file[pos].split('|')
                    valores = [str(prot),str(case),int(data[0]),float(data[1]),int(data[2]),float(data[3]),int(data[4]), Number_Atoms,float(data[3])/Number_Atoms ]
                    temp_df = pd.DataFrame([valores], columns=enegy_df.columns)
                    enegy_df = pd.concat([enegy_df, temp_df], ignore_index=True)
                    pos = pos + 1
                else:
                    flag = False

    return(enegy_df)





if __name__ == '__main__':

    Proteinas = ['1AKD' , '1IZO' , '1PHA' , '1Q5D' , '2CI0' , '2NZ5' , '3CV9' , '3DBZ' , '3G5H', '4DNJ' , '5GWE' , '5IKI', '5LI7' , '5OMU' , '6T0J' ]

    DF_Resumen_Data = pd.DataFrame(columns=['proteina','ligand','rank', 'lowest energy', 'run', 'mean_energy', 'cluster_pob' , 'Atom', 'normalizer_energy'])
    
    for prot in Proteinas:    
        
        DF_salida = pd.DataFrame(columns=['Ligand','# Atoms','Mean_Energy', 'Mean Energy Norm','Z score energy','Number in Cluster','Z score cluster'])
        archivos = buscar_archivos_por_nombre('DLG_GPU', prot)
        
        for DLG_File in (archivos):

            ### Obtengo # de atomos por ligando ###
            number_of_atom = 0
            lig = DLG_File[5:].split('.')[0]
            ligando = open('Ligandos/{}.pdb'.format(lig,'r')).readlines()
            for lines in ligando:
                if ('HETATM' in lines) or ('ATOM'in lines):
                    number_of_atom = number_of_atom + 1
            print('# atomos de {} : {}'.format(lig,number_of_atom))

            resultados = open('DLG_GPU/{}'.format(DLG_File),'r').readlines()
            enegy_df = Armado_DataFrame(resultados,lig,number_of_atom)
            DF_Resumen_Data = pd.concat([DF_Resumen_Data, enegy_df], ignore_index=True)

            

    DF_Resumen_Data.to_csv('Resume_data.csv')

