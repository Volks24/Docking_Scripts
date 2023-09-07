import os
import shutil
from os import path
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from statistics import mean


def Histograma_Energia(Num_in_clus , Mean_Energy, caso):
    sns.set(style="darkgrid")
    Datos = []
    for j in range(0,len(Num_in_clus)):
        for k in range(0,Num_in_clus[j]):
            Datos.append(Mean_Energy[j])
    Fig = sns.histplot(data=Datos, y = Datos , kde=False , bins = len(Num_in_clus))
    Fig.set(xlabel='Number in Cluster', ylabel='Binding Energy' , title='CLUSTERING HISTOGRAM {} '.format(caso))
    plt.savefig(path_complete+'/Graficos/Histograma_'+caso+'.jpg')
    #plt.show()
    plt.close()
    return()

def Z_energy_vs_Z_Cluster(Ze,Zc,caso):
    plt.scatter(Ze, Zc, marker='x')
    #plt.axline((0, 0.5), slope=0.25, color="black", linestyle=(0, (5, 5)))
    plt.axvline(color="grey")
    plt.xlabel("Z Energy")
    plt.ylabel("Z Cluster")
    plt.title(caso)
    plt.savefig(path_complete+'/Graficos/Z_'+caso+'.jpg')
    #plt.show()
    plt.close()
    return()


def Armado_DataFrame(resultados,DF_Resumen_Data,DF_salida):
    for j in range(0,len(resultados)):
        if 'CLUSTERING HISTOGRAM' in resultados[j]:
            flag = 0
            Mean_Energy = []
            Num_in_clus = []
            while flag == 0:
                if '____' in resultados[j+9][0:4]:
                    Histograma_Energia(Num_in_clus , Mean_Energy , '{}_{}'.format(Prot,lig))
                    os.chdir(path_complete)
                    flag = 1
                else :
                    datos = resultados[j+9].split('|')
                    Mean_Energy.append(float(datos[3]))
                    Num_in_clus.append(int(datos[4]))
                    DF_salida.loc[len(DF_salida.index)] = [lig,int(number_of_atom),float(datos[3]),(float(datos[3])/number_of_atom),0,int(datos[4]),0]
                    j = j + 1
            DF_Resumen_Data.loc[len(DF_Resumen_Data)] = [Prot,lig,float(Mean_Energy[0]),len(Num_in_clus)]
    return(DF_Resumen_Data,DF_salida)


#### Preparo Archivos a Usar ###
if __name__ == '__main__':

    path_complete = os.getcwd()
    Lista_Ligandos = os.listdir(path_complete+'/Ligandos/')
    
    Lista_Proteinas = [name for name in os.listdir(path_complete) if os.path.isdir(os.path.join(path_complete, name))]
    
    
    #### Carpetas ####
    
    Carpetas = ['DLG' , 'Graficos' ,'DataFrames' , 'Temp' , 'Ligandos']
    
    for carpeta in Carpetas:
        if not os.path.exists(carpeta):
            os.makedirs(carpeta)
        else:
            Lista_Proteinas.remove(carpeta)


    DF_Resumen_Data = pd.DataFrame(columns=['Protein','Ligand','Minimun Energy' , 'Number of Pose'])

    for Prot in Lista_Proteinas:
        Prot = Prot.split('.')[0]
        resume_data = open(path_complete+'/Temp/'+str(Prot),'w')

        ###### 'Mean Energy Norm' = 'Mean Energy' / '# Atoms'
        DF_salida = pd.DataFrame(columns=['Ligand','# Atoms','Mean Energy', 'Mean Energy Norm','Z score energy','Number in Cluster','Z score cluster'])
        for lig in Lista_Ligandos:
            ### Obtengo # de atomos por ligando ###
            number_of_atom = 0
            lig = lig.split('.')[0]
            ligando = open('Ligandos/{}.pdb'.format(lig,'r')).readlines()
            for lines in ligando:
                if ('HETATM' in lines) or ('ATOM'in lines):
                    number_of_atom = number_of_atom + 1
            print('# atomos de {} : {}'.format(lig,len(ligando)))
            ### Caso ###
            print('Receptor {} Ligando {} \n'.format(Prot,lig))
            ### Paso a la Carpeta de Interes ###
            os.chdir('{}/{}_{}'.format(Prot,Prot,lig))
            shutil.copy('ligand_dock.dlg'.format(Prot,lig) , path_complete+'/DLG/{}_{}.dlg'.format(Prot,lig))
            resultados = open('ligand_dock.dlg','r').readlines()
            DF_Resumen_Data , DF_salida = Armado_DataFrame(resultados,DF_Resumen_Data,DF_salida)
        
                   
        #### calculos , busco la media por proteina , no por ligando ####
        # Z energia
        Media = DF_salida['Mean Energy Norm'].mean()
        Desvio = DF_salida['Mean Energy Norm'].std()
        DF_salida['Z score energy'] = (DF_salida['Mean Energy Norm'] - Media) / Desvio
        # Z Cluster
        Media_Cluster = DF_salida['Number in Cluster'].mean()
        Desvio_Cluster = DF_salida['Number in Cluster'].std()
        DF_salida['Z score cluster'] = (DF_salida['Number in Cluster'] - Media_Cluster) / Desvio_Cluster
        
        # Grafico Z vs Z
        Z_Energy = DF_salida['Z score energy'].tolist()
        Z_Cluster = DF_salida['Z score cluster'].tolist()
        Z_energy_vs_Z_Cluster(Z_Energy,Z_Cluster,'{}'.format(Prot))
        
        ### Redondeo para salida ###
        DF_salida = DF_salida.round(3)
        DF_salida.to_csv('DataFrames/{}.csv'.format(Prot))
        
        ### Armo subset con rakings por orden ###
        Sub_Set = DF_salida[['Ligand','Z score energy','Z score cluster']]
        Sub_Set = Sub_Set.sort_values(by=['Z score energy'], ascending=True)
        Sub_Set = Sub_Set.round(3)
        Sub_Set.to_csv('DataFrames/Ranking_{}.csv'.format(Prot))
       
        
        ### Archivo temporal con datos de interes ###
        resume_data.write('Media Energia: {}\n'.format(Media))
        resume_data.write('Desvio Energia: {}\n'.format(Desvio))
        resume_data.write('Media Cluster: {}\n'.format(Media_Cluster))
        resume_data.write('Desvio Cluster: {}\n'.format(Desvio_Cluster))
        resume_data.close()
    print(DF_Resumen_Data)
    #### Archivo Resumen Info ####
    DF_Resumen_Data.to_csv('DataFrames/{}'.format('Resumen.csv'))