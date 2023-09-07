import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import zscore
import seaborn as sns
import dataframe_image as dfi


def scatter_energia_pob_all(Data_Frame , col_x , col_y):
    
    plt.figure(figsize=(14, 10))
    plot = sns.scatterplot(data=Data_Frame, x=col_x, y=col_y, palette='dark' ,s=70)

    

    plot.set_xlabel("Binding Energy (kcal/mol)")
    plot.set_ylabel("% Population")
    plot.set_title("All")
    
    #plot.set_facecolor("grey")  # Cambia el color del fondo
    plt.grid(True)
   
   
    plt.tight_layout()
    # Mostrar el gráfico
    plt.savefig('Plots/Docking_All.jpg', dpi=300)
    #plt.show()
    plt.close



def scatter_energia_pob(Data_Frame , col_x , col_y , TP ,Grupo):
    
    
    plt.figure(figsize=(14, 10))
    plot = sns.scatterplot(data=Data_Frame, x=col_x, y=col_y, hue='Grupo', palette='dark' ,s=70)

    # Resaltar un dato específico con un formato de punto diferente
    for lig in TP:
        dato_especifico = Data_Frame[Data_Frame['ligand'] == lig]  # Supongamos que queremos resaltar el punto con x=3
        sns.scatterplot(data=dato_especifico, x=col_x, y=col_y, color='red', s=200, marker='X')



    plot.set_xlabel("Binding Energy (kcal/mol)")
    plot.set_ylabel("% Population")
    plot.set_title("Grupo {}".format(Grupo))
    
    #plot.set_facecolor("grey")  # Cambia el color del fondo
    plt.grid(True)
   
    #Mostrar la leyenda y asignar el objeto de leyenda a una variable
    leyenda = plt.legend()
    # Cambiar el tamaño de la leyenda
    leyenda.set_bbox_to_anchor((1, 1))  # Ajusta la posición de la leyenda
    leyenda.set_title("Ligando")
    leyenda.get_title().set_fontsize('10')  # Cambia el tamaño del título de la leyenda
    for texto in leyenda.get_texts():
        texto.set_fontsize('8')  # Cambia el tamaño de los textos en la leyenda

    plt.tight_layout()
    # Mostrar el gráfico
    plt.savefig('Plots/Docking_Grupo_{}.jpg'.format(Grupo), dpi=300)
    #plt.show()
    plt.close


def Z_energia_pob(Data_Frame , col_x , col_y , TP ,Grupo):
    plt.figure(figsize=(14, 10))
    plot = sns.scatterplot(data=Data_Frame, x=col_x, y=col_y, hue='Grupo', palette='dark' ,s=70)

    # Resaltar un dato específico con un formato de punto diferente
    for lig in TP:
        dato_especifico = Data_Frame[Data_Frame['ligand'] == lig]  # Supongamos que queremos resaltar el punto con x=3
        sns.scatterplot(data=dato_especifico, x=col_x, y=col_y, color='red', s=200, marker='X')



    plot.set_xlabel("Z Energy")
    plot.set_ylabel("Z Population")
    plot.set_title("Grupo {} (ligandos: {})".format(Grupo,TP))
    
    # Ajustar tamaño de los números en los ejes
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    #plot.set_facecolor("grey")  # Cambia el color del fondo
    plt.grid(True)
   
   # Mostrar la leyenda y asignar el objeto de leyenda a una variable
    leyenda = plt.legend()
    # Cambiar el tamaño de la leyenda
    leyenda.set_bbox_to_anchor((1.05, 1))  # Ajusta la posición de la leyenda
    leyenda.set_title("Ligando")
    leyenda.get_title().set_fontsize('10')  # Cambia el tamaño del título de la leyenda
    for texto in leyenda.get_texts():
        texto.set_fontsize('8')  # Cambia el tamaño de los textos en la leyenda

    plt.tight_layout()
    # Mostrar el gráfico
    plt.savefig('Plots/Z_Grupo_{}.jpg'.format(Grupo), dpi=300)
    #plt.show()
    plt.close


def distribution_plot_energy(DF,Grupo):
    plt.figure(figsize=(8, 6))
    plot = sns.displot(DF, x="mean_energy" ,kind="kde")

    # Personalizar nombres de ejes
    plot.set_axis_labels("Binding Energy (kcal/mol)", "Frecuency")
    plot.set(title="Distribucion energia")
    
    plt.savefig('Plots/Dist_Grupo_{}.jpg'.format(Grupo), dpi=300)
    plt.close()


def scatter_energias(Data_Frame , col_x , col_y , TP ,Grupo):
    plt.figure(figsize=(14, 10))
    plot = sns.scatterplot(data=Data_Frame, x=col_x, y=col_y, hue='Grupo', palette='dark' ,s=70)

    # Resaltar un dato específico con un formato de punto diferente
    for lig in TP:
        dato_especifico = Data_Frame[Data_Frame['ligand'] == lig]  # Supongamos que queremos resaltar el punto con x=3
        sns.scatterplot(data=dato_especifico, x=col_x, y=col_y, color='red', s=200, marker='X')



    plot.set_xlabel("Binding Energy (kcal/mol)")
    plot.set_ylabel("Population Energy")
    plot.set_title("Grupo {}".format(Grupo))
    
    #plot.set_facecolor("grey")  # Cambia el color del fondo
    plt.grid(True)
   
    #Mostrar la leyenda y asignar el objeto de leyenda a una variable
    leyenda = plt.legend()
    # Cambiar el tamaño de la leyenda
    leyenda.set_bbox_to_anchor((1, 1))  # Ajusta la posición de la leyenda
    leyenda.set_title("Ligando")
    leyenda.get_title().set_fontsize('10')  # Cambia el tamaño del título de la leyenda
    for texto in leyenda.get_texts():
        texto.set_fontsize('8')  # Cambia el tamaño de los textos en la leyenda

    plt.tight_layout()
    # Mostrar el gráfico
    plt.savefig('Plots/Energias_Grupo_{}.jpg'.format(Grupo), dpi=300)
    #plt.show()
    plt.close


if __name__ == '__main__':

    Stats_DF = pd.DataFrame(columns=['Grupo' , 'Media_Energia' , 'Media_Cluster' , 'n'])

    Proteinas = ['1Q5D','1AKD','1PHA','2NZ5' , '5GWE']
    #Proteinas = ['1AKD' , '1IZO' , '1PHA' , '1Q5D' , '2CI0' , '2NZ5' , '3CV9' , '3DBZ' , '3G5H', '4DNJ' , '5GWE' , '5IKI', '5LI7' , '5OMU' , '6T0J' ]
    Proteinas_Grupo = {'1AKD':'3' , '1IZO':'1.1' , '1PHA':'5.1' , '1Q5D':'2' , '2CI0':'9' , '2NZ5':'8' , '3CV9':'4.2' , '3DBZ':'1.2' , '3G5H':'7', '4DNJ':'10.2' , '5GWE':'10.3' , '5IKI':'4.1', '5LI7':'5.2' , '5OMU':'10.1' , '6T0J':'1.3'}
    Lig_Original = {'1AKD':'CAM' , '1IZO':'VGJ' , '1PHA':'PFZ' , '1Q5D':'CHEBI_16089' , '2CI0':'1CM' , '2NZ5':'NQ' , '3CV9':'VDX' , '3DBZ':'MLI' , '3G5H':'YTT', '4DNJ':'FIV' , '5GWE':'88L' , '5IKI':'CHEMBL_112570', '5LI7':'6XD' , '5OMU':'V55' , '6T0J':'DXJ' }
    Lig_x_grupo = {'88L': '10.3', 'CHEBI_83951': '1.3', '1CM': '9', 'EPB': '2', 'ANN': '10.2', 'VGJ': '1.1', 'YTT': '7', 'CHEMBL_3804971': '6', 'M65': '5.2', 'CHEMBL_91': '5.1', 'ZMP': '1.4', 'V55': '10.1', 'CHEMBL_112570': '4.1', 'NQ': '8', 'JZ3': '10.1', '6XD': '5.2', 'CHEMBL_2048331': '4.1', 'CHEBI_16089': '2', 'FW6': '10.3', 'DXJ': '1.3', 'CAM': '3', 'CII': '9', 'FLV': '8', 'FIV': '10.2', 'CAH': '3', 'CHEBI_30807': '1.1', 'VDX': '4.2', 'CHEMBL_3805159': '6', '9LF': '2', 'MLI': '1.2', '7ZU': '4.2', 'MLA': '1.2', 'ZMO': '1.4', 'FJQ': '5.1', 'CHEMBL_4467901': '7', 'TWO': '10.2', '226': '8', 'PAM': '1.1', 'PFZ': '5.1', 'A9H': '4.1', 'GWM': '10.3', 'RWZ': '1.3', '3DM': '10.1'}
    Grupos = {'10.3': ['88L', 'FW6', 'GWM'], '1.3': ['CHEBI_83951', 'DXJ', 'RWZ'], '9': ['1CM', 'CII'], '2': ['EPB', 'CHEBI_16089', '9LF'], '10.2': ['ANN', 'FIV', 'TWO'], '1.1': ['VGJ', 'CHEBI_30807', 'PAM'], '7': ['YTT', 'CHEMBL_4467901'], '6': ['CHEMBL_3804971', 'CHEMBL_3805159'], '5.2': ['M65', '6XD'], '5.1': ['CHEMBL_91', 'FJQ', 'PFZ'], '1.4': ['ZMP', 'ZMO'], '10.1': ['V55', 'JZ3', '3DM'], '4.1': ['CHEMBL_112570', 'CHEMBL_2048331', 'A9H'], '8': ['NQ', 'FLV', '226'], '3': ['CAM', 'CAH'], '4.2': ['VDX', '7ZU'], '1.2': ['MLI', 'MLA']}

    columnas_ranks = list(Grupos.keys())
    columnas_ranks.insert(0, 'Grupo')
    
    Ranks_DF = pd.DataFrame(columns=columnas_ranks , index=Proteinas)

    DF = pd.read_csv('Resume_data.csv' , index_col='Unnamed: 0')
    ### '5LIB':'6XD' revisar

    DF_All = DF.copy()

    DF_All = DF_All.drop(DF_All[DF_All['mean_energy'] > 0].index)


    scatter_energia_pob_all(DF_All , 'mean_energy' , 'cluster_pob')

    ### Agrego Grupos de cada ligando ####
    DF_All['Grupo'] = ''
    for j in range(0,DF_All.shape[0]):
        DF_All.iloc[j,9]= Lig_x_grupo[DF_All.iloc[j,1]]

    ### Calculo de Poblacion como Energia ###

    DF_All['Energia_P'] = 100 - DF_All['cluster_pob']
    DF_All['Energia_P'] = DF_All['Energia_P'].replace(0, 0.0001)
    DF_All['Energia_P'] = np.log(DF_All['Energia_P'])
    DF_All['Energia_P'] = -0.6*(DF_All['Energia_P'])

    for k in range(0,len(Proteinas)):
     
    
        prot = Proteinas[k]

        sub_Set = DF_All.query('proteina == @prot').copy()
        Grupo = (Proteinas_Grupo[prot])

        ## Stats ##
        ### Filtro Datos ###
        sub_Set = sub_Set.drop(sub_Set[sub_Set['mean_energy'] > 0].index)
        sub_Set = sub_Set.drop(sub_Set[sub_Set['cluster_pob'] < 10].index)

        n = sub_Set.shape[0]
        Media_energia_gral = sub_Set['mean_energy'].mean()
        Std_energia_gral = sub_Set['mean_energy'].std()
        Media_cluster_gral = sub_Set['cluster_pob'].mean()
        Std_cluster_gral = sub_Set['cluster_pob'].std()
        print('Grupo {} Media Energia {} Media Poblacion {} Casos {}'.format(Grupo, Media_energia_gral,Media_cluster_gral,n))
        Stats_DF.loc[(Stats_DF.shape[0])] = [Grupo, Media_energia_gral,Media_cluster_gral,n]

        sub_Set['Z_energy'] = zscore(sub_Set['mean_energy'])
        sub_Set['Z_cluster'] = zscore(sub_Set['cluster_pob'])

        Ligandos = list(sub_Set['ligand'].unique())
        for lig in Ligandos:
            sub_set_z = sub_Set.query('ligand == @lig')
            Pos = (sub_set_z.index[0])
            if (sub_set_z.shape[0]) == 1:
                Z_e = sub_set_z.iloc[0,5] / Std_energia_gral
                sub_Set.loc[Pos,'Z_energy'] = Z_e
                Z_p = sub_set_z.iloc[0,6] / Std_cluster_gral
                sub_Set.loc[Pos,'Z_cluster'] = Z_p
            
               
        ### Distribucion Energia
        distribution_plot_energy(sub_Set,Grupo)
        
        ### Scatter  plots Dock
        scatter_energia_pob(sub_Set ,'mean_energy' , 'cluster_pob' ,Grupos[Grupo],Grupo)
        sub_Set = sub_Set.query('mean_energy < 0')
        scatter_energia_pob(sub_Set ,'mean_energy' , 'cluster_pob' ,Grupos[Grupo],Grupo)
    
    
    
        ### Scatter  plots Z
        Z_energia_pob(sub_Set ,'Z_energy' , 'Z_cluster' ,Grupos[Grupo],Grupo)
        sub_Set.to_csv('DF/Grupo_{}.csv'.format(Grupo))

        ### Scatter  plots Energias
        scatter_energias(sub_Set ,'mean_energy' , 'Energia_P' ,Grupos[Grupo],Grupo)


        #### Ranks ####

        sub_Set_sorted = sub_Set.sort_values(by=['Z_energy', 'Z_cluster'])

        print(sub_Set_sorted)
        casto_TP  =  Lig_Original[prot]
        Corte = (sub_Set_sorted.query('ligand == @casto_TP'))
        print(Corte)
        Umbral = (Corte.iloc[0,11])
        Umbral = (int(Umbral))

        Plot_res = (sub_Set_sorted.query('Z_energy <= @Umbral'))

        conteo = dict(Plot_res['Grupo'].value_counts())
        

        Ranks_DF.loc[prot,'Grupo'] = Proteinas_Grupo[prot]
        for llaves in conteo:
            Ranks_DF.loc[prot,llaves] = conteo[llaves]

    columnas_interes = [ "10.3", "1.3", "9", "2", "10.2", "1.1", "7", "6", "5.2", "5.1", "1.4", "10.1", "4.1", "8", "3", "4.2", "1.2"]
    Ranks_DF = Ranks_DF.fillna(0)
    
    Ranks_DF['Total'] = 0
    for pro in Proteinas:
        Ranks_DF.loc[pro,'Total'] = Ranks_DF.loc[pro, columnas_interes].sum()
           
    print(Ranks_DF)
    Ranks_DF.to_csv('DF/Ranks_DF.csv')
    Ranks_DF = Ranks_DF.style.background_gradient()
    dfi.export(Ranks_DF, 'Plots/Ranks.png')
    
    
    
    
    Stats_DF['n'] = Stats_DF['n'].astype(int)
    Stats_DF['Grupo'] = Stats_DF['Grupo'].astype(str)
    
    Stats_DF.to_csv('DF/Resumen_grupos.csv')
    Stats_DF = Stats_DF.style.background_gradient()
    dfi.export(Stats_DF, 'Plots/Resume_data_grupos.png')

    