import os
import pandas as pd
import shutil
from os import path


def Generar_Mapas_Bias(ruta_carpeta):
    
    file_name = 'receptor.maps.fld'
    file_bias = []
    file_normal = []

    # Lista los archivos en la ruta especificada y subdirectorios
    for carpeta_actual, subcarpetas, archivos in os.walk(ruta_carpeta):
        for archivo in archivos:
            if archivo.endswith('.biased.dpf'):
                ruta_archivo = os.path.join(carpeta_actual, archivo)
                input = open(ruta_archivo , 'r').readlines()
                for lines in input:
                    if 'biased' in lines:
                        file_bias.append(lines.split(' ')[1].rstrip())


    for str in file_bias:
        new_name = str
        new_name = new_name.replace(".biased", "")
        file_normal.append(new_name)


if __name__ == '__main__':

    ## Mod maps ##

    path_actual = os.getcwd()
    Generar_Mapas_Bias(path_actual)
