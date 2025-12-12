#!/bin/bash

# Carpeta que contiene los archivos .txt
carpeta="/home/your/file/coordenadas_genesTransc"

# Ir a la carpeta
cd "$carpeta" || exit

# Procesar cada archivo .txt en la carpeta
for archivo in *.txt; do
    # Definir el nombre del archivo de salida
    archivo_salida="${archivo%.txt}_posGen.txt"
    
    # Extraer las columnas 2, 3 y 5 y guardarlas en el archivo de salida
    awk -F'\t' '{print $2 "\t" $3 "\t" $5}' "$archivo" > "$archivo_salida"
done
