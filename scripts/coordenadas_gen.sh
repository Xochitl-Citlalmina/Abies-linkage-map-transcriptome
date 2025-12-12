#!/bin/bash

# Nombres de los archivos
archivo1="order12_10PS.mapped"
archivo2="resultadoslg12.txt"
archivo_salida="coordenadas_geneslg12.txt"

# Crear el archivo de salida sin encabezado
> "$archivo_salida"

# Leer el archivo2 y construir un diccionario
declare -A diccionario

while IFS=$'\t' read -r col1_2 col2_2 col3_2 col4_2; do
    diccionario["$col1_2"]="$col1_2\t$col2_2\t$col3_2\t$col4_2"
done < "$archivo2"

# Leer el archivo1 y buscar coincidencias en el diccionario
while IFS=$'\t' read -r col1 col2 col3 col4; do
    if [[ -n "${diccionario[$col1]}" ]]; then
        # Escribir las columnas 1, 3, y 4 del archivo1 y la coincidencia en el archivo de salida
        echo -e "$col1\t$col3\t$col4\t${diccionario[$col1]}" >> "$archivo_salida"
    fi
done < "$archivo1"
