#!/bin/bash

# Directorio que contiene los archivos
directory="/home/to/files/archivos_usar"

# Nombre del archivo de salida
output_file="genesorder_concatenado.txt"

# Verificar si el archivo de salida ya existe y eliminarlo si es así
[ -f $output_file ] && rm $output_file

# Obtener la lista de archivos en el directorio
files=("$directory"/*)

# Leer y concatenar archivos
for file in "${files[@]}"; do
    if [ "$file" == "${files[0]}" ]; then
        # Para el primer archivo, copiar todo incluyendo la primera línea
        cat "$file" >> "$output_file"
    else
        # Para los demás archivos, omitir la primera línea
        tail -n +2 "$file" >> "$output_file"
    fi
done
