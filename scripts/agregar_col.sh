#!/bin/bash

# Verifica si el nombre del archivo se ha pasado como argumento
if [ $# -eq 0 ]; then
    echo "Por favor, proporciona el nombre del archivo."
    exit 1
fi

# Nombre del archivo
archivo=$1

# Crea un archivo temporal
archivo_temp=$(mktemp)

# Lee el archivo original y agrega la nueva columna
awk -F'\t' -v OFS='\t' '{print $0, "SC05_15"}' "$archivo" > "$archivo_temp"

# Reemplaza el archivo original con el archivo temporal
mv "$archivo_temp" "$archivo"

echo "Columna agregada exitosamente a $archivo"
