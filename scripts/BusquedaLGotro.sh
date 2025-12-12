#!/bin/bash

# Nombres de los archivos
archivo1="locusLG12.txt"
archivo2="e-30.blastn.txt"
output="resultadoslg12.txt"

# Limpiar el archivo de resultados si ya existe
> $output

# Verificación de existencia de los archivos
if [ ! -f "$archivo1" ]; then
  echo "El archivo $archivo1 no existe."
  exit 1
fi

if [ ! -f "$archivo2" ]; then
  echo "El archivo $archivo2 no existe."
  exit 1
fi

# Leer el archivo1 línea por línea
while IFS= read -r id; do
  echo "Buscando ID: '$id'" # Mensaje de depuración
  # Buscar el ID en el archivo2 y extraer la primera y segunda columnas con coincidencia exacta
  awk -v id="$id" '$1 == id {print $1, $2}' "$archivo2" >> "$output"
done < "$archivo1"

# Confirmación de que el resultado se ha guardado
if [ -s "$output" ]; then
  echo "El resultado ha sido guardado en $output"
else
  echo "No se encontraron coincidencias."
fi
