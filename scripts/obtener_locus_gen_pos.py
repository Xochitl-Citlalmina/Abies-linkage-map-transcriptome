import pandas as pd
import glob
import os

# Carpeta donde están los archivos (ajusta la ruta)
carpeta = "/home/nancy/Transcriptomes/BLAST/coordenadas_genesTransc"  # Cambia esto por la ruta real

# Lista para almacenar los datos procesados
datos_finales = []

# Procesar cada archivo de texto en la carpeta
for archivo in glob.glob(os.path.join(carpeta, "*.txt")):
    print(f"Procesando: {archivo}")
    
    # Leer el archivo con tabulaciones como separador
    df = pd.read_csv(archivo, sep="\t", header=None)
    
    # Seleccionar las columnas necesarias y calcular el promedio de las posiciones
    df_filtrado = df[[0, 1, 2, 4, 5]].copy()
    df_filtrado.columns = ["Locus", "Pos1", "Pos2", "Gen", "Grupo"]
    
    # Calcular el promedio de las posiciones
    df_filtrado["Pos_Avg"] = df_filtrado[["Pos1", "Pos2"]].mean(axis=1)
    
    # Eliminar columnas innecesarias
    df_filtrado = df_filtrado[["Locus", "Pos_Avg", "Gen", "Grupo"]]
    
    # Agregar los datos procesados a la lista
    datos_finales.append(df_filtrado)

# Concatenar todos los datos en un solo DataFrame
df_final = pd.concat(datos_finales, ignore_index=True)

# Eliminar duplicados manteniendo solo la primera aparición
df_final = df_final.drop_duplicates(subset=["Locus", "Gen", "Grupo"])

# Guardar en un archivo de salida
output_file = os.path.join(carpeta, "resultado_final.txt")
df_final.to_csv(output_file, sep="\t", index=False)

print(f"Proceso terminado. Archivo guardado en: {output_file}")
